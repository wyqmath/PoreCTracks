from __future__ import annotations

import csv
import gzip
from pathlib import Path


def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8-sig") if path.suffix == ".gz" else path.open(encoding="utf-8-sig")


def read_sizes(path):
    sizes = {}
    with open_text(path) as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(f"chrom.sizes line {line_no}: expected chromosome and length")
            chrom, raw = fields
            try:
                length = int(raw)
            except ValueError as exc:
                raise ValueError(f"chrom.sizes line {line_no}: length must be an integer") from exc
            if chrom in sizes or not 0 < length < 2**32:
                raise ValueError(f"chrom.sizes line {line_no}: duplicate chromosome or invalid length")
            sizes[chrom] = length
    if not sizes:
        raise ValueError("chrom.sizes is empty")
    return sizes


def read_bed(path, sizes):
    """Read and merge intervals. Reference annotation may list unused chromosomes."""
    intervals = {chrom: [] for chrom in sizes}
    ignored = 0
    with open_text(path) as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"BED line {line_no}: expected chromosome, start, end")
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"BED line {line_no}: coordinates must be integers") from exc
            if start < 0 or end <= start:
                raise ValueError(f"BED line {line_no}: invalid interval")
            if chrom not in sizes:
                ignored += 1
                continue
            if end > sizes[chrom]:
                raise ValueError(f"BED line {line_no}: interval exceeds chromosome length")
            intervals[chrom].append((start, end))
    merged = {}
    for chrom, rows in intervals.items():
        result = []
        for start, end in sorted(rows):
            if result and start <= result[-1][1]:
                result[-1] = (result[-1][0], max(end, result[-1][1]))
            else:
                result.append((start, end))
        merged[chrom] = result
    return merged, ignored


def iter_sets(path, sizes, bin_size):
    """Yield (id, sorted distinct bins, count, input_member_count); stream grouped TSV."""
    with open_text(path) as stream:
        lines = (line for line in stream if line.strip() and not line.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        required = {"set_id", "chrom", "start", "end"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("Contacts header requires set_id, chrom, start, end; count is optional")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Contacts header contains duplicate columns")
        seen = set()
        current = None
        members = set()
        count = None
        n_members = 0
        for row in reader:
            if None in row or any(row.get(key) is None for key in required):
                raise ValueError(f"Contacts row {reader.line_num}: inconsistent column count")
            set_id, chrom = row["set_id"].strip(), row["chrom"].strip()
            if not set_id or chrom not in sizes:
                raise ValueError(f"Contacts row {reader.line_num}: empty set_id or unknown chromosome")
            try:
                start, end = int(row["start"]), int(row["end"])
                support = int(row["count"]) if "count" in reader.fieldnames else 1
            except (ValueError, TypeError) as exc:
                raise ValueError(f"Contacts row {reader.line_num}: coordinates/count must be integers") from exc
            if not 0 <= start < end <= sizes[chrom] or support <= 0:
                raise ValueError(f"Contacts row {reader.line_num}: invalid interval or nonpositive count")
            if set_id != current:
                if set_id in seen:
                    raise ValueError(f"Set {set_id!r} reappeared; rows of each set_id must be contiguous")
                if current is not None:
                    yield current, tuple(sorted(members)), count, n_members
                seen.add(set_id)
                current, members, count, n_members = set_id, set(), support, 0
            elif support != count:
                raise ValueError(f"Set {set_id!r}: count must be identical on every member row")
            members.add((chrom, ((start + end) // 2) // bin_size))
            n_members += 1
        if current is not None:
            yield current, tuple(sorted(members)), count, n_members
