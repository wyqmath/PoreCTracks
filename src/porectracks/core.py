from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import csv
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pyBigWig

from . import __version__
from .io import iter_sets, read_bed, read_sizes

THRESHOLDS = (2, 3, 4, 5)
CHANNELS = tuple(f"contact_{m}plus" for m in THRESHOLDS)


def _masks(sizes, bin_size, blacklist=None, regions=None):
    coverage, unused_regions = read_bed(regions, sizes) if regions else ({}, 0)
    blocked, unused_blacklist = read_bed(blacklist, sizes) if blacklist else ({}, 0)
    masks, summary = {}, {}
    for chrom, length in sizes.items():
        n = (length + bin_size - 1) // bin_size
        mask = np.ones(n, dtype=bool) if regions is None else np.zeros(n, dtype=bool)
        for start, end in coverage.get(chrom, []):
            first = (start + bin_size - 1) // bin_size
            stop = n if end == length else end // bin_size
            mask[first:stop] = True
        covered = int(mask.sum())
        for start, end in blocked.get(chrom, []):
            mask[start // bin_size:(end + bin_size - 1) // bin_size] = False
        masks[chrom] = mask
        summary[chrom] = {"total_bins": n, "covered_bins": covered,
                          "blacklisted_covered_bins": covered - int(mask.sum()),
                          "valid_bins": int(mask.sum())}
    if not any(m.any() for m in masks.values()):
        raise ValueError("No valid bins remain after regions and blacklist masking")
    return masks, {"chromosomes": summary, "ignored_region_intervals": unused_regions,
                   "ignored_blacklist_intervals": unused_blacklist}


def _write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def build_tracks(contacts, chrom_sizes, output, bin_size=1024, max_order=10,
                 blacklist=None, regions=None, support_total=None,
                 normalization="per-million"):
    """Build four BigWigs and diagnostics; output masks never alter source sets.

    Each accepted k-member set with support c contributes c/k to every member
    in each threshold channel m <= k. The denominator counts each set once.
    """
    if not isinstance(bin_size, int) or isinstance(bin_size, bool) or bin_size < 1:
        raise ValueError("bin_size must be a positive integer")
    if not isinstance(max_order, int) or isinstance(max_order, bool) or max_order < 5:
        raise ValueError("max_order must be an integer of at least 5")
    if normalization not in ("raw", "per-million"):
        raise ValueError("normalization must be raw or per-million")
    if support_total is not None and (not math.isfinite(support_total) or support_total <= 0):
        raise ValueError("support_total must be finite and positive")
    output = Path(output)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError(f"Output must be new or empty: {output}")
    sizes = read_sizes(chrom_sizes)
    masks, mask_report = _masks(sizes, bin_size, blacklist, regions)
    arrays = {chrom: np.zeros((4, len(mask)), dtype=np.float64) for chrom, mask in masks.items()}
    counts = dict.fromkeys(("input_sets", "input_support", "accepted_sets", "accepted_support",
                           "excluded_low_order_sets", "excluded_low_order_support",
                           "excluded_high_order_sets", "excluded_high_order_support",
                           "cis_support", "trans_support", "input_members", "collapsed_members"), 0)
    order_stats = {k: Counter() for k in range(2, max_order + 1)}
    spans = Counter()
    observed_orders = Counter()
    for _, members, support, n_members in iter_sets(contacts, sizes, bin_size):
        k = len(members)
        counts["input_sets"] += 1
        counts["input_support"] += support
        counts["input_members"] += n_members
        counts["collapsed_members"] += n_members - k
        observed_orders[k] += support
        if k < 2 or k > max_order:
            label = "low" if k < 2 else "high"
            counts[f"excluded_{label}_order_sets"] += 1
            counts[f"excluded_{label}_order_support"] += support
            continue
        counts["accepted_sets"] += 1
        counts["accepted_support"] += support
        cis = len({chrom for chrom, _ in members}) == 1
        kind = "cis_support" if cis else "trans_support"
        counts[kind] += support
        order_stats[k].update({"sets": 1, "support": support, kind: support})
        if cis:
            spans[(k, (members[-1][1] - members[0][1]) * bin_size)] += support
        for chrom, index in members:
            arrays[chrom][:min(k, 5) - 1, index] += support / k
    if counts["accepted_support"] == 0:
        raise ValueError("No sets have an accepted order after binning")
    denominator = counts["accepted_support"] if support_total is None else support_total
    if denominator < counts["accepted_support"]:
        raise ValueError("support_total must be at least the accepted input support")
    factor = 1e6 / denominator if normalization == "per-million" else 1.0
    expected = np.array([sum(order_stats[k]["support"] for k in range(m, max_order+1)) for m in THRESHOLDS])
    full_mass = sum((a.sum(axis=1) for a in arrays.values()), np.zeros(4))
    if not np.allclose(full_mass, expected, rtol=1e-10, atol=1e-8):
        raise ArithmeticError("Accumulated support does not match source support")
    exported_mass = sum((a[:, masks[c]].sum(axis=1) for c, a in arrays.items()), np.zeros(4))
    metadata = {
        "status": "COMPLETE", "version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "bin_size": bin_size, "max_order": max_order, "channels": list(CHANNELS),
        "coordinate_system": "0-based half-open", "chrom_sizes": sizes,
        "order_definition": "Distinct (chromosome, midpoint-bin) members per input set; adjacent bins remain distinct",
        "weighting": "For each threshold m<=k, a k-member set contributes count/k per member",
        "normalization": {"mode": normalization, "denominator": denominator, "factor": factor,
                          "denominator_source": "accepted-input-support" if support_total is None else "explicit-support-total",
                          "scope": "shared across all channels and chromosomes, before output masks"},
        "counts": counts, "input_order_support": dict(sorted(observed_orders.items())),
        "mask": {"policy": "Fully covered genomic bins; any blacklist overlap is masked; missing entries are unobserved",
                 "regions_supplied": regions is not None, "blacklist_supplied": blacklist is not None, **mask_report},
        "support": {"before_mask": dict(zip(CHANNELS, full_mass.tolist())),
                    "exported_before_scaling": dict(zip(CHANNELS, exported_mass.tolist()))},
        "statistics_scope": "All accepted input sets before the output mask; span statistics cover cis sets only",
        "inputs": {"contacts": Path(contacts).name, "chrom_sizes": Path(chrom_sizes).name,
                   "blacklist": Path(blacklist).name if blacklist else None,
                   "regions": Path(regions).name if regions else None},
        "dependencies": {name: importlib.metadata.version(name) for name in ("numpy", "pyBigWig", "matplotlib")},
        "outputs": [*(f"{name}.bw" for name in CHANNELS), "valid_regions.bed", "order_stats.tsv", "span_distribution.tsv", "metadata.json"],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        with (stage / "valid_regions.bed").open("w", encoding="utf-8") as bed:
            for chrom, mask in masks.items():
                transitions = np.diff(np.r_[False, mask, False].astype(np.int8))
                for first, stop in zip(np.flatnonzero(transitions == 1), np.flatnonzero(transitions == -1)):
                    bed.write(f"{chrom}\t{int(first)*bin_size}\t{min(int(stop)*bin_size,sizes[chrom])}\n")
        for channel, name in enumerate(CHANNELS):
            with pyBigWig.open(str(stage / (name + ".bw")), "w") as bw:
                bw.addHeader(list(sizes.items()))
                for chrom, mask in masks.items():
                    ids = np.flatnonzero(mask)
                    for offset in range(0, len(ids), 100000):
                        selected = ids[offset:offset+100000]
                        starts = selected * bin_size
                        values = arrays[chrom][channel, selected] * factor
                        if not np.all(np.isfinite(values)) or np.any(values > np.finfo(np.float32).max):
                            raise ValueError("Track values exceed the finite BigWig float32 range")
                        bw.addEntries([chrom] * len(selected), starts.tolist(),
                                      ends=np.minimum(starts + bin_size, sizes[chrom]).tolist(),
                                      values=values.tolist())
        fields = ("sets", "support", "cis_support", "trans_support")
        _write_tsv(stage / "order_stats.tsv", ("order", *fields),
                   ((k, *(order_stats[k][f] for f in fields)) for k in order_stats))
        _write_tsv(stage / "span_distribution.tsv", ("order", "span_bp", "support"),
                   ((k, span, support) for (k, span), support in sorted(spans.items())))
        (stage / "metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
        if output.exists():
            output.rmdir()
        os.replace(stage, output)
    except BaseException:
        shutil.rmtree(stage)
        raise
    return metadata
