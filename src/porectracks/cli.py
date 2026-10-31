"""Command-line entry points for track construction and inspection."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

from . import __version__


def genomic_region(value):
    match = re.fullmatch(r"([^:]+):(\d[\d,]*)-(\d[\d,]*)", value)
    if match is None:
        raise argparse.ArgumentTypeError("Use chromosome:start-end, in 0-based half-open bp")
    chrom, left, right = match.groups()
    start, end = int(left.replace(",", "")), int(right.replace(",", ""))
    if start >= end:
        raise argparse.ArgumentTypeError("Region start must be smaller than end")
    return chrom, start, end


def main(argv=None):
    parser = argparse.ArgumentParser(prog="porectracks", description="Pore-C contact distributions as genomic tracks")
    parser.add_argument("--version", action="version", version=f"PoreCTracks {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build", help="Convert a member TSV into four cumulative BigWigs")
    build.add_argument("--contacts", type=Path, required=True, help="Grouped set_id/chrom/start/end[/count] TSV, optionally .gz")
    build.add_argument("--chrom-sizes", type=Path, required=True)
    build.add_argument("--out", type=Path, required=True, help="New or empty output directory")
    build.add_argument("--bin-size", type=int, default=1024)
    build.add_argument("--max-order", type=int, default=10, help="Accepted source order cap (at least 5)")
    build.add_argument("--blacklist", type=Path, help="BED regions to mask from the output")
    build.add_argument("--regions", type=Path, help="BED observation territory; require full-bin coverage")
    build.add_argument("--support-total", type=float, help="Shared external denominator for a subset of a larger dataset")
    build.add_argument("--normalization", choices=("per-million", "raw"), default="per-million")
    plot = commands.add_parser("plot", help="Plot four Pore-C contact tracks")
    plot.add_argument("--tracks", type=Path, required=True, help="Completed build directory")
    plot.add_argument("--region", type=genomic_region, required=True)
    plot.add_argument("--out", type=Path, required=True, help="Figure prefix for PNG and PDF")
    plot.add_argument("--atac", type=Path, help="Add an ATAC annotation panel")
    plot.add_argument("--dpi", type=int, default=300, help="PNG resolution in dots per inch (default: 300)")
    stats = commands.add_parser("stats", help="Plot source-order composition and cis span distributions")
    stats.add_argument("--tracks", type=Path, required=True)
    stats.add_argument("--out", type=Path, required=True, help="Figure prefix for PNG and PDF")
    stats.add_argument("--dpi", type=int, default=300, help="PNG resolution in dots per inch (default: 300)")
    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            from .core import build_tracks
            result = build_tracks(args.contacts, args.chrom_sizes, args.out, args.bin_size,
                                  args.max_order, args.blacklist, args.regions,
                                  args.support_total, args.normalization)
            print(json.dumps({"status": result["status"], "output": str(args.out), "counts": result["counts"]}))
        elif args.command == "plot":
            from .plotting import plot_region
            outputs = plot_region(args.tracks, args.region, args.out, args.atac, dpi=args.dpi)
            print(json.dumps({"figures": [str(p) for p in outputs]}))
        else:
            from .plotting import plot_distributions
            outputs = plot_distributions(args.tracks, args.out, dpi=args.dpi)
            print(json.dumps({"figures": [str(p) for p in outputs]}))
    except (ValueError, OSError, RuntimeError, ArithmeticError) as exc:
        parser.exit(2, f"porectracks: {exc}\n")
    return 0
