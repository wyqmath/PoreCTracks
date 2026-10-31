from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable

import matplotlib.ticker as mticker
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy as np
import pyBigWig


_CHANNELS = ("contact_2plus", "contact_3plus", "contact_4plus", "contact_5plus")
_COLORS = ("#326A9D", "#3E8C78", "#CB9139", "#9A5889")
Interval = tuple[int, int, float]


def _metadata(run_dir: Path) -> dict:
    with (run_dir / "metadata.json").open(encoding="utf-8") as handle:
        return json.load(handle)


def _save(fig: Figure, output: Path, dpi: int = 300) -> list[Path]:
    """Save both publication and preview formats from an output prefix."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    paths = [Path(str(output) + suffix) for suffix in (".png", ".pdf")]
    for path in paths:
        fig.savefig(path, dpi=dpi, facecolor="white")
    return paths


def _intervals(path: Path, chrom: str, start: int, end: int) -> list[Interval]:
    """Read overlapping full bin intervals; the plot viewport clips their edges."""
    with pyBigWig.open(str(path)) as bw:
        length = bw.chroms(chrom)
        if length is None or start >= length:
            return []
        values = bw.intervals(chrom, start, min(end, length)) or ()
        return [(left, right, float(value))
                for left, right, value in values if np.isfinite(value)]


def _segments(intervals: Iterable[Interval]) -> tuple[list[float], list[float]]:
    """Join adjacent interval edges, inserting NaNs only across coverage gaps."""
    xs: list[float] = []
    ys: list[float] = []
    previous_end: int | None = None
    for left, right, value in intervals:
        if previous_end is not None and left != previous_end:
            xs.append(float("nan"))
            ys.append(float("nan"))
        xs.extend((left, right))
        ys.extend((value, value))
        previous_end = right
    return xs, ys


def _atac_intervals(path: Path, chrom: str, intervals: list[Interval]) -> list[Interval]:
    """Exact covered-base ATAC mean over each full valid contact bin."""
    result: list[Interval] = []
    with pyBigWig.open(str(path)) as bw:
        length = bw.chroms(chrom)
        if length is None:
            return result
        for left, right, _ in intervals:
            if left >= length:
                continue
            value = bw.stats(chrom, left, min(right, length), type="mean", exact=True)[0]
            if value is not None and np.isfinite(value):
                if value < 0:
                    raise ValueError("ATAC values must be nonnegative for log1p display")
                result.append((left, right, float(np.log1p(value))))
    return result


def plot_region(
    run_dir: str | Path,
    region: tuple[str, int, int],
    output: str | Path,
    atac: str | Path | None = None,
    dpi: int = 300,
) -> list[Path]:
    """Plot cumulative contact tracks for a 0-based, half-open genomic region.

    Output is a prefix: ``.png`` and ``.pdf`` are appended (300 dpi by default).
    Missing BigWig
    intervals remain gaps. ATAC uses the exact covered-base mean over each
    full valid bin, independent of the requested display boundaries.
    """
    run_dir = Path(run_dir)
    chrom, start, end = region
    if start < 0 or end <= start:
        raise ValueError("Region requires 0 <= start < end (0-based, half-open)")
    meta = _metadata(run_dir)
    chrom_sizes = meta["chrom_sizes"]
    if chrom not in chrom_sizes:
        raise ValueError(f"Unknown chromosome: {chrom}")
    if end > int(chrom_sizes[chrom]):
        raise ValueError(f"Region end exceeds chromosome length ({chrom}: {chrom_sizes[chrom]})")
    tracks = [_intervals(run_dir / f"{name}.bw", chrom, start, end)
              for name in _CHANNELS]
    panels = 5 if atac is not None else 4
    fig = Figure(figsize=(10.5, 1.2 * panels + 0.9))
    FigureCanvasAgg(fig)
    axes = fig.subplots(panels, 1, sharex=True, squeeze=False)[:, 0]
    ymax = max((value for track in tracks for _, _, value in track), default=0.0)
    ymax = ymax * 1.08 if ymax > 0 else 1.0
    normalized = meta.get("normalization", {}).get("mode") == "per-million"
    bin_label = f"{int(meta['bin_size']):,}-bp bin"
    units = (f"Contact intensity (per million sets / {bin_label})" if normalized
             else f"Contact support / {bin_label}")
    for index, (ax, intervals) in enumerate(zip(axes[:4], tracks)):
        xs, ys = _segments(intervals)
        ax.plot(xs, ys, color=_COLORS[index], linewidth=1.2)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(f"{index + 2}+", fontsize=12, rotation=0, labelpad=18, va="center")
        if not intervals:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center", color="0.45")
    if atac is not None:
        atac_values = _atac_intervals(Path(atac), chrom, tracks[0])
        xs, ys = _segments(atac_values)
        axes[-1].plot(xs, ys, color="#656565", linewidth=1.2)
        axes[-1].set_ylabel("ATAC\nlog1p(mean)", fontsize=9)
        axes[-1].set_ylim(bottom=0, top=max((v for _, _, v in atac_values), default=0) * 1.08 or 1)
        if not atac_values:
            axes[-1].text(0.5, 0.5, "No data",
                          transform=axes[-1].transAxes, ha="center", va="center", color="0.45")
    for ax in axes:
        ax.set_xlim(start, end)
        ax.grid(axis="y", alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
        ax.xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.0f}"))
    axes[-1].set_xlabel("Position (bp)")
    axes[0].set_title(units, loc="left", fontsize=9, pad=7)
    fig.suptitle(f"{chrom}:{start:,}–{end:,}", fontsize=12, y=0.98)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.88, bottom=0.12, hspace=0.13)
    return _save(fig, Path(output), dpi=dpi)


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def plot_distributions(run_dir: str | Path, output: str | Path, dpi: int = 300) -> list[Path]:
    """Plot count-weighted source-order composition and exact cis-span CDFs.

    Composition uses source-set support, before pair distribution. Span is the
    full cis set's binned genomic extent, never a span of a derived subset.
    """
    run_dir = Path(run_dir)
    meta = _metadata(run_dir)
    stats = _rows(run_dir / "order_stats.tsv")
    spans = _rows(run_dir / "span_distribution.tsv")
    maximum = int(meta["max_order"])
    accepted = [row for row in stats if 2 <= int(row["order"]) <= maximum]
    orders = sorted({int(row["order"]) for row in accepted})
    support = {order: sum(float(row["support"]) for row in accepted
                          if int(row["order"]) == order) for order in orders}
    cis = {order: sum(float(row["cis_support"]) for row in accepted
                      if int(row["order"]) == order) for order in orders}
    trans = {order: sum(float(row["trans_support"]) for row in accepted
                        if int(row["order"]) == order) for order in orders}
    total = sum(support.values())
    trans_total = sum(trans.values())
    fig = Figure(figsize=(11, 5.2))
    FigureCanvasAgg(fig)
    left, right = fig.subplots(1, 2)
    if total > 0:
        cis_pct = [100 * cis[order] / total for order in orders]
        trans_pct = [100 * trans[order] / total for order in orders]
        left.bar(orders, cis_pct, width=0.75, label="cis (one chromosome)", color="#326A9D")
        left.bar(orders, trans_pct, bottom=cis_pct, width=0.75,
                 label="trans (multiple chromosomes)", color="#A9B8C5")
        left.legend(frameon=False, fontsize=8)
        left.set_ylim(0, max(1, max(a + b for a, b in zip(cis_pct, trans_pct)) * 1.2))
    else:
        left.text(0.5, 0.5, "No accepted source-set support", transform=left.transAxes,
                  ha="center", va="center", color="0.45")
        left.set_ylim(0, 100)
    left.set_title("Source-order composition", fontsize=11)
    left.set_xlabel("Full source-set order")
    left.set_ylabel("Share of accepted support (%)")
    if orders:
        if len(orders) <= 18:
            left.set_xticks(orders)
        else:
            left.xaxis.set_major_locator(mticker.MaxNLocator(nbins=10, integer=True))
    curves = 0
    for group, color in zip((2, 3, 4, 5), _COLORS):
        combined: dict[int, float] = {}
        for row in spans:
            order = int(row["order"])
            if not 2 <= order <= maximum or (order != group if group < 5 else order < 5):
                continue
            weight = float(row["support"])
            if weight > 0:
                span = int(row["span_bp"])
                combined[span] = combined.get(span, 0.0) + weight
        if not combined:
            continue
        positions = sorted(combined)
        cumulative = np.cumsum([combined[position] for position in positions], dtype=float)
        cumulative *= 100.0 / cumulative[-1]
        xs = [0.0] + [position / 1000 for position in positions]
        ys = [0.0] + cumulative.tolist()
        right.step(xs, ys, where="post", color=color, linewidth=1.5,
                   label=str(group) if group < 5 else f"5+ (5–{maximum})")
        curves += 1
    if curves:
        right.legend(title="Full source-set order", frameon=False, fontsize=8, title_fontsize=8)
    else:
        right.text(0.5, 0.5, "No cis span support", transform=right.transAxes,
                   ha="center", va="center", color="0.45")
    right.set_title("Count-weighted cis span distribution", fontsize=11)
    right.set_xlabel("Full cis-set binned span (kb)")
    right.set_ylabel("Cumulative cis support within order group (%)")
    right.set_ylim(0, 102)
    right.set_xlim(left=0)
    right.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda value, _: f"{value:,.3f}".rstrip("0").rstrip(".")))
    for ax in (left, right):
        ax.grid(axis="y", alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
    trans_note = (f"Trans support: {100 * trans_total / total:.2f}% of accepted support; excluded from cis span CDFs."
                  if total > 0 else "Trans support: unavailable because accepted support is zero.")
    fig.suptitle(f"Pore-C source-set distributions · maximum accepted order {maximum}", fontsize=12, y=0.97)
    fig.text(0.08, 0.035,
             "Support is source-set count weight; percentages are not counts of distinct sets.\n"
             "Span uses the full cis set's binned extent; coordinates are 0-based, half-open.\n" + trans_note,
             fontsize=8, va="bottom")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.84, bottom=0.25, wspace=0.28)
    return _save(fig, Path(output), dpi=dpi)
