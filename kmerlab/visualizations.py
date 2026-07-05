"""Matplotlib-based visualizations, rendered to base64 PNG data URIs.

All figures are rendered with the non-interactive ``Agg`` backend so the module
works headlessly (no display server needed) and returns strings that can be
embedded directly into an ``<img src="...">`` tag.
"""

from __future__ import annotations

import base64
import io
from collections import Counter

import matplotlib

matplotlib.use("Agg")  # headless backend; must be set before pyplot import
import matplotlib.pyplot as plt  # noqa: E402

from .kmer_counter import KmerResult  # noqa: E402

# Palette harmonized with the app's nucleotide colour system.
_ACCENT = "#2f73e0"   # base-C blue: primary single-series colour
_ACCENT_B = "#e0483c"  # base-T red: second series in comparisons
_INK = "#0d1b24"
_GRID = "#dce6eb"


def _style(ax):
    """Apply a clean, instrument-like axis style shared by all figures."""
    ax.set_facecolor("white")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(_GRID)
    ax.tick_params(colors=_INK, labelsize=8, length=0)
    ax.title.set_color(_INK)
    ax.xaxis.label.set_color(_INK)
    ax.yaxis.label.set_color(_INK)
    ax.grid(axis="both", color=_GRID, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)


def _fig_to_data_uri(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def top_kmers_bar(result: KmerResult, n: int = 20) -> str:
    """Horizontal bar chart of the top ``n`` k-mers by count."""
    top = result.top(n)
    fig, ax = plt.subplots(figsize=(7, max(3, 0.32 * len(top))))
    if not top:
        ax.text(0.5, 0.5, "No k-mers to display", ha="center", va="center")
        ax.axis("off")
        return _fig_to_data_uri(fig)
    labels = [km for km, _ in top][::-1]
    values = [c for _, c in top][::-1]
    ax.barh(labels, values, color=_ACCENT, height=0.72)
    ax.set_xlabel("Count")
    ax.set_title(f"Top {len(top)} {result.k}-mers", loc="left", fontsize=11, fontweight="bold")
    _style(ax)
    ax.tick_params(axis="y", labelsize=8)
    fig.patch.set_alpha(0)
    return _fig_to_data_uri(fig)


def spectrum_data(result: KmerResult) -> tuple[list[int], list[int]]:
    """Return the exact k-mer spectrum as ``(multiplicities, distinct_counts)``.

    ``multiplicities[i]`` is a k-mer occurrence count that actually appears in
    the data; ``distinct_counts[i]`` is how many distinct k-mers have exactly
    that multiplicity. No binning is applied, so integer multiplicities are
    never merged.
    """
    multiplicity_counts = Counter(result.counts.values())
    xs = sorted(multiplicity_counts)
    ys = [multiplicity_counts[m] for m in xs]
    return xs, ys


def kmer_spectrum(result: KmerResult) -> str:
    """True k-mer spectrum: distinct-k-mer count per exact multiplicity.

    x-axis = k-mer multiplicity (how many times a k-mer occurs), y-axis =
    number of distinct k-mers with that multiplicity. Unlike a binned
    histogram, each integer multiplicity is its own bar, so multiplicities are
    never merged or misrepresented.
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    xs, ys = spectrum_data(result)
    if not xs:
        ax.text(0.5, 0.5, "No k-mers to display", ha="center", va="center")
        ax.axis("off")
        return _fig_to_data_uri(fig)
    # When the multiplicity range is wide, bars become too thin to see, so fall
    # back to a stem/line plot. (Bars on a log x-axis would render with
    # distorted widths, so we never combine bars with a log x-axis.)
    if max(xs) > 60:
        ax.plot(xs, ys, marker="o", markersize=3, linewidth=1, color=_ACCENT)
        ax.set_xscale("log")
    else:
        ax.bar(xs, ys, width=0.9, color=_ACCENT, edgecolor="white", linewidth=0.4)
    ax.set_xlabel("k-mer multiplicity")
    ax.set_ylabel("Number of distinct k-mers")
    ax.set_title(f"K-mer Spectrum (k={result.k})", loc="left", fontsize=11, fontweight="bold")
    _style(ax)
    if max(ys) > 10:
        ax.set_yscale("log")
    fig.patch.set_alpha(0)
    return _fig_to_data_uri(fig)


def _fcgr_matrix(counts: Counter, k: int) -> list[list[float]]:
    """Build a 2^k x 2^k Frequency Chaos Game Representation matrix.

    Corners (as drawn by :func:`fcgr_heatmap`, where row 0 is the top): A =
    top-left, C = bottom-left, G = bottom-right, T = top-right. Each corner is
    ``(cx, cy)`` with ``cx`` the horizontal half (0 = left, 1 = right) and
    ``cy`` the vertical half (0 = bottom, 1 = top). Each k-mer maps to a unique
    cell by iteratively subdividing the square (the classic CGR construction),
    so the corner map here stays consistent with the corner labels in the plot.
    """
    size = 2 ** k
    matrix = [[0.0] * size for _ in range(size)]
    corner = {"A": (0, 1), "C": (0, 0), "G": (1, 0), "T": (1, 1)}
    for kmer, value in counts.items():
        x0, y0, x1, y1 = 0.0, 0.0, 1.0, 1.0
        valid = True
        for base in kmer:
            if base not in corner:
                valid = False
                break
            cx, cy = corner[base]
            midx = (x0 + x1) / 2
            midy = (y0 + y1) / 2
            if cx == 0:
                x1 = midx
            else:
                x0 = midx
            if cy == 0:
                y1 = midy
            else:
                y0 = midy
        if not valid:
            continue
        col = min(size - 1, int(((x0 + x1) / 2) * size))
        row = min(size - 1, int((1 - (y0 + y1) / 2) * size))
        matrix[row][col] += value
    return matrix


def fcgr_heatmap(result: KmerResult, max_k: int = 8) -> str | None:
    """Frequency Chaos Game Representation heatmap for DNA k-mers.

    Cells are normalized to relative frequency (each raw count divided by the
    total number of counted k-mers), so the colour scale is comparable across
    datasets. If no k-mers were counted, an all-zero grid is drawn cleanly
    rather than crashing on a divide-by-zero.

    Returns ``None`` when ``k`` exceeds ``max_k`` (a ``2^k x 2^k`` image grows
    too large to be useful), otherwise a base64 PNG data URI.
    """
    if result.k > max_k:
        return None
    matrix = _fcgr_matrix(result.counts, result.k)
    total = sum(result.counts.values())
    if total > 0:
        matrix = [[cell / total for cell in row] for row in matrix]
    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(matrix, cmap="cividis", interpolation="nearest", vmin=0)
    ax.set_title(
        f"FCGR normalized frequency heatmap (k={result.k})",
        loc="left", fontsize=11, fontweight="bold", color=_INK,
    )
    fig.patch.set_alpha(0)
    ax.set_xticks([])
    ax.set_yticks([])
    # Label the four corners with their nucleotide.
    ax.text(-0.02, 1.02, "A", transform=ax.transAxes, ha="right", va="bottom", fontsize=10)
    ax.text(1.02, 1.02, "T", transform=ax.transAxes, ha="left", va="bottom", fontsize=10)
    ax.text(-0.02, -0.02, "C", transform=ax.transAxes, ha="right", va="top", fontsize=10)
    ax.text(1.02, -0.02, "G", transform=ax.transAxes, ha="left", va="top", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="normalized frequency")
    return _fig_to_data_uri(fig)


def comparison_bar(a: KmerResult, b: KmerResult, n: int = 15) -> str:
    """Grouped bar chart comparing counts of the top k-mers across two profiles."""
    top_kmers = [km for km, _ in (a.counts + b.counts).most_common(n)]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    if not top_kmers:
        ax.text(0.5, 0.5, "No k-mers to display", ha="center", va="center")
        ax.axis("off")
        return _fig_to_data_uri(fig)
    x = range(len(top_kmers))
    width = 0.4
    a_vals = [a.counts.get(km, 0) for km in top_kmers]
    b_vals = [b.counts.get(km, 0) for km in top_kmers]
    ax.bar([i - width / 2 for i in x], a_vals, width, label="File A", color=_ACCENT)
    ax.bar([i + width / 2 for i in x], b_vals, width, label="File B", color=_ACCENT_B)
    ax.set_xticks(list(x))
    ax.set_xticklabels(top_kmers, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Count")
    ax.set_title(f"Top {len(top_kmers)} {a.k}-mers: File A vs File B", loc="left", fontsize=11, fontweight="bold")
    _style(ax)
    ax.legend(frameon=False, fontsize=9)
    fig.patch.set_alpha(0)
    return _fig_to_data_uri(fig)
