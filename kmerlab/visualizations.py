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

from .kmer_counter import KmerResult, reverse_complement  # noqa: E402

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


def frequency_histogram(result: KmerResult, bins: int = 40) -> str:
    """Histogram of k-mer occurrence counts (the k-mer spectrum).

    The x-axis is how many times a k-mer occurs; the y-axis is how many distinct
    k-mers share that occurrence count. This is the classic KAT/Jellyfish
    k-mer spectrum plot.
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    counts = list(result.counts.values())
    if not counts:
        ax.text(0.5, 0.5, "No k-mers to display", ha="center", va="center")
        ax.axis("off")
        return _fig_to_data_uri(fig)
    ax.hist(counts, bins=min(bins, max(1, max(counts))), color=_ACCENT, edgecolor="white", linewidth=0.6)
    ax.set_xlabel("k-mer occurrence count")
    ax.set_ylabel("Number of distinct k-mers")
    ax.set_title(f"{result.k}-mer frequency spectrum", loc="left", fontsize=11, fontweight="bold")
    _style(ax)
    if max(counts) > 1:
        ax.set_yscale("log")
    fig.patch.set_alpha(0)
    return _fig_to_data_uri(fig)


def _fcgr_matrix(counts: Counter, k: int) -> list[list[float]]:
    """Build a 2^k x 2^k Frequency Chaos Game Representation matrix.

    Corners are assigned A=(0,0), C=(0,1), G=(1,1), T=(1,0). Each k-mer maps to
    a unique cell by iteratively subdividing the square, following the classic
    CGR construction.
    """
    size = 2 ** k
    matrix = [[0.0] * size for _ in range(size)]
    corner = {"A": (0, 0), "C": (0, 1), "G": (1, 1), "T": (1, 0)}
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

    Returns ``None`` when ``k`` is too large to render a sensible grid (a
    ``2^k x 2^k`` image), otherwise a base64 PNG data URI.
    """
    if result.k > max_k:
        return None
    matrix = _fcgr_matrix(result.counts, result.k)
    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(matrix, cmap="cividis", interpolation="nearest")
    ax.set_title(f"FCGR heatmap (k={result.k})", loc="left", fontsize=11, fontweight="bold", color=_INK)
    fig.patch.set_alpha(0)
    ax.set_xticks([])
    ax.set_yticks([])
    # Label the four corners with their nucleotide.
    ax.text(-0.02, 1.02, "A", transform=ax.transAxes, ha="right", va="bottom", fontsize=10)
    ax.text(1.02, 1.02, "T", transform=ax.transAxes, ha="left", va="bottom", fontsize=10)
    ax.text(-0.02, -0.02, "C", transform=ax.transAxes, ha="right", va="top", fontsize=10)
    ax.text(1.02, -0.02, "G", transform=ax.transAxes, ha="left", va="top", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="frequency")
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
