"""Export helpers: k-mer counts to CSV, summaries and comparisons to JSON."""

from __future__ import annotations

import csv
import io
import json

import pandas as pd

from .kmer_counter import KmerResult
from .metrics import ComparisonResult


def kmer_counts_dataframe(result: KmerResult, sort: bool = True) -> pd.DataFrame:
    """Build a pandas DataFrame of k-mer counts and relative frequencies."""
    total = result.total_kmers or 1
    items = result.counts.most_common() if sort else list(result.counts.items())
    df = pd.DataFrame(items, columns=["kmer", "count"])
    if df.empty:
        df["frequency"] = []
    else:
        df["frequency"] = df["count"] / total
    return df


def kmer_counts_to_csv(result: KmerResult, sort: bool = True) -> str:
    """Serialize k-mer counts to a CSV string with a ``kmer,count,frequency`` header.

    ``frequency`` is the relative frequency (count / total k-mer occurrences).
    """
    df = kmer_counts_dataframe(result, sort=sort)
    return df.to_csv(index=False, float_format="%.8f")


def summary_to_json(summary: dict, indent: int = 2) -> str:
    """Serialize an analysis summary dict to a JSON string."""
    return json.dumps(summary, indent=indent, sort_keys=False)


def comparison_to_json(comparison: ComparisonResult, indent: int = 2) -> str:
    """Serialize a comparison result to a JSON string."""
    return json.dumps(comparison.to_dict(), indent=indent, sort_keys=False)


def comparison_to_csv(comparison: ComparisonResult) -> str:
    """Serialize a comparison result to a flat ``metric,value`` CSV string."""
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["metric", "value"])
    for key, value in comparison.to_dict().items():
        writer.writerow([key, value])
    return buf.getvalue()
