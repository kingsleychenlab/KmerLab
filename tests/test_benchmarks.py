"""Validation against hand-computed benchmark fixtures in ``benchmarks/``.

These tests recompute k-mer counts from tiny deterministic FASTA files and
assert an exact match against expected CSVs whose numbers were worked out by
hand (see ``benchmarks/README.md``). This guards against silent regressions in
the counting, canonical, and ambiguous-base logic.
"""

import csv
import os

import pytest

from kmerlab.kmer_counter import count_kmers
from kmerlab.sequence_parser import parse_file

BENCH_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "benchmarks")


def _load_expected(name: str) -> dict[str, int]:
    """Read a ``kmer,count`` expected CSV into a dict (order-independent)."""
    path = os.path.join(BENCH_DIR, name)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        return {row["kmer"]: int(row["count"]) for row in reader}


def _counts(fasta: str, k: int, canonical: bool, include_ambiguous: bool) -> dict:
    parsed = parse_file(os.path.join(BENCH_DIR, fasta))
    result = count_kmers(
        parsed.records, k, canonical_mode=canonical, include_ambiguous=include_ambiguous
    )
    return dict(result.counts)


@pytest.mark.parametrize(
    "fasta,k,canonical,include_ambiguous,expected",
    [
        ("tiny.fasta", 3, False, False, "tiny.expected.k3.csv"),
        ("tiny.fasta", 3, True, False, "tiny.expected.k3.canonical.csv"),
        ("ambiguous.fasta", 3, False, False, "ambiguous.expected.k3.csv"),
        ("ambiguous.fasta", 3, False, True, "ambiguous.expected.k3.include.csv"),
    ],
)
def test_benchmark_counts(fasta, k, canonical, include_ambiguous, expected):
    got = _counts(fasta, k, canonical, include_ambiguous)
    want = _load_expected(expected)
    assert got == want


def test_benchmark_files_exist():
    for name in (
        "tiny.fasta",
        "tiny.expected.k3.csv",
        "tiny.expected.k3.canonical.csv",
        "ambiguous.fasta",
        "ambiguous.expected.k3.csv",
    ):
        assert os.path.isfile(os.path.join(BENCH_DIR, name)), name
