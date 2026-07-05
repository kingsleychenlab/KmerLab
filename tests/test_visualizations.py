"""Tests for the k-mer spectrum data and FCGR normalization."""

from kmerlab.kmer_counter import count_kmers
from kmerlab.sequence_parser import SeqRecord
from kmerlab.visualizations import (
    _fcgr_matrix,
    fcgr_heatmap,
    kmer_spectrum,
    spectrum_data,
    top_kmers_bar,
)


def test_spectrum_data_exact_multiplicities():
    # ACGTACGT k=3 -> ACG:2, CGT:2, GTA:1, TAC:1.
    # multiplicity 1 -> 2 distinct k-mers; multiplicity 2 -> 2 distinct k-mers.
    result = count_kmers([SeqRecord(id="r", sequence="ACGTACGT")], 3)
    xs, ys = spectrum_data(result)
    assert xs == [1, 2]
    assert ys == [2, 2]
    # y-values sum to the number of distinct k-mers.
    assert sum(ys) == result.unique_kmers


def test_spectrum_data_no_binning():
    # Distinct multiplicities must not be merged into bins.
    seq = "A" * 10  # k=1 -> 'A' occurs 10 times; only multiplicity is 10.
    result = count_kmers([SeqRecord(id="r", sequence=seq)], 1)
    xs, ys = spectrum_data(result)
    assert xs == [10]
    assert ys == [1]


def test_spectrum_renders_data_uri():
    result = count_kmers([SeqRecord(id="r", sequence="ACGTACGT")], 3)
    assert kmer_spectrum(result).startswith("data:image/png")


def test_fcgr_matrix_raw_counts():
    # k=1: A appears 3 times. Corner A=(0,0) -> row size-1, col 0.
    result = count_kmers([SeqRecord(id="r", sequence="AAA")], 1)
    matrix = _fcgr_matrix(result.counts, 1)
    total = sum(cell for row in matrix for cell in row)
    assert total == 3  # raw counts preserved internally


def test_fcgr_normalized_render():
    result = count_kmers([SeqRecord(id="r", sequence="ACGTACGT")], 2)
    uri = fcgr_heatmap(result)
    assert uri is not None
    assert uri.startswith("data:image/png")


def test_fcgr_zero_counts_no_crash():
    # No countable k-mers (all skipped) -> FCGR must render without dividing by 0.
    result = count_kmers([SeqRecord(id="r", sequence="NNNN")], 2)
    assert result.counted_kmers == 0
    uri = fcgr_heatmap(result)
    assert uri is not None
    assert uri.startswith("data:image/png")


def test_fcgr_returns_none_for_large_k():
    result = count_kmers([SeqRecord(id="r", sequence="ACGT" * 5)], 9)
    assert fcgr_heatmap(result) is None


def test_top_bar_empty():
    result = count_kmers([SeqRecord(id="r", sequence="N")], 3)
    assert top_kmers_bar(result).startswith("data:image/png")
