"""Tests for summary statistics and similarity metrics."""

import math
from collections import Counter

import pytest

from kmerlab.kmer_counter import count_kmers
from kmerlab.metrics import (
    compare_profiles,
    cosine_similarity,
    gc_content,
    jaccard_similarity,
    summarize,
)
from kmerlab.sequence_parser import SeqRecord


def test_gc_content():
    assert gc_content([SeqRecord(id="r", sequence="GGCC")]) == 1.0
    assert gc_content([SeqRecord(id="r", sequence="AATT")]) == 0.0
    assert gc_content([SeqRecord(id="r", sequence="ACGT")]) == 0.5
    # ambiguous bases excluded from denominator
    assert gc_content([SeqRecord(id="r", sequence="GCNN")]) == 1.0
    assert gc_content([]) == 0.0


def test_jaccard_identical():
    a = Counter({"AA": 3, "AC": 1})
    b = Counter({"AA": 1, "AC": 5})
    assert jaccard_similarity(a, b) == 1.0  # same sets


def test_jaccard_disjoint():
    a = Counter({"AA": 1})
    b = Counter({"GG": 1})
    assert jaccard_similarity(a, b) == 0.0


def test_jaccard_partial():
    a = Counter({"AA": 1, "AC": 1})
    b = Counter({"AC": 1, "GG": 1})
    # intersection {AC}=1, union {AA,AC,GG}=3
    assert jaccard_similarity(a, b) == pytest.approx(1 / 3)


def test_jaccard_both_empty():
    assert jaccard_similarity(Counter(), Counter()) == 1.0


def test_cosine_identical_direction():
    a = Counter({"AA": 2, "AC": 4})
    b = Counter({"AA": 1, "AC": 2})  # proportional -> cosine 1
    assert cosine_similarity(a, b) == pytest.approx(1.0)


def test_cosine_orthogonal():
    a = Counter({"AA": 1})
    b = Counter({"GG": 1})
    assert cosine_similarity(a, b) == 0.0


def test_cosine_known_value():
    a = Counter({"AA": 1, "AC": 1})
    b = Counter({"AA": 1})
    # dot=1, |a|=sqrt2, |b|=1 -> 1/sqrt2
    assert cosine_similarity(a, b) == pytest.approx(1 / math.sqrt(2))


def test_cosine_empty():
    assert cosine_similarity(Counter(), Counter({"AA": 1})) == 0.0


def test_summarize_shape():
    records = [SeqRecord(id="r", sequence="ACGTACGT")]
    result = count_kmers(records, 3)
    summary = summarize(result, records, top_n=5)
    assert summary["k"] == 3
    assert summary["total_sequences"] == 1
    assert summary["total_bases"] == 8
    assert "gc_content" in summary
    assert len(summary["top_kmers"]) <= 5
    assert summary["top_kmers"][0]["kmer"]


def test_compare_profiles():
    a = count_kmers([SeqRecord(id="a", sequence="ACGTACGT")], 3)
    b = count_kmers([SeqRecord(id="b", sequence="ACGTACGT")], 3)
    comp = compare_profiles(a, b)
    assert comp.jaccard == 1.0
    assert comp.cosine == pytest.approx(1.0)
    assert comp.unique_a == 0
    assert comp.unique_b == 0


def test_compare_profiles_mismatched_k():
    a = count_kmers([SeqRecord(id="a", sequence="ACGT")], 2)
    b = count_kmers([SeqRecord(id="b", sequence="ACGT")], 3)
    with pytest.raises(ValueError, match="different k"):
        compare_profiles(a, b)
