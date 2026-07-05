"""Tests for k-mer counting, reverse complement, and canonical logic."""

import pytest

from kmerlab.kmer_counter import (
    canonical,
    count_kmers,
    count_kmers_in_sequence,
    reverse_complement,
    validate_k,
)
from kmerlab.sequence_parser import SeqRecord
from collections import Counter


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("AAT") == "ATT"
    assert reverse_complement("GATTACA") == "TGTAATC"


def test_reverse_complement_iupac():
    # R(A/G)<->Y(C/T), S<->S, W<->W, K(G/T)<->M(A/C), B<->V, D<->H, N<->N.
    assert reverse_complement("R") == "Y"
    assert reverse_complement("Y") == "R"
    assert reverse_complement("S") == "S"
    assert reverse_complement("W") == "W"
    assert reverse_complement("K") == "M"
    assert reverse_complement("M") == "K"
    assert reverse_complement("B") == "V"
    assert reverse_complement("D") == "H"
    assert reverse_complement("N") == "N"
    # Combined: reverse then complement each code.
    assert reverse_complement("ACGTN") == "NACGT"
    assert reverse_complement("RYSWKM") == "KMWSRY"
    # Lower-case is preserved.
    assert reverse_complement("acgt") == "acgt"


def test_canonical():
    assert canonical("AAT") == "AAT"  # rc is ATT, AAT < ATT
    assert canonical("TTG") == "CAA"  # rc is CAA, CAA < TTG
    assert canonical("ACGT") == "ACGT"  # palindrome


def test_validate_k():
    assert validate_k("5") == 5
    assert validate_k(3) == 3
    with pytest.raises(ValueError):
        validate_k(0)
    with pytest.raises(ValueError):
        validate_k(-1)
    with pytest.raises(ValueError):
        validate_k("abc")
    with pytest.raises(ValueError):
        validate_k(100)  # exceeds default max_k


def test_count_kmers_basic():
    counts = Counter()
    counted, skipped, invalid = count_kmers_in_sequence("ACGTA", 3, counts)
    assert counted == 3
    assert skipped == 0
    assert invalid == 0
    assert counts == Counter({"ACG": 1, "CGT": 1, "GTA": 1})


def test_count_kmers_repeated():
    counts = Counter()
    count_kmers_in_sequence("AAAA", 2, counts)
    assert counts["AA"] == 3


def test_count_kmers_shorter_than_k():
    counts = Counter()
    counted, skipped, invalid = count_kmers_in_sequence("AC", 3, counts)
    assert counted == 0
    assert len(counts) == 0


def test_invalid_base_handling_excluded():
    counts = Counter()
    counted, skipped, invalid = count_kmers_in_sequence("ACNGT", 2, counts)
    # windows: AC, CN, NG, GT -> CN and NG dropped
    assert counted == 2
    assert skipped == 2
    assert invalid == 1
    assert counts == Counter({"AC": 1, "GT": 1})


def test_invalid_base_handling_included():
    counts = Counter()
    counted, skipped, invalid = count_kmers_in_sequence(
        "ACNGT", 2, counts, include_ambiguous=True
    )
    assert counted == 4
    assert skipped == 0
    assert counts["CN"] == 1
    assert counts["NG"] == 1


def test_canonical_counting():
    # AAA and its rc TTT should merge into canonical AAA.
    records = [SeqRecord(id="r", sequence="AAATTT")]
    plain = count_kmers(records, 3)
    canon = count_kmers(records, 3, canonical_mode=True)
    # plain: AAA, AAT, ATT, TTT
    assert plain.counts["AAA"] == 1
    assert plain.counts["TTT"] == 1
    # canonical: TTT -> AAA, ATT -> AAT ; so AAA=2, AAT=2
    assert canon.counts["AAA"] == 2
    assert canon.counts["AAT"] == 2
    assert "TTT" not in canon.counts


def test_count_kmers_aggregate():
    records = [
        SeqRecord(id="r1", sequence="ACGT"),
        SeqRecord(id="r2", sequence="ACGT"),
    ]
    result = count_kmers(records, 2)
    assert result.n_sequences == 2
    assert result.n_bases == 8
    assert result.counts["AC"] == 2
    assert result.unique_kmers == 3
    assert result.counted_kmers == 6


def test_renamed_metric_fields_present():
    """The result object exposes the clarified field names."""
    result = count_kmers([SeqRecord(id="r", sequence="ACNGT")], 2)
    assert hasattr(result, "counted_kmers")
    assert hasattr(result, "skipped_kmers")
    assert hasattr(result, "invalid_base_count")
    # ACNGT, k=2: AC counted, CN/NG skipped, GT counted; 1 invalid base (N).
    assert result.counted_kmers == 2
    assert result.skipped_kmers == 2
    assert result.invalid_base_count == 1


def test_canonical_with_ambiguous_included():
    """Canonical mode must handle IUPAC codes when ambiguous bases are counted."""
    # ACN and its reverse complement NGT -> canonical is min("ACN","NGT")="ACN".
    counts = Counter()
    count_kmers_in_sequence(
        "ACN", 3, counts, canonical_mode=True, include_ambiguous=True
    )
    assert list(counts) == ["ACN"]
    assert counts["ACN"] == 1


def test_ambiguous_excluded_skips_iupac():
    """With ambiguous excluded, any k-mer with a non-ACGT base is skipped."""
    counts = Counter()
    counted, skipped, invalid = count_kmers_in_sequence(
        "ACRGT", 2, counts, include_ambiguous=False
    )
    # windows AC, CR, RG, GT -> CR and RG dropped (contain R).
    assert counted == 2
    assert skipped == 2
    assert invalid == 1
    assert set(counts) == {"AC", "GT"}
