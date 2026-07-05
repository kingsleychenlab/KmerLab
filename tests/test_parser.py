"""Tests for FASTA/FASTQ parsing and format detection."""

import gzip

import pytest

from kmerlab.sequence_parser import (
    ParseError,
    detect_format_from_text,
    open_maybe_gzip,
    parse_fasta,
    parse_fastq,
    parse_text,
)

FASTA = ">s1 desc\nACGT\nACGT\n>s2\nGGGG\n"
FASTQ = "@r1\nACGT\n+\nIIII\n@r2\nGGCC\n+\nIIII\n"


def test_parse_fasta_joins_multiline():
    records = list(parse_fasta(FASTA.splitlines()))
    assert len(records) == 2
    assert records[0].id == "s1 desc"
    assert records[0].sequence == "ACGTACGT"
    assert records[1].sequence == "GGGG"


def test_parse_fasta_counts_via_parse_text():
    result = parse_text(FASTA, fmt="fasta")
    assert result.n_sequences == 2
    assert result.n_bases == 12
    assert result.fmt == "fasta"


def test_parse_fastq_basic():
    records = list(parse_fastq(FASTQ.splitlines()))
    assert len(records) == 2
    assert records[0].id == "r1"
    assert records[0].sequence == "ACGT"
    assert records[0].quality == "IIII"


def test_detect_format():
    assert detect_format_from_text(FASTA) == "fasta"
    assert detect_format_from_text(FASTQ) == "fastq"


def test_detect_format_unrecognized():
    with pytest.raises(ParseError):
        detect_format_from_text("hello world\nnot a sequence\n")


def test_detect_format_empty():
    with pytest.raises(ParseError):
        detect_format_from_text("\n\n  \n")


def test_malformed_fastq_length_mismatch():
    bad = "@r1\nACGTACGT\n+\nIII\n"  # seq len 8, qual len 3
    with pytest.raises(ParseError, match="does not match quality length"):
        list(parse_fastq(bad.splitlines()))


def test_malformed_fastq_truncated():
    bad = "@r1\nACGT\n+\n"  # missing quality line
    with pytest.raises(ParseError, match="truncated"):
        list(parse_fastq(bad.splitlines()))


def test_malformed_fastq_bad_header():
    bad = "r1\nACGT\n+\nIIII\n"  # header missing '@'
    with pytest.raises(ParseError, match="must start with '@'"):
        list(parse_fastq(bad.splitlines()))


def test_fasta_data_before_header():
    with pytest.raises(ParseError, match="before any '>' header"):
        list(parse_fasta("ACGT\n>s1\nACGT\n".splitlines()))


def test_no_sequences_raises():
    with pytest.raises(ParseError, match="No sequences"):
        parse_text("\n\n", fmt="fasta")


def test_empty_sequence_warns():
    # A header with no following bases yields a record with an empty sequence.
    result = parse_text(">s1\nACGT\n>s2\n", fmt="fasta")
    assert result.n_sequences == 2
    assert any("empty sequence" in w for w in result.warnings)


def test_gzip_roundtrip(tmp_path):
    path = tmp_path / "seq.fasta.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(FASTA)
    with open_maybe_gzip(str(path)) as fh:
        content = fh.read()
    assert content == FASTA
