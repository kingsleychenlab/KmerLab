"""Tests for the FASTQ quality summary (Phred+33)."""

from kmerlab.metrics import quality_summary
from kmerlab.sequence_parser import SeqRecord, parse_text


def test_quality_summary_fasta_returns_none():
    records = [SeqRecord(id="r", sequence="ACGT")]  # no quality
    assert quality_summary(records) is None


def test_quality_summary_basic():
    # 'I' = Phred 40, '#' = Phred 2 under Phred+33.
    records = [
        SeqRecord(id="r1", sequence="ACGT", quality="IIII"),
        SeqRecord(id="r2", sequence="ACGT", quality="II##"),
    ]
    q = quality_summary(records)
    assert q["reads"] == 2
    assert q["avg_read_length"] == 4.0
    assert q["max_quality"] == 40
    assert q["min_quality"] == 2
    # 6 bases at 40, 2 bases at 2 -> (6*40 + 2*2) / 8 = 30.5
    assert q["avg_quality"] == 30.5
    assert len(q["per_read_avg_quality"]) == 2
    assert q["per_read_avg_quality"][0] == 40.0


def test_quality_summary_from_parsed_fastq():
    text = "@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nI#\n"
    parsed = parse_text(text)
    q = quality_summary(parsed.records)
    assert q is not None
    assert q["reads"] == 2
    assert q["avg_read_length"] == 3.0  # (4 + 2) / 2


def test_quality_summary_mixed_ignores_fasta_records():
    records = [
        SeqRecord(id="a", sequence="ACGT"),  # no quality
        SeqRecord(id="b", sequence="ACGT", quality="IIII"),
    ]
    q = quality_summary(records)
    assert q["reads"] == 1  # only the record with quality counted
