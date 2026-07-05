"""Tests for CSV/JSON export helpers and the compare-two-files workflow."""

import io
import json
import csv

from kmerlab.exports import (
    comparison_to_csv,
    comparison_to_json,
    kmer_counts_to_csv,
    summary_to_json,
)
from kmerlab.kmer_counter import count_kmers
from kmerlab.metrics import compare_profiles, summarize
from kmerlab.sequence_parser import SeqRecord, parse_text


def _sample_result():
    records = [SeqRecord(id="r", sequence="ACGTACGT")]
    return count_kmers(records, 2), records


def test_kmer_counts_to_csv():
    result, _ = _sample_result()
    csv_text = kmer_counts_to_csv(result)
    reader = list(csv.reader(io.StringIO(csv_text)))
    assert reader[0] == ["kmer", "count", "frequency"]
    rows = {r[0]: int(r[1]) for r in reader[1:]}
    # ACGTACGT with k=2 -> AC:2 CG:2 GT:2 TA:1
    assert rows["AC"] == 2
    assert rows["TA"] == 1
    # frequencies sum to ~1
    freqs = sum(float(r[2]) for r in reader[1:])
    assert abs(freqs - 1.0) < 1e-6


def test_summary_to_json():
    result, records = _sample_result()
    summary = summarize(result, records, top_n=5)
    text = summary_to_json(summary)
    parsed = json.loads(text)
    assert parsed["k"] == 2
    assert parsed["total_bases"] == 8
    assert isinstance(parsed["top_kmers"], list)


def test_comparison_to_json():
    a = count_kmers([SeqRecord(id="a", sequence="ACGTACGT")], 2)
    b = count_kmers([SeqRecord(id="b", sequence="ACGTTTTT")], 2)
    comp = compare_profiles(a, b)
    parsed = json.loads(comparison_to_json(comp))
    assert "jaccard_similarity" in parsed
    assert "cosine_similarity" in parsed
    assert parsed["k"] == 2


def test_comparison_to_csv():
    a = count_kmers([SeqRecord(id="a", sequence="ACGTACGT")], 2)
    b = count_kmers([SeqRecord(id="b", sequence="ACGTTTTT")], 2)
    comp = compare_profiles(a, b)
    csv_text = comparison_to_csv(comp)
    reader = list(csv.reader(io.StringIO(csv_text)))
    assert reader[0] == ["metric", "value"]
    keys = {r[0] for r in reader[1:]}
    assert "jaccard_similarity" in keys


def test_compare_two_files_end_to_end():
    """Full workflow: parse two files, count, compare, export."""
    text_a = ">a\nACGTACGTACGT\n"
    text_b = ">b\nACGTACGTTTTT\n"
    parsed_a = parse_text(text_a)
    parsed_b = parse_text(text_b)
    result_a = count_kmers(parsed_a.records, 3)
    result_b = count_kmers(parsed_b.records, 3)
    comp = compare_profiles(result_a, result_b)
    assert 0.0 <= comp.jaccard <= 1.0
    assert 0.0 <= comp.cosine <= 1.0
    assert comp.shared > 0  # both share ACG/CGT etc.
    data = json.loads(comparison_to_json(comp))
    assert data["shared_kmers"] == comp.shared
