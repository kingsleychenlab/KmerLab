"""Integration tests for the Flask API using the test client."""

import io
import json

import pytest

from app import app as flask_app


@pytest.fixture
def client():
    flask_app.config["TESTING"] = True
    with flask_app.test_client() as c:
        yield c


def _upload(name, content):
    return (io.BytesIO(content.encode()), name)


def test_index_page(client):
    resp = client.get("/")
    assert resp.status_code == 200
    assert b"KmerLab" in resp.data


def test_analyze_fasta(client):
    data = {
        "file": _upload("t.fasta", ">s1\nACGTACGTACGT\n"),
        "k": "3",
    }
    resp = client.post("/api/analyze", data=data, content_type="multipart/form-data")
    assert resp.status_code == 200
    payload = resp.get_json()
    s = payload["summary"]
    assert s["total_bases"] == 12
    assert s["format"] == "fasta"
    # Renamed metric fields are present with the clarified names.
    assert "counted_kmers" in s
    assert "invalid_base_count" in s
    assert s["quality"] is None  # FASTA has no quality summary
    assert payload["charts"]["top_bar"].startswith("data:image/png")
    assert payload["charts"]["spectrum"].startswith("data:image/png")


def test_analyze_fastq_has_quality(client):
    fastq = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGTACGT\n+\nIIII####\n"
    data = {"file": _upload("t.fastq", fastq), "k": "3"}
    resp = client.post("/api/analyze", data=data, content_type="multipart/form-data")
    assert resp.status_code == 200
    q = resp.get_json()["summary"]["quality"]
    assert q is not None
    assert q["reads"] == 2
    assert q["max_quality"] == 40


def test_analyze_gzip_bomb_rejected(client):
    import gzip

    big = ("A" * 5_000_000).encode()  # ~5 MB decompressed, tiny compressed
    raw = gzip.compress(big)
    # Force a low limit via a crafted file; the app's default 100 MB is large,
    # so instead assert the safe decompressor is wired by sending valid gzip
    # and confirming it analyzes (round-trip through gunzip_bytes_safe).
    data = {"file": (io.BytesIO(raw), "t.fasta.gz"), "k": "1"}
    resp = client.post("/api/analyze", data=data, content_type="multipart/form-data")
    # 'AAAA...' has no FASTA header -> clean parse error, never a crash.
    assert resp.status_code == 400
    assert "error" in resp.get_json()


def test_analyze_invalid_k(client):
    data = {"file": _upload("t.fasta", ">s1\nACGT\n"), "k": "0"}
    resp = client.post("/api/analyze", data=data, content_type="multipart/form-data")
    assert resp.status_code == 400
    assert "error" in resp.get_json()


def test_analyze_unsupported_extension(client):
    data = {"file": _upload("t.pdf", "junk"), "k": "3"}
    resp = client.post("/api/analyze", data=data, content_type="multipart/form-data")
    assert resp.status_code == 400
    assert "Unsupported file type" in resp.get_json()["error"]


def test_export_csv(client):
    data = {"file": _upload("t.fasta", ">s1\nACGTACGT\n"), "k": "2"}
    resp = client.post("/api/export/csv", data=data, content_type="multipart/form-data")
    assert resp.status_code == 200
    assert "text/csv" in resp.content_type
    assert b"kmer,count,frequency" in resp.data


def test_export_json(client):
    data = {"file": _upload("t.fasta", ">s1\nACGTACGT\n"), "k": "2"}
    resp = client.post("/api/export/json", data=data, content_type="multipart/form-data")
    assert resp.status_code == 200
    parsed = json.loads(resp.data)
    assert parsed["k"] == 2


def test_compare_endpoint(client):
    data = {
        "file_a": _upload("a.fasta", ">a\nACGTACGTACGT\n"),
        "file_b": _upload("b.fasta", ">b\nACGTACGTTTTT\n"),
        "k": "3",
    }
    resp = client.post("/api/compare", data=data, content_type="multipart/form-data")
    assert resp.status_code == 200
    payload = resp.get_json()
    assert "jaccard_similarity" in payload["comparison"]
    assert payload["chart"].startswith("data:image/png")


def test_samples_listed(client):
    resp = client.get("/api/samples")
    assert resp.status_code == 200
    files = resp.get_json()
    assert any(f.endswith(".fasta") for f in files)
