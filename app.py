"""KmerLab Flask application.

A local-only web UI for k-mer frequency analysis of FASTA/FASTQ files. No
network access, database, authentication, or cloud services are used. Uploaded
files are processed in memory and never persisted to disk.

Run with:  python app.py   (then open http://127.0.0.1:5001)
"""

from __future__ import annotations

import os

from flask import (
    Flask,
    Response,
    jsonify,
    render_template,
    request,
    send_from_directory,
)

from kmerlab import __version__
from kmerlab.exports import (
    comparison_to_csv,
    comparison_to_json,
    kmer_counts_to_csv,
    summary_to_json,
)
from kmerlab.kmer_counter import count_kmers, validate_k
from kmerlab.metrics import compare_profiles, summarize
from kmerlab.sequence_parser import (
    MAX_DECOMPRESSED_BYTES,
    ParseError,
    gunzip_bytes_safe,
    parse_text,
)
from kmerlab.visualizations import (
    comparison_bar,
    fcgr_heatmap,
    kmer_spectrum,
    top_kmers_bar,
)

# 50 MB upload ceiling — this keeps the in-memory model honest. Larger files
# should use a future streaming mode (see README > Future improvements).
MAX_CONTENT_LENGTH = 50 * 1024 * 1024
ALLOWED_EXTENSIONS = (".fa", ".fasta", ".fq", ".fastq", ".gz", ".txt")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLES_DIR = os.path.join(BASE_DIR, "samples")

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = MAX_CONTENT_LENGTH


class ApiError(Exception):
    """A user-facing error with an HTTP status code."""

    def __init__(self, message: str, status: int = 400):
        super().__init__(message)
        self.message = message
        self.status = status


def _read_upload(field_name: str) -> str:
    """Read an uploaded file (or a pasted-text fallback) into a decoded string.

    Handles gzip transparently by inspecting the magic bytes. Raises
    :class:`ApiError` with a clean message on any problem.
    """
    file = request.files.get(field_name)
    if file and file.filename:
        filename = file.filename.lower()
        if not filename.endswith(ALLOWED_EXTENSIONS):
            raise ApiError(
                f"Unsupported file type: '{file.filename}'. Allowed: "
                f".fa, .fasta, .fq, .fastq (optionally .gz)."
            )
        raw = file.read()
        if not raw:
            raise ApiError(f"Uploaded file '{file.filename}' is empty.")
        if raw[:2] == b"\x1f\x8b":  # gzip magic bytes
            try:
                raw = gunzip_bytes_safe(raw, limit=MAX_DECOMPRESSED_BYTES)
            except ParseError as exc:
                raise ApiError(str(exc))
        return raw.decode("utf-8", errors="replace")

    # Fallback: pasted text (used by the sample loader and manual entry).
    text = request.form.get(field_name + "_text", "").strip()
    if text:
        return text
    raise ApiError(f"No file or text provided for '{field_name}'.")


def _parse_bool(value: object) -> bool:
    return str(value).lower() in ("1", "true", "on", "yes")


def _analyze_text(text: str, k: int, canonical: bool, include_ambiguous: bool):
    """Parse text and count k-mers, returning ``(parse_result, kmer_result)``."""
    try:
        parsed = parse_text(text)
    except ParseError as exc:
        raise ApiError(str(exc))
    result = count_kmers(
        parsed.records,
        k,
        canonical_mode=canonical,
        include_ambiguous=include_ambiguous,
    )
    return parsed, result


# --------------------------------------------------------------------------- #
# Page routes
# --------------------------------------------------------------------------- #
@app.route("/")
def index():
    return render_template("index.html", version=__version__)


@app.route("/compare")
def compare_page():
    return render_template("compare.html", version=__version__)


@app.route("/about")
def about_page():
    return render_template("about.html", version=__version__)


@app.route("/samples/<path:filename>")
def sample_file(filename: str):
    """Serve bundled sample files so the UI can offer one-click demos."""
    return send_from_directory(SAMPLES_DIR, filename)


@app.route("/api/samples")
def list_samples():
    """List the bundled sample files available for demoing."""
    if not os.path.isdir(SAMPLES_DIR):
        return jsonify([])
    files = sorted(
        f
        for f in os.listdir(SAMPLES_DIR)
        if f.lower().endswith(ALLOWED_EXTENSIONS)
    )
    return jsonify(files)


# --------------------------------------------------------------------------- #
# API routes
# --------------------------------------------------------------------------- #
@app.route("/api/analyze", methods=["POST"])
def api_analyze():
    try:
        k = validate_k(request.form.get("k", "4"))
        canonical = _parse_bool(request.form.get("canonical", "false"))
        include_ambiguous = _parse_bool(request.form.get("include_ambiguous", "false"))
        top_n = min(200, max(1, int(request.form.get("top_n", "20"))))
        want_fcgr = _parse_bool(request.form.get("fcgr", "true"))

        text = _read_upload("file")
        parsed, result = _analyze_text(text, k, canonical, include_ambiguous)

        summary = summarize(result, parsed.records, top_n=top_n)
        summary["format"] = parsed.fmt
        summary["warnings"] = parsed.warnings

        charts = {
            "top_bar": top_kmers_bar(result, top_n),
            "spectrum": kmer_spectrum(result),
        }
        if want_fcgr:
            charts["fcgr"] = fcgr_heatmap(result)

        return jsonify({"summary": summary, "charts": charts})
    except ApiError as exc:
        return jsonify({"error": exc.message}), exc.status
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.route("/api/export/csv", methods=["POST"])
def api_export_csv():
    try:
        k = validate_k(request.form.get("k", "4"))
        canonical = _parse_bool(request.form.get("canonical", "false"))
        include_ambiguous = _parse_bool(request.form.get("include_ambiguous", "false"))
        text = _read_upload("file")
        _, result = _analyze_text(text, k, canonical, include_ambiguous)
        csv_data = kmer_counts_to_csv(result)
        return Response(
            csv_data,
            mimetype="text/csv",
            headers={"Content-Disposition": "attachment; filename=kmer_counts.csv"},
        )
    except ApiError as exc:
        return jsonify({"error": exc.message}), exc.status
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.route("/api/export/json", methods=["POST"])
def api_export_json():
    try:
        k = validate_k(request.form.get("k", "4"))
        canonical = _parse_bool(request.form.get("canonical", "false"))
        include_ambiguous = _parse_bool(request.form.get("include_ambiguous", "false"))
        top_n = min(1000, max(1, int(request.form.get("top_n", "50"))))
        text = _read_upload("file")
        parsed, result = _analyze_text(text, k, canonical, include_ambiguous)
        summary = summarize(result, parsed.records, top_n=top_n)
        summary["format"] = parsed.fmt
        json_data = summary_to_json(summary)
        return Response(
            json_data,
            mimetype="application/json",
            headers={"Content-Disposition": "attachment; filename=kmer_summary.json"},
        )
    except ApiError as exc:
        return jsonify({"error": exc.message}), exc.status
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.route("/api/compare", methods=["POST"])
def api_compare():
    try:
        k = validate_k(request.form.get("k", "4"))
        canonical = _parse_bool(request.form.get("canonical", "false"))
        include_ambiguous = _parse_bool(request.form.get("include_ambiguous", "false"))

        text_a = _read_upload("file_a")
        text_b = _read_upload("file_b")
        parsed_a, result_a = _analyze_text(text_a, k, canonical, include_ambiguous)
        parsed_b, result_b = _analyze_text(text_b, k, canonical, include_ambiguous)

        comparison = compare_profiles(result_a, result_b)
        chart = comparison_bar(result_a, result_b)

        return jsonify(
            {
                "comparison": comparison.to_dict(),
                "file_a": {
                    "format": parsed_a.fmt,
                    "sequences": result_a.n_sequences,
                    "bases": result_a.n_bases,
                    "unique_kmers": result_a.unique_kmers,
                },
                "file_b": {
                    "format": parsed_b.fmt,
                    "sequences": result_b.n_sequences,
                    "bases": result_b.n_bases,
                    "unique_kmers": result_b.unique_kmers,
                },
                "chart": chart,
            }
        )
    except ApiError as exc:
        return jsonify({"error": exc.message}), exc.status
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.route("/api/compare/export", methods=["POST"])
def api_compare_export():
    try:
        fmt = request.form.get("format", "json").lower()
        k = validate_k(request.form.get("k", "4"))
        canonical = _parse_bool(request.form.get("canonical", "false"))
        include_ambiguous = _parse_bool(request.form.get("include_ambiguous", "false"))
        text_a = _read_upload("file_a")
        text_b = _read_upload("file_b")
        _, result_a = _analyze_text(text_a, k, canonical, include_ambiguous)
        _, result_b = _analyze_text(text_b, k, canonical, include_ambiguous)
        comparison = compare_profiles(result_a, result_b)
        if fmt == "csv":
            return Response(
                comparison_to_csv(comparison),
                mimetype="text/csv",
                headers={
                    "Content-Disposition": "attachment; filename=comparison.csv"
                },
            )
        return Response(
            comparison_to_json(comparison),
            mimetype="application/json",
            headers={"Content-Disposition": "attachment; filename=comparison.json"},
        )
    except ApiError as exc:
        return jsonify({"error": exc.message}), exc.status
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.errorhandler(413)
def too_large(_error):
    limit_mb = MAX_CONTENT_LENGTH // (1024 * 1024)
    return (
        jsonify({"error": f"File too large. Maximum upload size is {limit_mb} MB."}),
        413,
    )


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5001, debug=True)
