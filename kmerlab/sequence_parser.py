"""Parsing utilities for FASTA and FASTQ files.

The parsers are intentionally dependency-free (no Biopython required) and are
written as generators so that large files are streamed rather than loaded fully
into memory. Both plain-text and gzip-compressed (``.gz``) files are supported.

Public API
----------
- ``open_maybe_gzip``      : open a path transparently, decompressing ``.gz``.
- ``detect_format``        : sniff FASTA vs FASTQ from the first non-blank line.
- ``parse_fasta``          : yield ``SeqRecord`` for each FASTA entry.
- ``parse_fastq``          : yield ``SeqRecord`` for each FASTQ entry.
- ``parse_file``           : auto-detect format and yield ``SeqRecord``.
- ``ParseError``           : raised for malformed / unsupported input.
"""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass, field
from typing import Iterable, Iterator, TextIO


class ParseError(ValueError):
    """Raised when a file cannot be parsed as FASTA or FASTQ."""


@dataclass
class SeqRecord:
    """A single parsed sequence record."""

    id: str
    sequence: str
    quality: str | None = None  # present only for FASTQ


@dataclass
class ParseResult:
    """Aggregate outcome of parsing an entire file."""

    records: list[SeqRecord] = field(default_factory=list)
    fmt: str = "unknown"          # "fasta" | "fastq"
    n_sequences: int = 0
    n_bases: int = 0
    warnings: list[str] = field(default_factory=list)


# Bytes that mark a gzip file: 0x1f 0x8b.
_GZIP_MAGIC = b"\x1f\x8b"


def _looks_gzip(path: str) -> bool:
    try:
        with open(path, "rb") as handle:
            return handle.read(2) == _GZIP_MAGIC
    except OSError:
        return False


def open_maybe_gzip(path: str) -> TextIO:
    """Open ``path`` as text, transparently decompressing gzip files.

    Detection is by magic bytes (not just the extension) so a mislabelled file
    still opens correctly.
    """
    if path.endswith(".gz") or _looks_gzip(path):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def detect_format_from_text(text: str) -> str:
    """Detect ``"fasta"`` or ``"fastq"`` from a chunk of text.

    Rules mirror the seqkit/FastQC style of sniffing: the first non-blank line
    starting with ``>`` is FASTA, with ``@`` is FASTQ. Raises ``ParseError`` if
    neither marker is found.
    """
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            return "fasta"
        if stripped.startswith("@"):
            return "fastq"
        raise ParseError(
            "Unrecognized format: first content line does not start with "
            "'>' (FASTA) or '@' (FASTQ)."
        )
    raise ParseError("File is empty or contains only blank lines.")


def detect_format(path: str) -> str:
    """Detect the format of a file on disk."""
    with open_maybe_gzip(path) as handle:
        head = handle.read(4096)
    return detect_format_from_text(head)


def parse_fasta(lines: Iterable[str]) -> Iterator[SeqRecord]:
    """Yield :class:`SeqRecord` objects from an iterable of FASTA lines.

    Multi-line sequences are joined. Blank lines are ignored. A header with no
    following sequence yields an empty-sequence record (callers may warn).
    """
    header: str | None = None
    chunks: list[str] = []

    def flush() -> Iterator[SeqRecord]:
        if header is not None:
            yield SeqRecord(id=header, sequence="".join(chunks))

    for raw in lines:
        line = raw.rstrip("\n").rstrip("\r")
        if not line.strip():
            continue
        if line.startswith(">"):
            yield from flush()
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ParseError("FASTA sequence data found before any '>' header.")
            chunks.append(line.strip())
    yield from flush()


def parse_fastq(lines: Iterable[str]) -> Iterator[SeqRecord]:
    """Yield :class:`SeqRecord` objects from an iterable of FASTQ lines.

    Each record is four lines: ``@id``, sequence, ``+`` separator, quality.
    Malformed records raise :class:`ParseError` describing the record number.
    """
    it = iter(lines)
    record_no = 0
    while True:
        try:
            header = next(it)
        except StopIteration:
            return
        # Skip stray blank lines between records.
        if not header.strip():
            continue
        record_no += 1
        if not header.startswith("@"):
            raise ParseError(
                f"FASTQ record {record_no}: header must start with '@' "
                f"(got {header[:20]!r})."
            )
        try:
            seq = next(it).rstrip("\n").rstrip("\r")
            plus = next(it)
            qual = next(it).rstrip("\n").rstrip("\r")
        except StopIteration:
            raise ParseError(
                f"FASTQ record {record_no}: truncated (expected 4 lines)."
            )
        if not plus.startswith("+"):
            raise ParseError(
                f"FASTQ record {record_no}: third line must start with '+' "
                f"(got {plus[:20]!r})."
            )
        if len(seq) != len(qual):
            raise ParseError(
                f"FASTQ record {record_no}: sequence length ({len(seq)}) does "
                f"not match quality length ({len(qual)})."
            )
        yield SeqRecord(id=header[1:].strip(), sequence=seq, quality=qual)


def parse_text(text: str, fmt: str | None = None) -> ParseResult:
    """Parse an in-memory FASTA/FASTQ string into a :class:`ParseResult`."""
    if fmt is None:
        fmt = detect_format_from_text(text)
    lines = text.splitlines()
    parser = parse_fasta if fmt == "fasta" else parse_fastq
    result = ParseResult(fmt=fmt)
    for rec in parser(lines):
        if not rec.sequence:
            result.warnings.append(f"Record '{rec.id}' has an empty sequence.")
        result.records.append(rec)
        result.n_sequences += 1
        result.n_bases += len(rec.sequence)
    if result.n_sequences == 0:
        raise ParseError("No sequences were found in the input.")
    return result


def parse_file(path: str, fmt: str | None = None) -> ParseResult:
    """Parse a file on disk (auto-detecting format) into a :class:`ParseResult`."""
    if fmt is None:
        fmt = detect_format(path)
    with open_maybe_gzip(path) as handle:
        parser = parse_fasta if fmt == "fasta" else parse_fastq
        result = ParseResult(fmt=fmt)
        for rec in parser(handle):
            if not rec.sequence:
                result.warnings.append(f"Record '{rec.id}' has an empty sequence.")
            result.records.append(rec)
            result.n_sequences += 1
            result.n_bases += len(rec.sequence)
    if result.n_sequences == 0:
        raise ParseError("No sequences were found in the input.")
    return result
