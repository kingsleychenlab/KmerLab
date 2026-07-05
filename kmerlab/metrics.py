"""Summary statistics and similarity metrics for k-mer profiles."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable

from .kmer_counter import KmerResult
from .sequence_parser import SeqRecord


def gc_content(records: Iterable[SeqRecord]) -> float:
    """Return the GC fraction (0..1) over all records.

    Only A/C/G/T bases are counted in the denominator so ambiguous characters
    don't distort the ratio. Returns 0.0 for empty input.
    """
    gc = 0
    at = 0
    for rec in records:
        for base in rec.sequence.upper():
            if base in ("G", "C"):
                gc += 1
            elif base in ("A", "T"):
                at += 1
    total = gc + at
    return gc / total if total else 0.0


def base_composition(records: Iterable[SeqRecord]) -> dict:
    """Return per-base counts and fractions for A/C/G/T and 'other'.

    'other' aggregates every non-ACGT character (N, gaps, IUPAC ambiguity
    codes). Fractions are over all characters so they always sum to ~1.0.
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "other": 0}
    for rec in records:
        for base in rec.sequence.upper():
            if base in counts:
                counts[base] += 1
            else:
                counts["other"] += 1
    total = sum(counts.values()) or 1
    fractions = {b: counts[b] / total for b in counts}
    return {"counts": counts, "fractions": fractions}


def sequence_length_stats(records: list[SeqRecord]) -> dict:
    """Return min/max/mean sequence length (0 if empty)."""
    lengths = [len(r.sequence) for r in records]
    if not lengths:
        return {"min": 0, "max": 0, "mean": 0.0}
    return {
        "min": min(lengths),
        "max": max(lengths),
        "mean": sum(lengths) / len(lengths),
    }


def quality_summary(records: Iterable[SeqRecord], phred_offset: int = 33) -> dict | None:
    """Summarize FASTQ Phred quality scores, or ``None`` for FASTA input.

    Assumes Phred+33 encoding (Sanger / Illumina 1.8+). Returns read count,
    average read length, and average/min/max quality across all bases, plus the
    per-read average quality distribution. Records without a quality string
    (i.e. FASTA) are ignored; if none have quality, ``None`` is returned.
    """
    quals = [r.quality for r in records if r.quality is not None]
    if not quals:
        return None
    read_count = len(quals)
    total_len = sum(len(q) for q in quals)
    total_score = 0
    total_bases = 0
    min_q: int | None = None
    max_q: int | None = None
    per_read_avg: list[float] = []
    for q in quals:
        if not q:
            per_read_avg.append(0.0)
            continue
        scores = [ord(c) - phred_offset for c in q]
        s = sum(scores)
        total_score += s
        total_bases += len(scores)
        lo, hi = min(scores), max(scores)
        min_q = lo if min_q is None else min(min_q, lo)
        max_q = hi if max_q is None else max(max_q, hi)
        per_read_avg.append(round(s / len(scores), 3))
    return {
        "reads": read_count,
        "avg_read_length": round(total_len / read_count, 3),
        "avg_quality": round(total_score / total_bases, 3) if total_bases else 0.0,
        "min_quality": min_q if min_q is not None else 0,
        "max_quality": max_q if max_q is not None else 0,
        "per_read_avg_quality": per_read_avg,
    }


def summarize(result: KmerResult, records: list[SeqRecord], top_n: int = 20) -> dict:
    """Build a JSON-serializable summary dict for a single analysis."""
    return {
        "k": result.k,
        "canonical": result.canonical,
        "include_ambiguous": result.include_ambiguous,
        "total_sequences": result.n_sequences,
        "total_bases": result.n_bases,
        "counted_kmers": result.counted_kmers,
        "unique_kmers": result.unique_kmers,
        "skipped_kmers": result.skipped_kmers,
        "invalid_base_count": result.invalid_base_count,
        "gc_content": round(gc_content(records), 6),
        "base_composition": base_composition(records),
        "length_stats": sequence_length_stats(records),
        "quality": quality_summary(records),
        "top_kmers": [
            {"kmer": km, "count": c} for km, c in result.top(top_n)
        ],
    }


def jaccard_similarity(a: Counter, b: Counter) -> float:
    """Jaccard similarity over the *sets* of k-mers.

    ``|A ∩ B| / |A ∪ B|``. Returns 1.0 if both are empty (identical), 0.0 if
    exactly one is empty.
    """
    set_a = set(a)
    set_b = set(b)
    if not set_a and not set_b:
        return 1.0
    union = set_a | set_b
    if not union:
        return 0.0
    intersection = set_a & set_b
    return len(intersection) / len(union)


def cosine_similarity(a: Counter, b: Counter) -> float:
    """Cosine similarity over k-mer *frequency vectors*.

    ``(A · B) / (||A|| * ||B||)``. Returns 0.0 if either vector is all zeros.
    """
    if not a or not b:
        return 0.0
    # Iterate over the smaller dict for the dot product.
    small, large = (a, b) if len(a) <= len(b) else (b, a)
    dot = sum(count * large.get(kmer, 0) for kmer, count in small.items())
    norm_a = math.sqrt(sum(v * v for v in a.values()))
    norm_b = math.sqrt(sum(v * v for v in b.values()))
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return dot / (norm_a * norm_b)


@dataclass
class ComparisonResult:
    """Outcome of comparing two k-mer profiles."""

    k: int
    canonical: bool
    shared: int
    unique_a: int
    unique_b: int
    total_a: int
    total_b: int
    jaccard: float
    cosine: float

    def to_dict(self) -> dict:
        return {
            "k": self.k,
            "canonical": self.canonical,
            "shared_kmers": self.shared,
            "unique_to_a": self.unique_a,
            "unique_to_b": self.unique_b,
            "unique_kmers_a": self.total_a,
            "unique_kmers_b": self.total_b,
            "jaccard_similarity": round(self.jaccard, 6),
            "cosine_similarity": round(self.cosine, 6),
        }


def compare_profiles(a: KmerResult, b: KmerResult) -> ComparisonResult:
    """Compare two :class:`KmerResult` profiles.

    Both profiles must share the same ``k`` and canonical setting for a
    meaningful comparison; a ``ValueError`` is raised otherwise.
    """
    if a.k != b.k:
        raise ValueError(
            f"Cannot compare profiles with different k values ({a.k} vs {b.k})."
        )
    if a.canonical != b.canonical:
        raise ValueError(
            "Cannot compare profiles with different canonical settings."
        )
    set_a = set(a.counts)
    set_b = set(b.counts)
    shared = set_a & set_b
    return ComparisonResult(
        k=a.k,
        canonical=a.canonical,
        shared=len(shared),
        unique_a=len(set_a - set_b),
        unique_b=len(set_b - set_a),
        total_a=len(set_a),
        total_b=len(set_b),
        jaccard=jaccard_similarity(a.counts, b.counts),
        cosine=cosine_similarity(a.counts, b.counts),
    )
