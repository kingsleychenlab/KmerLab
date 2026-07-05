"""k-mer counting logic.

Implements a straightforward sliding-window counter (inspired by lh3/kmer-cnt
and Jellyfish's count model) with support for:

- canonical k-mers (a k-mer and its reverse complement counted together),
  with full IUPAC ambiguity-code support in the reverse complement;
- optional inclusion/exclusion of k-mers containing ambiguous/invalid bases
  (e.g. ``N`` and the other IUPAC codes);
- clean tracking of skipped windows and ambiguous/invalid base characters.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable

from .sequence_parser import SeqRecord

# The four unambiguous DNA bases. Everything else is ambiguous/invalid.
VALID_BASES = frozenset("ACGT")

# The full IUPAC nucleotide alphabet (unambiguous + ambiguity codes).
IUPAC_BASES = frozenset("ACGTRYSWKMBDHVN")

# IUPAC-aware complement, upper- and lower-case. Reverse complement of an
# ambiguity code is the complement of the set of bases it represents:
#   R (A/G) <-> Y (C/T), S (C/G) <-> S, W (A/T) <-> W, K (G/T) <-> M (A/C),
#   B (C/G/T) <-> V (A/C/G), D (A/G/T) <-> H (A/C/T), N <-> N.
IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def reverse_complement(seq: str) -> str:
    """Return the IUPAC-aware reverse complement of a nucleotide string.

    Unambiguous bases follow ``A<->T`` / ``C<->G``; ambiguity codes are mapped
    to the complement of the base set they represent (e.g. ``R<->Y``,
    ``K<->M``, ``B<->V``, ``D<->H``, ``S/W/N`` self-complementary). Characters
    outside the IUPAC alphabet are passed through unchanged. Case is preserved.
    """
    return seq.translate(IUPAC_COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    """Return the canonical form: the lexicographically smaller of the k-mer and
    its reverse complement."""
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


@dataclass
class KmerResult:
    """Result of counting k-mers over one or more sequences."""

    k: int
    counts: Counter = field(default_factory=Counter)
    canonical: bool = False
    include_ambiguous: bool = False
    n_sequences: int = 0
    n_bases: int = 0
    counted_kmers: int = 0       # number of k-mer occurrences actually counted
    skipped_kmers: int = 0       # sliding windows skipped (ambiguous/invalid base)
    invalid_base_count: int = 0  # individual non-ACGT characters seen in input

    @property
    def unique_kmers(self) -> int:
        return len(self.counts)

    def top(self, n: int = 20) -> list[tuple[str, int]]:
        """Return the ``n`` most common k-mers as ``(kmer, count)`` pairs."""
        return self.counts.most_common(n)


def validate_k(k: object, max_k: int = 31) -> int:
    """Validate and normalize the ``k`` parameter.

    Returns the integer value of ``k`` or raises ``ValueError`` with a clear
    message. ``k`` must be an integer in ``[1, max_k]``.
    """
    try:
        k_int = int(k)
    except (TypeError, ValueError):
        raise ValueError(f"k must be an integer, got {k!r}.")
    if k_int < 1:
        raise ValueError(f"k must be at least 1, got {k_int}.")
    if k_int > max_k:
        raise ValueError(
            f"k must be at most {max_k}, got {k_int}. Large k values explode "
            f"memory usage; raise max_k explicitly if you really need it."
        )
    return k_int


def count_kmers_in_sequence(
    seq: str,
    k: int,
    counts: Counter,
    canonical_mode: bool = False,
    include_ambiguous: bool = False,
) -> tuple[int, int, int]:
    """Count k-mers in a single sequence (case-normalized to upper).

    Returns ``(counted, skipped, invalid_base_count)`` where ``counted`` is the
    number of k-mer occurrences added to ``counts``, ``skipped`` is the number
    of sliding windows dropped for containing a non-ACGT base (only when
    ``include_ambiguous`` is False), and ``invalid_base_count`` is the number of
    individual non-ACGT characters in the sequence.

    When ``include_ambiguous`` is False, any k-mer containing a non-ACGT
    character is skipped. When True, every window is counted and the canonical
    form (if requested) is computed with IUPAC-aware reverse complement.
    """
    seq = seq.upper()
    n = len(seq)
    invalid_base_count = sum(1 for b in seq if b not in VALID_BASES)
    if n < k:
        return 0, 0, invalid_base_count

    counted = 0
    skipped = 0
    for i in range(n - k + 1):
        kmer = seq[i : i + k]
        if not include_ambiguous and any(b not in VALID_BASES for b in kmer):
            skipped += 1
            continue
        if canonical_mode:
            kmer = canonical(kmer)
        counts[kmer] += 1
        counted += 1
    return counted, skipped, invalid_base_count


def count_kmers(
    records: Iterable[SeqRecord],
    k: int,
    canonical_mode: bool = False,
    include_ambiguous: bool = False,
) -> KmerResult:
    """Count k-mers across an iterable of :class:`SeqRecord`."""
    k = validate_k(k)
    result = KmerResult(
        k=k, canonical=canonical_mode, include_ambiguous=include_ambiguous
    )
    for rec in records:
        result.n_sequences += 1
        result.n_bases += len(rec.sequence)
        counted, skipped, invalid = count_kmers_in_sequence(
            rec.sequence, k, result.counts, canonical_mode, include_ambiguous
        )
        result.counted_kmers += counted
        result.skipped_kmers += skipped
        result.invalid_base_count += invalid
    return result
