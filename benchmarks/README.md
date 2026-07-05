# Benchmark validation fixtures

Tiny, deterministic, **manually inspectable** inputs and their expected k-mer
counts. `tests/test_benchmarks.py` recomputes counts from these FASTA files and
asserts an exact match against the expected CSVs, so the counting logic is
validated against numbers a human worked out by hand.

Each expected CSV has two columns, `kmer,count` (order-independent).

| Input | Mode | Expected |
|-------|------|----------|
| `tiny.fasta` (`ACGTACGT`) | k=3, non-canonical | `tiny.expected.k3.csv` |
| `tiny.fasta` | k=3, canonical | `tiny.expected.k3.canonical.csv` |
| `ambiguous.fasta` (`ACGTNACGT`) | k=3, exclude ambiguous (default) | `ambiguous.expected.k3.csv` |
| `ambiguous.fasta` | k=3, include ambiguous | `ambiguous.expected.k3.include.csv` |

### Worked example — `tiny.fasta`, k=3, non-canonical

`ACGTACGT` has 6 length-3 windows: `ACG CGT GTA TAC ACG CGT`
→ `ACG=2, CGT=2, GTA=1, TAC=1`.

### Worked example — canonical

Each window is replaced by the lexicographically smaller of it and its reverse
complement: `ACG↔CGT` → `ACG`, `GTA↔TAC` → `GTA`.
→ `ACG=4, GTA=2`.

### Worked example — `ambiguous.fasta`, k=3

`ACGTNACGT` has 7 windows: `ACG CGT GTN TNA NAC ACG CGT`.
- **Exclude** (default): the three windows containing `N` are skipped →
  `ACG=2, CGT=2`.
- **Include**: all windows counted → `ACG=2, CGT=2, GTN=1, TNA=1, NAC=1`.
