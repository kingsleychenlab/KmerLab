# KmerLab

**Local k-mer analysis dashboard for FASTA/FASTQ files.**

KmerLab is a local-first, educational and validation-focused bioinformatics
utility for small-to-medium FASTA/FASTQ datasets. Upload (or select) a sequence
file, choose a length *k*, and get k-mer counts, base composition, a true k-mer
spectrum, a normalized FCGR heatmap, and two-file similarity — computed entirely
on your own machine.

> **Runs 100% locally.** No cloud, no database, no login, no API keys, no
> internet requirement after installation, no paid services. Uploaded files are
> processed in memory and are never written to disk.

**Status — read this honestly:** KmerLab is *validated* against hand-computed
benchmark fixtures and covered by an automated test suite, and it uses correct,
standard definitions for every metric it reports. It is **not** a research-grade
counter for large genomes: it does exact in-memory counting (no streaming, no
disk-backed hashing, no probabilistic sketches) and is intended for teaching,
exploration, and validating small datasets — not for whole-genome-scale runs.
See [Limitations](#limitations) and
[Is it research-grade?](#is-it-research-grade).

---

## What it does

- Parses **FASTA** and **FASTQ** (auto-detected), plain or **gzip** (`.gz`).
- Counts k-mers with a sliding window (exact `collections.Counter`).
- **Canonical reverse-complement mode** with full **IUPAC** ambiguity support.
- Include or exclude k-mers containing ambiguous/invalid bases.
- Reports: sequences, bases, **counted k-mers**, **unique k-mers**, **skipped
  windows**, **ambiguous/invalid base count**, GC content (over A/C/G/T), base
  composition, and sequence-length stats.
- **True k-mer spectrum** (distinct k-mers per exact multiplicity — no binning).
- **FCGR normalized frequency heatmap** for k ≤ 8.
- Optional **FASTQ quality summary** (Phred+33): reads, avg length, avg/min/max
  quality, per-read average quality.
- **Compare two files**: shared / unique k-mers, **Jaccard** and **cosine**
  similarity, side-by-side chart.
- **Exports**: k-mer counts → CSV, summary → JSON, comparison → JSON/CSV.

## What it does *not* do

- No large-file / streaming mode — everything is held in memory.
- No probabilistic sketching (MinHash, HyperLogLog) or disk-backed counting.
- No protein / amino-acid k-mers (nucleotide only).
- No alignment, assembly, taxonomic classification, or variant calling.
- No network access, database, accounts, or telemetry of any kind.
- FASTQ quality assumes **Phred+33**; legacy Phred+64 is not auto-detected.

---

## Installation

```bash
git clone <your-fork-url> KmerLab
cd KmerLab

python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install -r requirements.txt          # flexible versions
# or, for an exact reproducible environment:
pip install -r requirements-lock.txt
```

## Run

```bash
python app.py
```

Open **http://127.0.0.1:5001** (binds to localhost only).

## Test

```bash
pytest
```

The suite covers parsing, strict FASTQ error cases, k-mer counting, IUPAC
reverse complement, canonical + ambiguous modes, similarity metrics, the k-mer
spectrum, FCGR normalization (including the zero-count case), the gzip
decompression limit, the FASTQ quality summary, CSV/JSON exports, the Flask API,
and exact-match **benchmark validation** (`benchmarks/`).

---

## Example workflow

1. Open the **Analyzer** page.
2. Load the sample `example.fasta` (or upload your own file).
3. Set **k = 4**, optionally enable **Canonical reverse-complement mode**.
4. Click **Analyze** and inspect the cards, base composition, spectrum, FCGR
   heatmap, and frequency table.
5. Export counts (CSV) or the summary (JSON).
6. Use the **Compare** page with `compare_a.fasta` / `compare_b.fasta` for
   Jaccard/cosine similarity.

The core library also works standalone:

```python
from kmerlab.sequence_parser import parse_file
from kmerlab.kmer_counter import count_kmers
from kmerlab.metrics import summarize

parsed = parse_file("samples/example.fasta")
result = count_kmers(parsed.records, k=4, canonical_mode=True)
print(summarize(result, parsed.records))
```

---

## Method

**k-mer.** A length-*k* substring extracted with a one-base sliding window.
Counting k-mers gives an alignment-free fingerprint of a dataset.

**Canonical k-mers.** DNA is double-stranded, so a k-mer and its reverse
complement usually describe the same physical sequence. The reverse complement
swaps A↔T and C↔G and reverses the string; the **canonical** form is the
lexicographically smaller of the two. Canonical mode merges both orientations.

**IUPAC ambiguous bases.** Sequence files contain more than A/C/G/T. IUPAC codes
encode uncertainty (`N`=any, `R`=A/G, `Y`=C/T, …). KmerLab knows the full
alphabet `ACGTRYSWKMBDHVN` and complements ambiguity codes correctly (`R↔Y`,
`K↔M`, `B↔V`, `D↔H`; `S`, `W`, `N` self-complementary). By default, any k-mer
window containing a non-ACGT base is **skipped** (counted as a *skipped
window*); enable *Include ambiguous IUPAC bases* to count them.

**Strict FASTQ parser.** KmerLab uses a strict 4-line FASTQ parser: every record
must be `@header` / sequence / `+` separator / quality, and the sequence and
quality lengths must match. Malformed records raise a clear error naming the
offending record (incomplete record, missing `@`, missing `+`, or length
mismatch) instead of being silently skipped.

**k-mer spectrum.** Plots the number of *distinct* k-mers at each exact
occurrence count (x = multiplicity, y = number of distinct k-mers). KmerLab uses
exact integer multiplicities — **not histogram bins** — so counts are never
merged.

**FCGR normalized frequency heatmap.** Maps k-mer frequencies onto a
2ᵏ × 2ᵏ grid (corners A/C/G/T). Each cell is **normalized** to a relative
frequency (raw count ÷ total counted k-mers) so scales are comparable across
datasets; the zero-count case is handled without dividing by zero. Rendered for
**k ≤ 8** because the image is a 2ᵏ × 2ᵏ grid.

**Jaccard similarity.** `|A ∩ B| / |A ∪ B|` over the *sets* of k-mers.

**Cosine similarity.** `(A · B) / (‖A‖ · ‖B‖)` over k-mer *frequency vectors*.

---

## Folder structure

```
KmerLab/
├── app.py                    # Flask app: routes + API endpoints
├── requirements.txt          # flexible dependencies
├── requirements-lock.txt     # exact pinned dependencies
├── pytest.ini
├── LICENSE                   # MIT
├── README.md
├── kmerlab/                  # core library (importable, Flask-independent)
│   ├── __init__.py
│   ├── sequence_parser.py    # FASTA/FASTQ parsing, format detection, safe gzip
│   ├── kmer_counter.py       # counting, IUPAC reverse complement, canonical
│   ├── metrics.py            # GC, composition, quality, Jaccard, cosine
│   ├── visualizations.py     # matplotlib charts -> base64 PNG data URIs
│   └── exports.py            # CSV / JSON export helpers (pandas)
├── templates/                # Jinja2 HTML (base, index, compare, about)
├── static/                   # style.css, main.js (vanilla JS)
├── samples/                  # example / invalid / compare sample files
├── benchmarks/               # tiny deterministic files + expected counts
├── tests/                    # pytest suite
└── .github/workflows/        # CI: pytest on Python 3.10 / 3.11 / 3.12
```

---

## Limitations

- **In-memory, exact counting.** Memory scales with the number of distinct
  k-mers and with *k*. Not suitable for whole genomes or deep sequencing runs.
- Web UI caps uploads at **50 MB compressed** and **100 MB decompressed**
  (gzip-bomb protection); larger inputs are rejected cleanly.
- *k* is capped at **31**; the FCGR heatmap at **k ≤ 8**.
- Nucleotide sequences only (no protein k-mers).
- FASTQ quality assumes **Phred+33**.
- Similarity metrics are only meaningful when both files use the same *k* and
  canonical setting (enforced by the app).
- `debug=True` is enabled for local development; disable it if you ever expose
  the app beyond localhost.

## Future improvements

- Streaming / bounded-memory counting mode.
- Broader compression support (bzip2, zstd).
- Protein / amino-acid k-mers.
- More FASTQ QC (per-position quality, adapter scan).
- MinHash sketching for fast approximate comparison.
- A standalone CLI and Docker image.
- PDF report export.

## Is it research-grade?

**Partially, and only for what it claims.** The metric definitions are correct
and standard, the results are deterministic and validated against hand-computed
benchmarks, the code is tested across Python 3.10–3.12 in CI, and the tool is
honest about what it computes. That makes it a solid **validated educational and
validation-focused** utility.

It is **not** research-grade as a high-performance counter: it does not stream,
shard, or sketch, so it cannot handle genome-scale data, and it has not been
benchmarked for speed or memory against tools like KMC or Jellyfish. Calling it
a "fully research-grade k-mer counter" would be inaccurate, and this README
deliberately does not.

## License

[MIT License](LICENSE).
