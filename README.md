# KmerLab

**Local k-mer analysis dashboard for FASTA/FASTQ files.**

KmerLab is a local-first, educational and validation-focused bioinformatics
utility for small-to-medium FASTA/FASTQ datasets. Load a sequence file, choose a
k-mer length *k*, and get k-mer counts, base composition, a true k-mer spectrum,
a normalized FCGR heatmap, optional FASTQ quality metrics, and two-file
similarity — computed entirely on your own machine, in a clean web dashboard.

> **Runs 100% locally.** No cloud, no database, no login, no API keys, no
> internet requirement after installation, no paid services. Uploaded files are
> processed in memory and are never written to disk.

**Status — read this honestly.** KmerLab is *validated* against hand-computed
benchmark fixtures and covered by an automated test suite, and it uses correct,
standard definitions for every metric it reports. It is **not** a research-grade
counter for large genomes: it does exact in-memory counting (no streaming, no
disk-backed hashing, no probabilistic sketches) and is meant for teaching,
exploration, and validating small datasets — not whole-genome-scale runs. See
[Limitations](#limitations) and [Is it research-grade?](#is-it-research-grade).

---

## Table of contents

- [What it does / doesn't](#what-it-does)
- [Install](#install)
- [Run](#run)
- [Quick start (2 minutes)](#quick-start-2-minutes)
- [Using the Analyzer](#using-the-analyzer)
- [Using Compare](#using-compare)
- [Understanding every output](#understanding-every-output)
- [What you'd actually use it for](#what-youd-actually-use-it-for)
- [Methods & math](#methods--math)
- [Use it as a Python library](#use-it-as-a-python-library)
- [HTTP API reference](#http-api-reference)
- [Export formats](#export-formats)
- [Testing & benchmarks](#testing--benchmarks)
- [Folder structure](#folder-structure)
- [Limitations](#limitations)
- [Future improvements](#future-improvements)
- [Is it research-grade?](#is-it-research-grade)
- [License](#license)

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
  similarity, with a side-by-side chart.
- **Exports**: k-mer counts → CSV, summary → JSON, comparison → JSON/CSV.

### What it does *not* do

- No large-file / streaming mode — everything is held in memory.
- No probabilistic sketching (MinHash, HyperLogLog) or disk-backed counting.
- No protein / amino-acid k-mers (nucleotide only).
- No alignment, assembly, taxonomic classification, or variant calling.
- No network access, database, accounts, or telemetry of any kind.
- FASTQ quality assumes **Phred+33**; legacy Phred+64 is not auto-detected.

---

## Install

Requires **Python 3.10+**.

```bash
git clone https://github.com/kingsleychenlab/KmerLab.git
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

Then open **http://127.0.0.1:5001** in your browser. The server binds to
`127.0.0.1` (localhost) only, so it is not reachable from other machines.

> Port 5001 is used by default (port 5000 is taken by AirPlay Receiver on
> macOS). To change it, edit the `app.run(...)` call at the bottom of `app.py`.

---

## Quick start (2 minutes)

1. Run `python app.py` and open http://127.0.0.1:5001.
2. On the **Analyzer** page, use the **"Load a sample"** dropdown and pick
   `example.fasta` (no file needed — sample files are bundled).
3. Leave **k = 4**, click **Analyze sequences**.
4. Read the summary cards, the A/C/G/T composition bar, the k-mer spectrum, the
   FCGR heatmap, and the top-k-mer table.
5. Click **Export k-mer counts (CSV)** to download the full count table.
6. Open the **Compare** page, load `compare_a.fasta` and `compare_b.fasta`,
   click **Compare profiles**, and read the Jaccard / cosine scores.

---

## Using the Analyzer

The Analyzer page turns one file into a full k-mer profile.

**1. Provide input.** Either drag-and-drop a file onto the dropzone / click to
browse, or pick a bundled file from **Load a sample**. Accepted:
`.fa`, `.fasta`, `.fq`, `.fastq`, optionally gzip-compressed (`.gz`). Format
(FASTA vs FASTQ) is detected automatically from the first line.

**2. Set the options.**

| Control | What it does | When to change it |
|---------|--------------|-------------------|
| **k (length)** | Length of each k-mer window (1–31). | Small *k* (2–4) for composition/GC-style views; larger *k* (11–21) for genomic uniqueness. See [Choosing k](#choosing-k). |
| **Top N** | How many of the most frequent k-mers to chart and table (1–200). | Raise it to see more of the head of the distribution. |
| **Canonical reverse-complement mode** | Counts each k-mer together with its reverse complement (IUPAC-aware). | On for double-stranded DNA where strand shouldn't matter (the usual case). Off to keep strands separate. |
| **Include ambiguous IUPAC bases** | Counts k-mers containing non-ACGT codes (N, R, Y, …) instead of skipping them. | On when you want ambiguity codes represented; off (default) for clean A/C/G/T-only counts. |
| **FCGR normalized frequency heatmap** | Renders the FCGR image (k ≤ 8 only). | Turn off to skip the heatmap for a faster run. |

**3. Click Analyze.** Results appear below; the page scrolls to them. Any
problem (bad file, invalid *k*, empty/malformed input, oversized upload) shows a
single clear red error panel instead of a crash.

**4. Export.** *Export k-mer counts (CSV)* downloads every k-mer with its count
and relative frequency. *Export summary (JSON)* downloads the full metric
summary. Exports use the same options currently set in the form.

---

## Using Compare

The Compare page measures how similar two datasets are by their k-mers.

1. Provide **File A** and **File B** (upload or sample). Both are analyzed with
   the **same** *k* and canonical setting so the comparison is meaningful.
2. Choose *k* and options, click **Compare profiles**.
3. Read the result cards: **Jaccard** and **cosine** similarity, plus **shared**,
   **unique-to-A**, and **unique-to-B** k-mer counts, and per-file stats.
4. The side-by-side bar chart shows the most frequent k-mers in each file.
5. *Export comparison (JSON/CSV)* downloads the numbers.

**Reading the two scores together:** Jaccard is presence/absence (set overlap);
cosine is composition (frequency vectors). High Jaccard **and** high cosine →
the files share k-mers *and* in similar proportions. High Jaccard but low cosine
→ they share many k-mers but at very different frequencies.

---

## Understanding every output

**Summary cards**

| Card | Meaning |
|------|---------|
| Format | FASTA or FASTQ (auto-detected). |
| Sequences | Number of records parsed. |
| Total bases | Total characters across all sequences. |
| Counted k-mers | Number of k-mer *occurrences* actually counted. |
| Unique k-mers | Number of *distinct* k-mers. |
| GC content (A/C/G/T) | Fraction of G/C among A/C/G/T bases only. |
| Skipped windows | Sliding windows dropped for containing a non-ACGT base (only when *Include ambiguous* is off). |
| Ambiguous / invalid bases | Count of individual non-ACGT characters seen in the input. |

**Base composition bar** — a stacked A/C/G/T/other bar with exact percentages,
using KmerLab's nucleotide colour system (A green, C blue, G amber, T red).

**FASTQ quality summary** (FASTQ only) — reads, average read length, and
average / min / max Phred quality (Phred+33). Hidden for FASTA.

**K-mer spectrum** — x-axis is k-mer multiplicity (how many times a k-mer
occurs), y-axis is how many *distinct* k-mers have that multiplicity. Exact
integer multiplicities are used (no histogram binning), so counts are never
merged. The `y` values sum to the number of unique k-mers.

**FCGR normalized frequency heatmap** — k-mer frequencies mapped onto a
2ᵏ × 2ᵏ grid (corners A top-left, C bottom-left, G bottom-right, T top-right),
each cell normalized to a relative frequency. Rendered for k ≤ 8.

**Top k-mer frequency table** — the most frequent k-mers with counts; each base
letter is coloured by nucleotide.

---

## What you'd actually use it for

- **Teaching / learning** what k-mers, canonical k-mers, reverse complement,
  and k-mer spectra are, with immediate visual feedback.
- **Quick QC of a small FASTA/FASTQ** before a larger pipeline: check GC%, base
  composition, read-length and quality ranges, and spot obvious contamination or
  low-complexity regions.
- **Comparing two small datasets** (e.g. two amplicons, two plasmids, two short
  assemblies) for k-mer overlap without running a full aligner.
- **Sanity-checking another tool's counts** on a tiny input where you can verify
  the numbers by hand (this is exactly what the `benchmarks/` fixtures do).
- **Exploring composition signatures** via the FCGR heatmap and spectrum shape.

It is **not** for whole genomes, deep sequencing runs, or production pipelines —
see [Limitations](#limitations).

---

## Methods & math

### k-mer

A **k-mer** is a length-*k* substring taken with a one-base sliding window. For
`ACGTA` with *k* = 3 the k-mers are `ACG`, `CGT`, `GTA`. Counting k-mers gives an
alignment-free "fingerprint" of a dataset.

### Choosing k

There is no universally correct *k*; it depends on the question.

- **Small k (1–4):** captures base/short-motif composition; nearly every k-mer
  appears, so it is good for GC/composition-style summaries and FCGR textures.
- **Medium k (5–10):** motif-level structure; a reasonable default for
  exploration (the app defaults to 4).
- **Large k (11–21+):** k-mers become mostly unique across a genome, which is
  what assemblers and comparison tools rely on. Requires more memory (distinct
  k-mers grow quickly). KmerLab caps *k* at 31.

A practical heuristic: pick the smallest *k* at which most k-mers are unique for
your purpose. The k-mer spectrum helps — a large single-copy peak means most
k-mers occur about the same number of times.

### Reverse complement & canonical k-mers

DNA is double-stranded, so a k-mer and its reverse complement usually describe
the same physical sequence read from opposite strands. The **reverse
complement** swaps A↔T and C↔G and reverses the string (`AAT` → `ATT`). The
**canonical** k-mer is the lexicographically smaller of a k-mer and its reverse
complement; canonical mode merges both strand orientations into one count, as in
KMC and Jellyfish.

### IUPAC ambiguous bases

Real files contain more than A/C/G/T. IUPAC codes encode uncertainty: `N` = any
base, `R` = A/G, `Y` = C/T, and so on. KmerLab knows the full alphabet
`ACGTRYSWKMBDHVN` and reverse-complements ambiguity codes correctly
(`R↔Y`, `K↔M`, `B↔V`, `D↔H`; `S`, `W`, `N` are self-complementary). By default,
any window containing a non-ACGT base is **skipped** (and reported as a *skipped
window*); enable *Include ambiguous IUPAC bases* to count them instead.

### Strict FASTQ parser

KmerLab uses a strict 4-line FASTQ parser: every record must be `@header` /
sequence / `+` separator / quality, and the sequence and quality lengths must
match. Malformed records raise a clear error naming the offending record
(incomplete record, missing `@`, missing `+`, or length mismatch) rather than
being silently skipped.

### FASTQ quality (Phred+33)

Quality scores are decoded as `ord(char) - 33` (Sanger / Illumina 1.8+). The
summary reports read count, average read length, and average / min / max quality
across all bases, plus each read's average quality.

### k-mer spectrum

Plots the number of *distinct* k-mers at each exact occurrence count. Its shape
is diagnostic: a peak at multiplicity 1 that dominates usually means sequencing
errors or high diversity; a clear peak at higher multiplicity indicates typical
coverage; secondary peaks can indicate repeats or ploidy. KmerLab uses exact
integer multiplicities, never bins.

### FCGR normalized frequency heatmap

Frequency Chaos Game Representation maps k-mer frequencies onto a 2ᵏ × 2ᵏ grid
by recursively subdividing a square whose corners are A/C/G/T. Each cell is
**normalized** to a relative frequency (raw count ÷ total counted k-mers) so the
colour scale is comparable across datasets; the zero-count case is handled
without dividing by zero. Rendered for **k ≤ 8** because the image is 2ᵏ × 2ᵏ.

### Jaccard similarity

`J(A, B) = |A ∩ B| / |A ∪ B|` over the *sets* of k-mers (presence/absence).
Ranges 0 (nothing shared) to 1 (identical sets).

### Cosine similarity

`cos(A, B) = (A · B) / (‖A‖ · ‖B‖)` over the k-mer *frequency vectors* (counts).
Ranges 0 (orthogonal) to 1 (proportional composition). Captures similarity even
when absolute depths differ.

---

## Use it as a Python library

The `kmerlab` package is importable and has no Flask dependency:

```python
from kmerlab.sequence_parser import parse_file
from kmerlab.kmer_counter import count_kmers, reverse_complement, canonical
from kmerlab.metrics import summarize, compare_profiles
from kmerlab.exports import kmer_counts_to_csv

# Parse (FASTA/FASTQ, plain or .gz auto-detected)
parsed = parse_file("samples/example.fasta")

# Count k-mers
result = count_kmers(parsed.records, k=4, canonical_mode=True,
                     include_ambiguous=False)
print(result.counted_kmers, result.unique_kmers)
print(result.top(5))                      # [(kmer, count), ...]

# Full metric summary (JSON-serializable dict)
print(summarize(result, parsed.records, top_n=20))

# Building blocks
reverse_complement("AAT")                 # -> "ATT"
canonical("TTG")                          # -> "CAA"

# Compare two profiles (same k and canonical setting required)
a = count_kmers(parse_file("samples/compare_a.fasta").records, k=4)
b = count_kmers(parse_file("samples/compare_b.fasta").records, k=4)
print(compare_profiles(a, b).to_dict())

# Export
open("counts.csv", "w").write(kmer_counts_to_csv(result))
```

---

## HTTP API reference

All endpoints are local (`http://127.0.0.1:5001`) and accept
`multipart/form-data`. Provide input either as a file field or, for a bundled
sample, as `<field>_text` containing the file contents.

| Method | Path | Purpose |
|--------|------|---------|
| GET | `/` | Analyzer page |
| GET | `/compare` | Compare page |
| GET | `/about` | Method / about page |
| GET | `/api/samples` | JSON list of bundled sample filenames |
| GET | `/samples/<name>` | Download a bundled sample file |
| POST | `/api/analyze` | Analyze one file → summary + chart data URIs |
| POST | `/api/export/csv` | k-mer counts as CSV download |
| POST | `/api/export/json` | Summary as JSON download |
| POST | `/api/compare` | Compare two files → comparison + chart |
| POST | `/api/compare/export` | Comparison as JSON or CSV (`format=json|csv`) |

**`/api/analyze` fields:** `file` (or `file_text`), `k` (default 4),
`top_n` (default 20), `canonical` (bool), `include_ambiguous` (bool),
`fcgr` (bool). Returns `{ "summary": {...}, "charts": { "top_bar", "spectrum",
"fcgr" } }` where each chart is a base64 PNG data URI (`fcgr` is `null` for k > 8).

**`/api/compare` fields:** `file_a`/`file_a_text`, `file_b`/`file_b_text`, `k`,
`canonical`, `include_ambiguous`.

Errors return `{ "error": "..." }` with HTTP 400 (bad request) or 413 (upload
too large). Example:

```bash
curl -s -F "file=@samples/example.fasta" -F "k=4" -F "canonical=true" \
     http://127.0.0.1:5001/api/analyze | python -m json.tool
```

---

## Export formats

**k-mer counts CSV** (`/api/export/csv`) — one row per distinct k-mer, sorted by
count descending:

```
kmer,count,frequency
ACG,2,0.33333333
CGT,2,0.33333333
GTA,1,0.16666667
```

`frequency` is `count / counted_kmers`.

**Summary JSON** (`/api/export/json`) — the full metric summary: `k`,
`canonical`, `include_ambiguous`, `total_sequences`, `total_bases`,
`counted_kmers`, `unique_kmers`, `skipped_kmers`, `invalid_base_count`,
`gc_content`, `base_composition` (`counts` + `fractions`), `length_stats`,
`quality` (or `null` for FASTA), and `top_kmers`.

**Comparison JSON/CSV** (`/api/compare/export`) — `k`, `canonical`,
`shared_kmers`, `unique_to_a`, `unique_to_b`, `unique_kmers_a`,
`unique_kmers_b`, `jaccard_similarity`, `cosine_similarity`.

---

## Testing & benchmarks

```bash
pytest
```

The suite (76 tests) covers parsing, strict FASTQ error cases, k-mer counting,
IUPAC reverse complement, canonical + ambiguous modes, similarity metrics, the
k-mer spectrum, FCGR normalization and corner placement (including the
zero-count case), the gzip decompression limit, the FASTQ quality summary,
CSV/JSON exports, and the Flask API.

**Benchmark validation.** `benchmarks/` holds tiny, deterministic FASTA files
and hand-computed expected counts (`benchmarks/README.md` shows the worked
arithmetic). `tests/test_benchmarks.py` recomputes counts and asserts an exact
match against those numbers, so the counting logic is validated against values a
human verified by hand.

CI (`.github/workflows/tests.yml`) runs the suite on Python 3.10, 3.11, and 3.12.

---

## Folder structure

```
KmerLab/
├── app.py                    # Flask app: routes + API endpoints (port 5001)
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
- Similarity requires both files to use the same *k* and canonical setting
  (enforced by the app).
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
