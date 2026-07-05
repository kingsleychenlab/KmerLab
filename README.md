# KmerLab

**Local k-mer analysis dashboard for FASTA/FASTQ files.**

KmerLab is a small, self-contained web app for exploring the k-mer composition
of DNA sequence files. Upload (or select) a FASTA/FASTQ file, choose a value of
*k*, and get counts, summary statistics, visualizations, comparisons, and
exports — all computed **locally on your own machine**.

> **No cloud. No database. No login. No API keys. No internet required.**
> Your sequence data never leaves your computer. Uploaded files are processed in
> memory and are never written to disk.

---

## The problem it solves

k-mer counting is a foundational step in genomics — used for genome assembly,
error and contamination detection, taxonomic classification, and comparing
datasets. Serious tools like [KMC](https://github.com/refresh-bio/KMC) and
[Jellyfish](https://github.com/gmarcais/Jellyfish) are fast but command-line
only and can be intimidating for students. Web tools usually mean uploading your
data to someone else's server.

KmerLab fills the gap: a friendly, visual, **completely local** utility that a
biology or bioinformatics student can run in a minute and actually use for
real (small-to-medium) files — without giving up privacy or fighting a CLI.

---

## Features

- **Input**: FASTA and FASTQ, with automatic format detection.
- **Extensions**: `.fa`, `.fasta`, `.fq`, `.fastq`, plus gzip-compressed `.gz`.
- **k-mer counting** with a sliding window and `collections.Counter`.
- **Canonical k-mers** — merge each k-mer with its reverse complement.
- **Ambiguous bases** — include or exclude k-mers containing `N` (or other
  non-ACGT characters); invalid bases are reported cleanly, never silently.
- **Summary stats**: total sequences, total bases, valid k-mers, unique k-mers,
  GC content, skipped k-mers, invalid bases, sequence-length stats.
- **Base composition** bar (A/C/G/T/other) using the sequencing-chromatogram
  nucleotide colour system that runs through the whole UI.
- **Visualizations**: top-k-mer bar chart, k-mer frequency spectrum histogram,
  and an optional **FCGR** (Frequency Chaos Game Representation) heatmap.
- **Compare two files**: shared k-mers, unique-to-A, unique-to-B, **Jaccard**
  similarity, and **cosine** similarity, with a side-by-side bar chart.
- **Exports**: k-mer counts to **CSV**, analysis summary to **JSON**, and
  comparison results to JSON or CSV.
- **Clean error handling**: invalid files, invalid *k*, empty sequences,
  malformed FASTQ, unsupported formats, and oversized uploads all produce clear
  messages instead of stack traces.
- **Bundled sample files** for instant demoing.

---

## Tech stack

- **Python 3.9+**
- **Flask** — web framework
- **pandas** — CSV/table export
- **matplotlib** — charts (headless `Agg` backend)
- **HTML / CSS / vanilla JavaScript** — no React, no Tailwind, no Bootstrap
- **pytest** — tests

No database, no ORM, no external APIs, no build step.

---

## Screenshots

_Add screenshots here once you run it locally:_

- `docs/analyzer.png` — the Analyzer page with summary cards and charts
- `docs/compare.png` — the Compare page with similarity scores
- `docs/fcgr.png` — an FCGR heatmap

---

## Installation

```bash
git clone <your-fork-url> KmerLab
cd KmerLab

python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install -r requirements.txt
```

## How to run locally

```bash
python app.py
```

Then open **http://127.0.0.1:5000** in your browser. The app binds to
`127.0.0.1` (localhost) only.

## How to run the tests

```bash
pytest
```

(48 tests covering parsing, counting, canonical k-mers, reverse complement,
invalid-base handling, similarity metrics, exports, the compare workflow, and
the Flask API.)

---

## Example usage

1. Open the **Analyzer** page.
2. Click **load a sample → `example.fasta`** (or upload your own file).
3. Set **k = 4**, optionally tick **Canonical k-mers**.
4. Click **Analyze**.
5. Inspect the summary cards, top-k-mer bar chart, frequency spectrum, FCGR
   heatmap, and the frequency table.
6. Click **Export k-mer counts (CSV)** or **Export summary (JSON)**.

To compare two datasets, open the **Compare** page, load `compare_a.fasta` and
`compare_b.fasta`, choose *k*, and click **Compare** to see Jaccard/cosine
similarity and shared/unique k-mer counts.

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

- **k-mer**: a length-*k* substring extracted with a one-base sliding window.
- **Reverse complement**: A↔T, C↔G, reversed.
- **Canonical k-mer**: the lexicographically smaller of a k-mer and its reverse
  complement — merges the two DNA strand orientations.
- **GC content**: fraction of G/C among A/C/G/T bases only.
- **Jaccard similarity**: `|A ∩ B| / |A ∪ B|` over the *sets* of k-mers.
- **Cosine similarity**: `(A · B) / (‖A‖ · ‖B‖)` over k-mer *frequency vectors*.
- **FCGR**: k-mer frequencies mapped onto a 2ᵏ × 2ᵏ grid (corners A/C/G/T).

See the in-app **About / Method** page for fuller explanations.

---

## Folder structure

```
KmerLab/
├── app.py                    # Flask app: routes + API endpoints
├── requirements.txt
├── pytest.ini
├── README.md
├── kmerlab/                  # core library (importable, no Flask dependency)
│   ├── __init__.py
│   ├── sequence_parser.py    # FASTA/FASTQ parsing, format detection, gzip
│   ├── kmer_counter.py       # counting, reverse complement, canonical, validate_k
│   ├── metrics.py            # GC, summaries, Jaccard, cosine, comparison
│   ├── visualizations.py     # matplotlib charts -> base64 PNG data URIs
│   └── exports.py            # CSV / JSON export helpers (pandas)
├── templates/                # Jinja2 HTML (base, index, compare, about)
├── static/                   # style.css, main.js (vanilla JS)
├── samples/                  # example.fasta/.fastq, invalid.fastq, compare_a/b.fasta
└── tests/                    # pytest suite
```

---

## Limitations

- Files are processed **in memory**; very large genomes may exceed available
  RAM. The web UI caps uploads at **50 MB**.
- **Exact** counting is used (no probabilistic sketching), so memory scales with
  the number of distinct k-mers and with *k*. The UI caps *k* at 31.
- **Nucleotide sequences only** — no protein k-mer support yet.
- The FCGR heatmap renders only for **k ≤ 8** (a 2ᵏ × 2ᵏ image).
- Similarity metrics are only meaningful when both files use the same *k* and
  canonical setting (the app enforces this).
- `debug=True` is on for local development convenience — turn it off if you ever
  expose the app beyond localhost.

---

## Future improvements

- Larger-file **streaming** counting mode (bounded memory).
- Broader compressed-file support (bzip2, zstd).
- **Protein** sequence / k-mer support.
- More sequence **QC** metrics (per-position quality, length distributions).
- **MinHash** sketching for fast approximate comparison of large datasets.
- A standalone **command-line interface**.
- **Docker** support.
- **PDF** report export.

---

## License

Open source — released for educational and research use. Add your preferred
license (e.g. MIT) here.
