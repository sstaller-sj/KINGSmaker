# KINGSmaker: FACS-SELEX RNA Sensor Analysis Pipeline

A Snakemake pipeline for processing paired-end sequencing data from FACS-based SELEX experiments targeting fluorescent RNA sensor aptamers. It trims adapters, merges paired reads, quantifies unique-sequence abundance, ranks candidates across enrichment metrics, clusters them, and produces a final synthesis order list.

## Pipeline overview

Six stages, raw reads → synthesis candidates:

1. **Trim adapters** (`cutadapt`) — removes Illumina adapters from R1 and R2.
2. **Merge paired reads** (`PEAR`) — stitches R1/R2 overlap into single consensus reads.
3. **Count unique sequences** — collapses identical reads and tallies per sample.
4. **Enrichment analysis** — computes per-sample RPM, partition ratios, cumulative enrichment, aggregation indices, and lineage; ranks every sequence on every metric; emits top-100 lists.
5. **Cluster and pick "kings"** — UMAP + HDBSCAN per top-100 list; identifies four representative sequences per cluster (metric winner, enrichment winner, centroid king, consensus king).
6. **Aggregate** — metric-overlap matrix + final deduplicated synthesis order list.

## Repository layout

```
KINGSmaker/
├── Snakefile                # Workflow definition
├── config/
│   └── config.yaml          # Adapter/constant region sequences
├── input/
│   ├── samples.tsv          # Sample metadata + FASTQ paths
│   └── fastq/               # (optional) Local raw reads
├── scripts/
│   ├── analyze_selex_enrichment.py
│   ├── cluster_candidates.py
│   ├── compile_order_list.py
│   └── generate_overlap_table.py
├── trimmed/                 # Adapter-trimmed reads (auto-deleted after use)
├── merged/                  # PEAR-merged reads (auto-deleted after use)
├── counts/                  # Per-sample sequence counts (persisted)
├── qc/                      # FastQC reports
├── results/                 # Final outputs (see below)
└── logs/                    # Per-rule logs
```

## Requirements

Create the conda environment from the provided spec:

```bash
conda env create -f environment.yaml
conda activate kingsmaker
```

(Or with `mamba env create -f environment.yaml` if you have mamba — much faster.)

| Tool / library | Purpose |
|----------------|---------|
| snakemake | Workflow orchestration |
| cutadapt | Adapter trimming |
| PEAR | Paired-end read merging |
| FastQC | QC reports |
| pandas, numpy | Data wrangling |
| rapidfuzz, python-Levenshtein | Fuzzy matching for lineage + cluster consensus |
| umap-learn, hdbscan | Dimensionality reduction and clustering |
| scikit-learn | k-mer vectorization |
| matplotlib, seaborn | Plots |
| openpyxl | Excel output with conditional formatting |

Requires a POSIX shell (bash). Runs natively on macOS/Linux. On Windows, use WSL.

## Input setup

### `input/samples.tsv`

Tab-separated, one row per sample. Required columns:

| Column | Description |
|--------|-------------|
| `Sample` | Unique sample ID. Used in all downstream filenames. |
| `R1` | Path to R1 FASTQ.gz. Multiple lane files can be comma-separated. |
| `R2` | Path to R2 FASTQ.gz. Multiple lane files can be comma-separated. |
| `Round` | SELEX round number (e.g. `1`, `6`). Drives RPM column naming. |
| `Fraction` | Sort fraction — typically `Top`, `Bot`, or `All` (for unsorted inputs). |
| `Group` | Sample subgroup — typically `Singlet` or `Doublet`. |
| `Experiment` | Free-text experiment tag (provenance only). |
| `Stage` | `First`, `Last`, or blank. Tags which round serves as the input baseline and which as the endpoint for enrichment metrics. Only one round can be `First`, only one `Last`. Both are optional — metrics that need them are filled with `x` if missing. |

Example:

```
Sample   R1                                       R2                                       Round  Fraction  Group     Experiment  Stage
3412244  /Volumes/Z/.../L008_R1_001.fastq.gz      /Volumes/Z/.../L008_R2_001.fastq.gz      1      All       Singlet   F11R        First
3412250  /Volumes/Z/.../L008_R1_001.fastq.gz      /Volumes/Z/.../L008_R2_001.fastq.gz      6      Top       Singlet   F11R        Last
```

#### How `Stage` drives the enrichment math

- **`Cumulative_Enrichment`** = `R{Last}_Top / R{First}_Global` — requires both `First` and `Last` to be tagged, plus `Top` samples in the last round.
- **`Partition_Ratio_R{Last}`** = `R{Last}_Top / R{Last}_Bot` — requires `Last` plus both `Top` and `Bot` samples in that round.
- **`Aggregation_Index`** = `R{Last}_Singlet_Top / R{Last}_Doublet_Top` — requires `Last` plus both groups with `Top` samples.

When a metric can't be computed (the data isn't there), its column still appears in the master report but every cell is literally `x`, and no top-100 list is generated for it. This makes it obvious in the output which analyses ran and which were skipped.

#### Multi-lane samples

If a sample was sequenced across multiple lanes (e.g. `L001` + `L002`), list all lane files comma-separated within the same column:

```
3412300  /Volumes/Z/.../L001_R1.fastq.gz,/Volumes/Z/.../L002_R1.fastq.gz   /Volumes/Z/.../L001_R2.fastq.gz,/Volumes/Z/.../L002_R2.fastq.gz   ...
```

The trim rule concatenates lane streams during cutadapt — no pre-concatenation step needed. Single-lane samples just have one path per column.

### `config/config.yaml`

Defines the library's adapter and constant-region sequences:

| Key | Meaning |
|-----|---------|
| `five_prime_adapter_R{1,2}`, `three_prime_adapter_R{1,2}` | Illumina adapter sequences trimmed from each read direction. |
| `*_constant_*` | Library constant regions (informational only; not currently consumed by the trim rule). |

## Running the pipeline

All commands run from the repo root, with the conda env activated.

### Dry run (check the DAG)
```bash
snakemake -n
```

### Single-sample smoke test
```bash
snakemake --cores 4 --resources disk_mb=12000 counts/<sample_id>_merged_counts.txt
```

### Full run
```bash
snakemake --cores 8 --resources disk_mb=12000 --keep-going --retries 3
```

### Disk-constrained / streaming-from-network setup

When raw FASTQs live on a mounted network share (or any remote-feel filesystem) and local disk is tight:

- The pipeline marks `trimmed/` and `merged/` outputs as `temp()`. Snakemake deletes each intermediate the moment its downstream consumer finishes, so only `counts/`, `qc/`, `results/`, and `logs/` persist.
- `--resources disk_mb=N` caps how many samples Snakemake processes concurrently. Each in-flight sample reserves `disk_mb = 6000` (≈6 GB). A budget of `12000` allows ~2 samples in flight.
- `--retries 3` re-attempts a rule if the share hiccups.

Tune `disk_mb` once you've measured a real sample's peak intermediate size.

## Outputs

| Path | Description |
|------|-------------|
| `counts/<sample>_merged_counts.txt` | Per-sample unique-sequence counts (count, sequence). Persisted. |
| `qc/<sample>_R{1,2}_fastqc.html` | FastQC reports on trimmed reads. |
| `results/selex_enrichment_report.csv` | Master table: every unique sequence × every metric, ranked. |
| `results/top_100_lists/top_100_<metric>.csv` | Top 100 sequences per metric. |
| `results/top_100_lists/master_top_*_winners.csv` | Union pools across metrics. |
| `results/clustering/data/<metric>_clustered.csv` | All top-100 sequences with cluster IDs and UMAP coords. |
| `results/clustering/plots/<metric>_umap.png` | Annotated UMAP scatter per metric. |
| `results/clustering/kings/<metric>_clustered_kings.csv` | Four representative sequences per cluster. |
| `results/metric_overlap_table.xlsx` | Two-sheet matrix: which sequences won which metrics. |
| `results/final_synthesis_order_list.csv` | Deduplicated candidate list across king files — what you'd actually order. |

## Computed metrics

Built from per-sample RPM columns (`R{round}_{group}_{fraction}_RPM`). Round numbers in metric names come from the `Stage` tags — `R{First}` and `R{Last}` substitute the round you marked in `samples.tsv`.

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| `R{First}_Global_RPM` | Σ RPM across all First-round samples | Starting pool abundance |
| `R{Last}_Total_Top_RPM` | Σ RPM across Last-round Top samples | Endpoint fluorescent-fraction abundance |
| `R{Last}_Total_Bot_RPM` | Σ RPM across Last-round Bot samples | Endpoint non-fluorescent-fraction abundance |
| `Cumulative_Enrichment` | (R{Last}_Top + 0.1) / (R{First}_Global + 0.1) | Fold-enrichment, start → finish |
| `Partition_Ratio_R{Last}` | (R{Last}_Top + 0.1) / (R{Last}_Bot + 0.1) | Preference for fluorescent fraction at endpoint |
| `Aggregation_Index` | (R{Last}_Singlet_Top + 0.1) / (R{Last}_Doublet_Top + 0.1) | Tendency to behave monomerically |
| `Lineage` | Fuzzy match to known library scaffolds | `Lib_6`, `Lib_7a/b`, `Lib_9a/b`, etc. |

Any metric whose inputs aren't present in `samples.tsv` is filled with the literal string `x` in the master report (no top-100 list, no clustering).

`Seq_ID` is an md5 hash of the sequence — stable and identical across every output file, so any two CSVs can be joined by `Seq_ID`.

## Notes

- **Counts are persisted; intermediates are not.** Re-running after `counts/<sample>_merged_counts.txt` exists skips trim/merge/QC for that sample. To force a re-run, delete the relevant count file (or `--forceall`).
- **The clustering stage is a Snakemake checkpoint.** The set of top-100 CSVs isn't known until the enrichment script runs — Snakemake re-evaluates the DAG afterwards and fans out clustering jobs accordingly.
- **Raw inputs are read-only.** Nothing in `samples.tsv` paths is modified.
