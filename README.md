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

Install via conda:

```bash
conda create -n kingsmaker -c conda-forge -c bioconda \
    snakemake cutadapt pear fastqc \
    pandas pyyaml rapidfuzz umap-learn hdbscan python-Levenshtein \
    scikit-learn matplotlib seaborn openpyxl numpy
conda activate kingsmaker
```

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

Example:

```
Sample   R1                                       R2                                       Round  Fraction  Group     Experiment
3412244  /Volumes/Z/.../L008_R1_001.fastq.gz      /Volumes/Z/.../L008_R2_001.fastq.gz      1      All       Singlet   F11R
3412250  /Volumes/Z/.../L008_R1_001.fastq.gz      /Volumes/Z/.../L008_R2_001.fastq.gz      6      Top       Singlet   F11R
```

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

Built from per-sample RPM columns (`R{round}_{group}_{fraction}_RPM`):

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| `R1_Global_RPM` | Σ R1 RPM | Starting pool abundance |
| `R6_Total_Top_RPM` | Σ R6 Top RPM | Final-round fluorescent-fraction abundance |
| `R6_Total_Bot_RPM` | Σ R6 Bot RPM | Final-round non-fluorescent-fraction abundance |
| `Cumulative_Enrichment` | (R6_Top + 0.1) / (R1 + 0.1) | Fold-enrichment, start → finish |
| `Partition_Ratio_R6` | (R6_Top + 0.1) / (R6_Bot + 0.1) | Preference for fluorescent fraction at endpoint |
| `Aggregation_Index` | (Singlet_Top + 0.1) / (Doublet_Top + 0.1) | Tendency to behave monomerically |
| `Lineage` | Fuzzy match to known library scaffolds | `Lib_6`, `Lib_7a/b`, `Lib_9a/b`, etc. |

`Seq_ID` is an md5 hash of the sequence — stable and identical across every output file, so any two CSVs can be joined by `Seq_ID`.

## Notes

- **Counts are persisted; intermediates are not.** Re-running after `counts/<sample>_merged_counts.txt` exists skips trim/merge/QC for that sample. To force a re-run, delete the relevant count file (or `--forceall`).
- **The clustering stage is a Snakemake checkpoint.** The set of top-100 CSVs isn't known until the enrichment script runs — Snakemake re-evaluates the DAG afterwards and fans out clustering jobs accordingly.
- **Raw inputs are read-only.** Nothing in `samples.tsv` paths is modified.
