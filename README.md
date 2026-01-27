# FACS-SELEX RNA Sensor Analysis Pipeline

This repository contains a Snakemake-based analysis pipeline for processing
high-throughput sequencing data from a FACS-based SELEX experiment targeting
fluorescent RNA sensor aptamers.

The pipeline trims constant regions and adapters, performs quality control,
and quantifies sequence abundances for paired-end sequencing data (R1/R2).

---

## Project Structure

# FACS-SELEX RNA Sensor Analysis Pipeline

This repository contains a Snakemake-based analysis pipeline for processing
high-throughput sequencing data from a FACS-based SELEX experiment targeting
fluorescent RNA sensor aptamers.

The pipeline trims constant regions and adapters, performs quality control,
and quantifies sequence abundances for paired-end sequencing data (R1/R2).

---

## Project Structure

SELEX/
├── README.md
├── Snakefile # Snakemake workflow
├── config/ # Pipeline configuration files
│ └── config.yaml
├── input/
│ ├── fastq/ # Raw sequencing reads (immutable)
│ └── samples.tsv # Sample table
├── trimmed/ # Trimmed FASTQ files
├── qc/ # FastQC reports
├── counts/ # Sequence abundance tables
├── logs/ # Snakemake and tool logs
├── metadata/ # Run provenance and reports
├── archive/ # Frozen snapshots of completed runs
├── .snakemake/ # Snakemake internals (auto-generated)
├── .vscode/ # Editor configuration (optional)
└── .github/ # GitHub configuration (optional)


---

## Requirements

- Snakemake (≥ 7.x)
- Python (≥ 3.10)
- cutadapt
- fastqc
- Standard UNIX tools (`awk`, `sort`, `uniq`, `gzip`)

All tools are expected to be available in the active conda environment.

---

## Input Files

### `input/samples.tsv`

A tab-separated file describing samples and their paired-end sequencing files.

Example:

sample R1 R2
sample1 input/fastq/sample1_R1.fastq.gz input/fastq/sample1_R2.fastq.gz

### Raw sequencing reads

All raw FASTQ files should be placed in:

input/fastq


These files are treated as immutable inputs and are never modified by the
pipeline.

---

## Configuration

Pipeline parameters such as adapter sequences, constant regions, trimming
behavior, and quality thresholds are defined in:

config/config.yaml


This configuration file is actively edited between runs. A snapshot of the
configuration and sample table used for each run is copied into `metadata/`
to ensure reproducibility.

---

## Running the Pipeline

Perform a dry run to verify the workflow:

```bash
snakemake -n

Run the pipeline using 8 CPU cores:

snakemake --cores 8

snakemake --report metadata/snakemake_report.html

snakemake --report metadata/snakemake_report.html


Outputs
trimmed/ — trimmed paired-end FASTQ files
qc/ — FastQC HTML reports
counts/ — per-sample sequence abundance tables
logs/ — logs for each rule and sample
metadata/ — run configuration snapshots and reports

Archiving a Run
Once an analysis run is complete and finalized, results can be frozen into an
archive for long-term storage and reproducibility:

snakemake archive

Archived runs are stored in:

archive/

Archived data should not be modified.

Notes
This pipeline is designed for SELEX-style libraries consisting of a variable
region flanked by constant regions and sequenced using paired-end Illumina
technology. It is optimized for FACS-based selection of fluorescent RNA
sensor aptamers but can be adapted for related amplicon-based workflows.