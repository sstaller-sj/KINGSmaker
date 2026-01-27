import yaml
import pandas as pd
import os

# 1. Configuration and Metadata
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

metadata = pd.read_csv("input/samples.tsv", sep="\t").set_index("Sample", drop=False)
SAMPLES = [str(s) for s in metadata.index.tolist()]

# 2. Dynamic Input Function
# This function waits for the analyze_enrichment checkpoint to finish,
# then discovers every CSV in the results folder to request UMAPs.
def get_all_clustering_outputs(wildcards):
    # This line creates the dependency on the checkpoint
    checkpoint_output = checkpoints.analyze_enrichment.get().output.top_dir
    
    # List all CSVs in the results/top_100_lists folder
    files = [f for f in os.listdir(checkpoint_output) if f.endswith('.csv')]
    metrics = [f.replace(".csv", "") for f in files]
    
    # Return the expected plot paths for every file found
    return expand("results/clustering/plots/{m}_umap.png", m=metrics)


rule all:
    input:
        # QC reports
        expand("qc/{sample}_R1_fastqc.html", sample=SAMPLES),
        # Merged/Assembled reads (this triggers PEAR)
        expand("merged/{sample}.assembled.fastq.gz", sample=SAMPLES),
        # Final merged counts
        expand("counts/{sample}_merged_counts.txt", sample=SAMPLES),
        # Enrichment report
        "results/selex_enrichment_report.csv",
        get_all_clustering_outputs, # Trigger the dynamic discovery
        "results/metric_overlap_table.xlsx",
        "results/final_synthesis_order_list.csv"
rule trim:
    input:
        R1 = lambda wc: SAMPLES[wc.sample]["R1"],
        R2 = lambda wc: SAMPLES[wc.sample]["R2"]
    output:
        R1 = "trimmed/{sample}_R1.fastq.gz",
        R2 = "trimmed/{sample}_R2.fastq.gz"
    threads: 8
    params:
        constant_5p_R1 = config["five_prime_constant_R1"],
        constant_3p_R1 = config["three_prime_constant_R1"],
        adapter_5p_R1 = config["five_prime_adapter_R1"],
        adapter_3p_R1 = config["three_prime_adapter_R1"],

        constant_5p_R2 = config["five_prime_constant_R2"],
        constant_3p_R2 = config["three_prime_constant_R2"],
        adapter_5p_R2 = config["five_prime_adapter_R2"],
        adapter_3p_R2 = config["three_prime_adapter_R2"],
        minlen = 45
    shell:
        """
        mkdir -p trimmed
        cutadapt \
            -g {params.adapter_5p_R1} \
            -G {params.adapter_5p_R2} \
            -a {params.adapter_3p_R1} \
            -A {params.adapter_3p_R2} \
            --overlap 10 \
            --error-rate 0.1 \
            --trim-n \
            -m {params.minlen} \
            -o {output.R1} -p {output.R2} \
            {input.R1} {input.R2}
        """

rule merge_reads:
    input:
        R1 = "trimmed/{sample}_R1.fastq.gz",
        R2 = "trimmed/{sample}_R2.fastq.gz"
    output:
        assembled = "merged/{sample}.assembled.fastq.gz"
    threads: 8
    shell:
        """
        # PEAR produces several files, we'll use a temp prefix
        pear -f {input.R1} -r {input.R2} -o merged/{wildcards.sample} -j {threads}
        
        # Gzip the main result and cleanup if necessary
        gzip -c merged/{wildcards.sample}.assembled.fastq > {output.assembled}
        rm merged/{wildcards.sample}.assembled.fastq
        """

rule fastqc:
    input:
        R1 = "trimmed/{sample}_R1.fastq.gz",
        R2 = "trimmed/{sample}_R2.fastq.gz"
    output:
        R1_html = "qc/{sample}_R1_fastqc.html",
        R2_html = "qc/{sample}_R2_fastqc.html"
    threads: 2
    shell:
        """
        mkdir -p qc
        fastqc -o qc {input.R1} {input.R2}
        """

rule count_sequences:
    input:
        merged = "merged/{sample}.assembled.fastq.gz"
    output:
        counts = "counts/{sample}_merged_counts.txt"
    shell:
        """
        mkdir -p counts
        # Count only the merged consensus sequences
        gunzip -c {input.merged} | awk 'NR%4==2' | sort | uniq -c | sort -nr > {output.counts}
        """

# Changed to 'checkpoint' to allow dynamic downstream rules
checkpoint analyze_enrichment:
    input:
        counts = expand("counts/{sample}_merged_counts.txt", sample=SAMPLES),
        metadata = "input/samples.tsv"
    output:
        report = "results/selex_enrichment_report.csv",
        top_dir = directory("results/top_100_lists")
    script:
        "scripts/analyze_selex_enrichment.py"

rule cluster_any_csv:
    input:
        csv = "results/top_100_lists/{metric}.csv"
    output:
        clustered_csv = "results/clustering/data/{metric}_clustered.csv",
        umap_plot = "results/clustering/plots/{metric}_umap.png",
        kings_csv = "results/clustering/kings/{metric}_clustered_kings.csv"
    script:
        "scripts/cluster_candidates.py"

rule generate_overlap_excel:
    input:
        # This ensures the checkpoint has run and the folder is full
        top_dir = "results/top_100_lists"
    output:
        excel = "results/metric_overlap_table.xlsx"
    script:
        "scripts/generate_overlap_table.py"

rule generate_order_list:
    input:
        report = "results/selex_enrichment_report.csv",
        kings = [
            "results/clustering/kings/top_100_r6_singlet_top_rpm_clustered_kings.csv",
            "results/clustering/kings/top_100_r6_doublet_top_rpm_clustered_kings.csv",
            "results/clustering/kings/top_100_cumulative_enrichment_clustered_kings.csv",
            "results/clustering/kings/master_top_selection_winners_clustered_kings.csv"
        ]
    output:
        order_list = "results/final_synthesis_order_list.csv"
    script:
        "scripts/compile_order_list.py"