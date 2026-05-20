import yaml
import pandas as pd
import os

with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

metadata = pd.read_csv("input/samples.tsv", sep="\t").set_index("Sample", drop=False)
SAMPLES = [str(s) for s in metadata.index.tolist()]


def _stage_round(label):
    if 'Stage' not in metadata.columns:
        return None
    tagged = metadata.loc[metadata['Stage'].astype(str).str.strip() == label, 'Round']
    rounds = sorted({int(r) for r in tagged.unique()})
    if len(rounds) > 1:
        raise ValueError(f"Multiple rounds tagged as '{label}' in samples.tsv: {rounds}. Only one round can be {label}.")
    return rounds[0] if rounds else None


FIRST_ROUND = _stage_round('First')
LAST_ROUND = _stage_round('Last')


def get_all_clustering_outputs(wildcards):
    checkpoint_output = checkpoints.analyze_enrichment.get().output.top_dir
    files = [f for f in os.listdir(checkpoint_output) if f.endswith('.csv')]
    metrics = [f.replace(".csv", "") for f in files]
    return expand("results/clustering/plots/{m}_umap.png", m=metrics)


rule all:
    input:
        expand("qc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("merged/{sample}.assembled.fastq.gz", sample=SAMPLES),
        expand("counts/{sample}_merged_counts.txt", sample=SAMPLES),
        "results/selex_enrichment_report.csv",
        get_all_clustering_outputs,
        "results/metric_overlap_table.xlsx",
        "results/final_synthesis_order_list.csv"


rule trim:
    input:
        R1 = lambda wc: metadata.loc[wc.sample, "R1"].split(","),
        R2 = lambda wc: metadata.loc[wc.sample, "R2"].split(",")
    output:
        R1 = temp("trimmed/{sample}_R1.fastq.gz"),
        R2 = temp("trimmed/{sample}_R2.fastq.gz")
    threads: 8
    resources:
        disk_mb = 6000
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
            <(cat {input.R1}) <(cat {input.R2})
        """


rule merge_reads:
    input:
        R1 = "trimmed/{sample}_R1.fastq.gz",
        R2 = "trimmed/{sample}_R2.fastq.gz"
    output:
        assembled = temp("merged/{sample}.assembled.fastq.gz")
    threads: 8
    resources:
        disk_mb = 6000
    shell:
        """
        pear -f {input.R1} -r {input.R2} -o merged/{wildcards.sample} -j {threads}
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
        gunzip -c {input.merged} | awk 'NR%4==2' | sort | uniq -c | sort -nr > {output.counts}
        """


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
        top_dir = "results/top_100_lists"
    output:
        excel = "results/metric_overlap_table.xlsx"
    script:
        "scripts/generate_overlap_table.py"


def get_order_list_kings(wildcards):
    checkpoint_output = checkpoints.analyze_enrichment.get().output.top_dir
    available = set(os.listdir(checkpoint_output))

    candidates = []
    if LAST_ROUND is not None:
        for fraction in ('singlet', 'doublet'):
            fname = f"top_100_r{LAST_ROUND}_{fraction}_top_rpm.csv"
            if fname in available:
                candidates.append(f"results/clustering/kings/{fname[:-4]}_clustered_kings.csv")
    for always_try in ("top_100_cumulative_enrichment.csv", "master_top_selection_winners.csv"):
        if always_try in available:
            candidates.append(f"results/clustering/kings/{always_try[:-4]}_clustered_kings.csv")
    return candidates


rule generate_order_list:
    input:
        report = "results/selex_enrichment_report.csv",
        kings = get_order_list_kings
    output:
        order_list = "results/final_synthesis_order_list.csv"
    script:
        "scripts/compile_order_list.py"
