import os

# Limit BLAS threads before numpy/UMAP/HDBSCAN import — prevents deadlocks when these libraries share the threadpool
os.environ.update({
    'OPENBLAS_NUM_THREADS': '1',
    'MKL_NUM_THREADS': '1',
    'OMP_NUM_THREADS': '1'
})

import pandas as pd
import numpy as np
import umap
import hdbscan
from sklearn.feature_extraction.text import CountVectorizer
import matplotlib.pyplot as plt
import seaborn as sns
from Levenshtein import distance

input_path = snakemake.input.csv
df = pd.read_csv(input_path)

if 'Seq_ID' not in df.columns:
    raise KeyError(f"Seq_ID column not found in {input_path}. Check analyze_selex_enrichment.py")

base_name = os.path.basename(input_path).replace(".csv", "").replace("top_100_", "")

if len(df) < 5:
    df['Cluster_ID'] = -1
    df['UMAP_1'] = 0
    df['UMAP_2'] = 0
    df.to_csv(snakemake.output.clustered_csv, index=False)
    os.makedirs(os.path.dirname(snakemake.output.kings_csv), exist_ok=True)
    pd.DataFrame(columns=['Cluster_ID', 'King_Type', 'Seq_ID', 'Sequence']).to_csv(snakemake.output.kings_csv, index=False)
    plt.figure().savefig(snakemake.output.umap_plot)
else:
    def get_kmers(seq, k=3):
        return [seq[i:i+k] for i in range(len(seq) - k + 1)]

    vectorizer = CountVectorizer(analyzer=lambda x: get_kmers(x, k=3))
    x_counts = vectorizer.fit_transform(df['Sequence'])

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
    embedding = reducer.fit_transform(x_counts)

    df['UMAP_1'] = embedding[:, 0]
    df['UMAP_2'] = embedding[:, 1]

    clusterer = hdbscan.HDBSCAN(min_cluster_size=3, gen_min_span_tree=True)
    df['Cluster_ID'] = clusterer.fit_predict(embedding)

    df.to_csv(snakemake.output.clustered_csv, index=False)

    kings = []
    enrich_col = 'Cumulative_Enrichment' if 'Cumulative_Enrichment' in df.columns else df.columns[1]

    for cid in df['Cluster_ID'].unique():
        if cid == -1: continue

        cluster_subset = df[df['Cluster_ID'] == cid].copy()

        metric_king = cluster_subset.iloc[0]
        enrich_king = cluster_subset.sort_values(by=enrich_col, ascending=False).iloc[0]

        c1, c2 = cluster_subset['UMAP_1'].mean(), cluster_subset['UMAP_2'].mean()
        cluster_subset['dist'] = np.sqrt((cluster_subset['UMAP_1']-c1)**2 + (cluster_subset['UMAP_2']-c2)**2)
        centroid_king = cluster_subset.sort_values(by='dist').iloc[0]

        sample_seqs = cluster_subset['Sequence'].tolist()
        def avg_lev(s): return sum(distance(s, x) for x in sample_seqs) / len(sample_seqs)
        cluster_subset['lev_score'] = cluster_subset['Sequence'].apply(avg_lev)
        consensus_king = cluster_subset.sort_values(by='lev_score').iloc[0]

        for k_type, k_data in [('Metric_Winner', metric_king),
                               ('Enrichment_Winner', enrich_king),
                               ('Centroid_King', centroid_king),
                               ('Consensus_King', consensus_king)]:
            kings.append({
                'Cluster_ID': cid,
                'King_Type': k_type,
                'Seq_ID': k_data['Seq_ID'],
                'Sequence': k_data['Sequence'],
                'Enrichment': k_data.get('Cumulative_Enrichment', 'N/A'),
                'Lineage': k_data.get('Lineage', 'N/A')
            })

    os.makedirs(os.path.dirname(snakemake.output.kings_csv), exist_ok=True)
    pd.DataFrame(kings).to_csv(snakemake.output.kings_csv, index=False)

    plt.figure(figsize=(12, 10))
    hue_col = 'Lineage' if 'Lineage' in df.columns else 'Cluster_ID'
    size_col = 'Cumulative_Enrichment' if 'Cumulative_Enrichment' in df.columns else None
    size_data = np.log10(df[size_col] + 1) if size_col else None

    scatter = sns.scatterplot(
        x='UMAP_1', y='UMAP_2', hue=hue_col,
        size=size_data, sizes=(40, 400), palette='Set1',
        data=df, alpha=0.6, edgecolor='w', linewidth=0.5
    )

    for cid in df['Cluster_ID'].unique():
        if cid == -1: continue
        c_df = df[df['Cluster_ID'] == cid]
        label = f"C{cid}"
        if 'Lineage' in c_df.columns:
            top_lin = c_df['Lineage'].mode()[0]
            label = f"C{cid}\n({top_lin})"

        plt.text(
            c_df['UMAP_1'].mean(), c_df['UMAP_2'].mean(), label,
            fontsize=9, fontweight='bold', ha='center',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2')
        )

    plt.title(f"UMAP Clustering: {base_name.replace('_', ' ').upper()}")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Groups")
    plt.tight_layout()
    plt.savefig(snakemake.output.umap_plot, dpi=300)
    plt.close()
