import hashlib
import os
import re

import pandas as pd
from rapidfuzz import fuzz

input_metadata = snakemake.input.metadata
output_file = snakemake.output.report

X = "x"  # marker placed in metric columns that can't be computed for the current Stage configuration


def identify_lineage(seq):
    seq = seq.replace('U', 'T').upper()
    ref_9a = "GACGACACCAACCGGCACUGGCGCCGGAGGGCUUCCUCCCAAUCGUGGCGUGUCGGCGGGUCUGUAAUGGAUGUCGUC".replace('U', 'T')
    ref_9b = "GAUCGCACUGGCGCCGGACACCAACCGGUAUUUCGGGUCUGUAAUGGAUGUCCCAAUCGUGGCGUGUCGGCGAUC".replace('U', 'T')
    ref_6 = "CTGCTTCGGCAG"

    if fuzz.partial_ratio(ref_9a, seq) > 80: return "Lib_9a"
    if fuzz.partial_ratio(ref_9b, seq) > 80: return "Lib_9b"
    if fuzz.partial_ratio(ref_6, seq) > 80: return "Lib_6"
    if re.search(r"TGGCGCC.{7,11}CGTGGCGTGT", seq): return "Lib_7b"
    if re.search(r"TGGC.{18,22}GTG", seq): return "Lib_7a"
    if 70 <= len(seq) <= 110: return "Lib_1/ Unknown"
    return "Unknown_Artifact"


def consistent_id(seq):
    return "seq_" + hashlib.md5(seq.encode()).hexdigest()[:8]


def resolve_stage_round(meta, stage_label):
    """Return the unique round tagged with stage_label, or None. Raise if more than one round is tagged."""
    if 'Stage' not in meta.columns:
        raise ValueError("samples.tsv is missing required 'Stage' column (values: First, Last, or blank).")
    tagged = meta.loc[meta['Stage'].astype(str).str.strip() == stage_label, 'Round']
    rounds = sorted({int(r) for r in tagged.unique()})
    if len(rounds) > 1:
        raise ValueError(f"Multiple rounds tagged as '{stage_label}': {rounds}. Only one round can be {stage_label}.")
    return rounds[0] if rounds else None


meta = pd.read_csv(input_metadata, sep="\t").set_index("Sample")
first_round = resolve_stage_round(meta, 'First')
last_round = resolve_stage_round(meta, 'Last')

all_dfs = []
for sample_id in meta.index:
    path = f"counts/{sample_id}_merged_counts.txt"
    if not os.path.exists(path):
        continue
    df = pd.read_csv(path, sep=r'\s+', names=['count', 'sequence'], header=None)
    if df.empty:
        continue
    total = df['count'].sum()
    df['RPM'] = (df['count'] / total) * 1e6
    row = meta.loc[sample_id]
    col_name = f"R{row['Round']}_{row['Group']}_{row['Fraction']}_RPM"
    df = df[['sequence', 'RPM']].set_index('sequence')
    df.columns = [col_name]
    all_dfs.append(df)

if not all_dfs:
    raise RuntimeError("No count files found. Run trim/merge/count first.")

master = pd.concat(all_dfs, axis=1).fillna(0)

# Stage-driven aggregates
first_prefix = f"R{first_round}_" if first_round is not None else None
last_prefix = f"R{last_round}_" if last_round is not None else None

first_global_col = None
if first_round is not None:
    first_cols = [c for c in master.columns if c.startswith(first_prefix)]
    first_global_col = f'R{first_round}_Global_RPM'
    master[first_global_col] = master[first_cols].sum(axis=1)

last_top_cols, last_bot_cols = [], []
last_top_col_name = last_bot_col_name = None
if last_round is not None:
    last_top_cols = [c for c in master.columns if c.startswith(last_prefix) and 'Top' in c]
    last_bot_cols = [c for c in master.columns if c.startswith(last_prefix) and 'Bot' in c]
    if last_top_cols:
        last_top_col_name = f'R{last_round}_Total_Top_RPM'
        master[last_top_col_name] = master[last_top_cols].sum(axis=1)
    if last_bot_cols:
        last_bot_col_name = f'R{last_round}_Total_Bot_RPM'
        master[last_bot_col_name] = master[last_bot_cols].sum(axis=1)
    if last_top_cols and last_bot_cols:
        master[f'R{last_round}_Global_RPM'] = master[last_top_col_name] + master[last_bot_col_name]

# Cumulative enrichment: First-round baseline → Last-round Top
if first_global_col and last_top_col_name:
    master['Cumulative_Enrichment'] = (master[last_top_col_name] + 0.1) / (master[first_global_col] + 0.1)
else:
    master['Cumulative_Enrichment'] = X

# Partition ratio: Last-round Top vs Bot
pr_col = f'Partition_Ratio_R{last_round}' if last_round is not None else 'Partition_Ratio'
if last_top_col_name and last_bot_col_name:
    master[pr_col] = (master[last_top_col_name] + 0.1) / (master[last_bot_col_name] + 0.1)
else:
    master[pr_col] = X

# Aggregation index: Singlet_Top vs Doublet_Top at Last round
if last_round is not None:
    s_cols = [c for c in master.columns if c.startswith(last_prefix) and 'Singlet' in c and 'Top' in c]
    d_cols = [c for c in master.columns if c.startswith(last_prefix) and 'Doublet' in c and 'Top' in c]
else:
    s_cols, d_cols = [], []
if s_cols and d_cols:
    master['Aggregation_Index'] = (master[s_cols[0]] + 0.1) / (master[d_cols[0]] + 0.1)
else:
    master['Aggregation_Index'] = X

master['Grand_Total_RPM'] = master.select_dtypes(include='number').filter(like='_RPM').sum(axis=1)
master['Sequence_Length'] = master.index.str.len()

master.index.name = 'Sequence'
master = master.reset_index()
master['Lineage'] = master['Sequence'].map(identify_lineage)
master['Seq_ID'] = master['Sequence'].apply(consistent_id)

# Sort by the strongest available signal
sort_key = 'Cumulative_Enrichment' if pd.api.types.is_numeric_dtype(master['Cumulative_Enrichment']) else 'Grand_Total_RPM'
master = master.sort_values(by=sort_key, ascending=False)

# Rank only numeric metrics
info_cols = {'Seq_ID', 'Sequence', 'Sequence_Length', 'Lineage'}
metric_cols = [c for c in master.columns
               if c not in info_cols
               and not c.endswith('_Rank')
               and pd.api.types.is_numeric_dtype(master[c])]
x_cols = [c for c in master.columns
          if c not in info_cols
          and not c.endswith('_Rank')
          and c not in metric_cols]

for metric in metric_cols:
    master[f'{metric}_Rank'] = master[metric].rank(ascending=False, method='min').astype(int)

rank_cols = [c for c in master.columns if c.endswith('_Rank')]
ordered = ['Seq_ID', 'Sequence', 'Lineage', 'Sequence_Length'] + metric_cols + x_cols + rank_cols
master = master[ordered]

# Per-metric top-100 lists
top_lists_dir = os.path.join(os.path.dirname(output_file), "top_100_lists")
os.makedirs(top_lists_dir, exist_ok=True)

all_winner_indices = set()
selected_winner_indices = set()

for metric in metric_cols:
    top_df = master.sort_values(by=metric, ascending=False).head(100)
    all_winner_indices.update(top_df.index)

    # Selection winners exclude First-round columns and Grand_Total — these track abundance, not enrichment
    is_first_round_metric = first_prefix is not None and metric.startswith(first_prefix)
    is_grand_total = 'Grand_Total' in metric
    if not (is_first_round_metric or is_grand_total):
        selected_winner_indices.update(top_df.index)

    top_df.to_csv(os.path.join(top_lists_dir, f"top_100_{metric.lower()}.csv"), index=False)

master_all = master.loc[list(all_winner_indices)].sort_values(by=sort_key, ascending=False)
master_all.to_csv(os.path.join(top_lists_dir, "master_top_all_winners.csv"), index=False)

master_selected = master.loc[list(selected_winner_indices)].sort_values(by=sort_key, ascending=False)
master_selected.to_csv(os.path.join(top_lists_dir, "master_top_selection_winners.csv"), index=False)

master.to_csv(output_file, index=False)
