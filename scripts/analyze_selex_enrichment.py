import pandas as pd
import os
import re
from rapidfuzz import fuzz

# Snakemake passes input/output paths automatically
input_metadata = snakemake.input.metadata
output_file = snakemake.output.report

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

# 1. Load Metadata
meta = pd.read_csv(input_metadata, sep="\t").set_index("Sample")

# 2. Aggregation Logic
all_dfs = []
for sample_id in meta.index:
    path = f"counts/{sample_id}_merged_counts.txt"
    if os.path.exists(path):
        df = pd.read_csv(path, sep='\s+', names=['count', 'sequence'], header=None)
        if df.empty: continue
        
        total = df['count'].sum()
        df['RPM'] = (df['count'] / total) * 1e6
        
        row = meta.loc[sample_id]
        col_name = f"R{row['Round']}_{row['Group']}_{row['Fraction']}_RPM"
        
        df = df[['sequence', 'RPM']].set_index('sequence')
        df.columns = [col_name]
        all_dfs.append(df)

# 3. Merge and Calculate Metrics
if all_dfs:
    master = pd.concat(all_dfs, axis=1).fillna(0)
    
    # --- GLOBAL SUMS ---
    r1_cols = [c for c in master.columns if 'R1' in c]
    master['R1_Global_RPM'] = master[r1_cols].sum(axis=1)
    
    r6_top_cols = [c for c in master.columns if 'R6' in c and 'Top' in c]
    master['R6_Total_Top_RPM'] = master[r6_top_cols].sum(axis=1)
    
    r6_bot_cols = [c for c in master.columns if 'R6' in c and 'Bot' in c]
    master['R6_Total_Bot_RPM'] = master[r6_bot_cols].sum(axis=1)
    
    master['R6_Global_RPM'] = master['R6_Total_Top_RPM'] + master['R6_Total_Bot_RPM']
    master['Grand_Total_RPM'] = master.filter(like='_RPM').sum(axis=1)

    # --- CALCULATED RATIOS ---
    master['Partition_Ratio_R6'] = (master['R6_Total_Top_RPM'] + 0.1) / (master['R6_Total_Bot_RPM'] + 0.1)
    master['Cumulative_Enrichment'] = (master['R6_Total_Top_RPM'] + 0.1) / (master['R1_Global_RPM'] + 0.1)
    
    # --- ROBUST AGGREGATION INDEX CALCULATION ---
    # Search for Singlet vs Double using keywords to avoid case/typo issues
    s_col = [c for c in master.columns if 'R6' in c and 'Singlet' in c and 'Top' in c]
    d_col = [c for c in master.columns if 'R6' in c and 'Double' in c and 'Top' in c]

    if s_col and d_col:
        master['Aggregation_Index'] = (master[s_col[0]] + 0.1) / (master[d_col[0]] + 0.1)
        print(f"DEBUG: Calculated Aggregation Index using {s_col[0]} / {d_col[0]}")
    else:
        print(f"DEBUG: Skipping Aggregation Index. Found Singlet: {s_col}, Found Double: {d_col}")

    master['Sequence_Length'] = master.index.str.len()

    # --- LINEAGE ASSIGNMENT ---
    master.index.name = 'Sequence'
    master = master.reset_index()
    print("Assigning lineages...")
    master['Lineage'] = master['Sequence'].map(identify_lineage)

# --- 1. GLOBAL ID GENERATION (Consistent across all files) ---
    import hashlib
    def generate_consistent_id(seq):
        # Creates a unique 8-character ID based on the sequence string
        return "seq_" + hashlib.md5(seq.encode()).hexdigest()[:8]

    print("Generating Global IDs...")
    master['Seq_ID'] = master['Sequence'].apply(generate_consistent_id)
    
    # Sort by a primary metric for the main file, but the ID stays with the sequence
    master = master.sort_values(by='Cumulative_Enrichment', ascending=False)
    # --- 2. DYNAMIC RANKING ---
    # Metrics are any numeric column we want to rank (excluding IDs and info)
    all_metrics = [c for c in master.columns if c not in ['Seq_ID', 'Sequence', 'Sequence_Length', 'Lineage'] and not c.endswith('_Rank')]
    
    for metric in all_metrics:
        master[f'{metric}_Rank'] = master[metric].rank(ascending=False, method='min').astype(int)

    # --- 3. REORDER COLUMNS ---
    rank_cols = [c for c in master.columns if '_Rank' in c]
    non_rank_cols = [c for c in master.columns if c not in (['Seq_ID', 'Sequence', 'Sequence_Length', 'Lineage'] + rank_cols)]
    final_cols = ['Seq_ID', 'Sequence', 'Lineage', 'Sequence_Length'] + non_rank_cols + rank_cols
    master = master[final_cols]

    # --- 4. EXPORT INDIVIDUAL TOP 100 LISTS ---
    top_lists_dir = os.path.join(os.path.dirname(output_file), "top_100_lists")
    os.makedirs(top_lists_dir, exist_ok=True)

    selected_winner_indices = set()
    all_winner_indices = set()

    for metric in all_metrics:
        top_df = master.sort_values(by=metric, ascending=False).head(100)
        all_winner_indices.update(top_df.index)
        
        # Track specific selection winners (ignore R1 and Grand Total for the Master Selected pool)
        if not any(x in metric for x in ['R1_', 'Grand_Total']):
            selected_winner_indices.update(top_df.index)
        
        list_path = os.path.join(top_lists_dir, f"top_100_{metric.lower()}.csv")
        top_df.to_csv(list_path, index=False)
        print(f"Saved: {list_path}")

    # --- 5. EXPORT MASTER POOLS ---
    base_path = os.path.dirname(output_file)
    master_all = master.loc[list(all_winner_indices)].sort_values(by='Cumulative_Enrichment', ascending=False)
    master_all.to_csv(os.path.join(top_lists_dir, "master_top_all_winners.csv"), index=False)

    master_selected = master.loc[list(selected_winner_indices)].sort_values(by='Cumulative_Enrichment', ascending=False)
    master_selected.to_csv(os.path.join(top_lists_dir, "master_top_selection_winners.csv"), index=False)

    # --- 6. FINAL MASTER SAVE ---
    master.to_csv(output_file, index=False)
    print(f"Report complete. Selection winners: {len(master_selected)}")