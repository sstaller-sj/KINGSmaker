import pandas as pd
import os
from openpyxl.styles import PatternFill
from openpyxl.formatting.rule import CellIsRule

top_100_dir = "results/top_100_lists/"
output_path = "results/metric_overlap_table.xlsx"

targets = {
    "Selection_Winners": os.path.join(top_100_dir, "master_top_selection_winners.csv"),
    "All_Winners": os.path.join(top_100_dir, "master_top_all_winners.csv")
}

files = [f for f in os.listdir(top_100_dir) if f.startswith("top_100_") and f.endswith(".csv")]

with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
    for sheet_name, file_path in targets.items():
        if not os.path.exists(file_path):
            continue

        master = pd.read_csv(file_path)
        overlap_df = master[['Seq_ID', 'Sequence', 'Lineage']].copy()

        for f in files:
            metric_label = f.replace("top_100_", "").replace(".csv", "").upper()
            temp_df = pd.read_csv(os.path.join(top_100_dir, f))
            top_ids = set(temp_df['Seq_ID'])
            overlap_df[metric_label] = overlap_df['Seq_ID'].apply(lambda x: 1 if x in top_ids else 0)

        metric_cols = [c for c in overlap_df.columns if c not in ['Seq_ID', 'Sequence', 'Lineage']]
        overlap_df['Total_Hits'] = overlap_df[metric_cols].sum(axis=1)
        overlap_df = overlap_df.sort_values(by='Total_Hits', ascending=False)

        overlap_df.to_excel(writer, index=False, sheet_name=sheet_name)

        worksheet = writer.sheets[sheet_name]
        green_fill = PatternFill(start_color='C6EFCE', end_color='C6EFCE', fill_type='solid')

        last_col_idx = len(overlap_df.columns)
        last_col_letter = chr(64 + last_col_idx) if last_col_idx <= 26 else "ZZ"
        cell_range = f"D2:{last_col_letter}{len(overlap_df) + 1}"

        worksheet.conditional_formatting.add(
            cell_range,
            CellIsRule(operator='equal', formula=['1'], fill=green_fill)
        )

print(f"Multi-sheet overlap table saved to {output_path}")
