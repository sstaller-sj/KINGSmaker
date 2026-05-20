import pandas as pd
import os

master = pd.read_csv(snakemake.input.report)

all_king_data = []
for path in snakemake.input.kings:
    origin_name = os.path.basename(path).replace("_clustered_kings.csv", "")
    df_king = pd.read_csv(path)
    df_king['Origin_List'] = origin_name
    all_king_data.append(df_king)

combined_kings = pd.concat(all_king_data).drop_duplicates(subset=['Seq_ID'])

final_order = pd.merge(
    combined_kings[['Seq_ID', 'Cluster_ID', 'King_Type', 'Origin_List']],
    master,
    on='Seq_ID',
    how='left'
)

# Cumulative_Enrichment becomes the "x" placeholder string when First/Last Stage tags aren't set — fall back to total abundance for the sort
sort_key = 'Cumulative_Enrichment' if pd.api.types.is_numeric_dtype(final_order['Cumulative_Enrichment']) else 'Grand_Total_RPM'
final_order = final_order.sort_values(by=sort_key, ascending=False)

final_order.to_csv(snakemake.output.order_list, index=False)
print(f"Synthesis order list created with {len(final_order)} unique sequences.")
