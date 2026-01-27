import pandas as pd
import os

# 1. Load the master enrichment data
master = pd.read_csv(snakemake.input.report)

all_king_data = []

# 2. Process each king file
for path in snakemake.input.kings:
    # Extract the name of the list from the filename
    origin_name = os.path.basename(path).replace("_clustered_kings.csv", "")
    
    df_king = pd.read_csv(path)
    # Add a column to track where this sequence came from
    df_king['Origin_List'] = origin_name
    all_king_data.append(df_king)

# 3. Combine all kings and remove duplicates
# If a sequence is a king in multiple lists, we keep the first one found
combined_kings = pd.concat(all_king_data).drop_duplicates(subset=['Seq_ID'])

# 4. Merge with Master data to get full metrics (Ratios, Lineage, etc.)
# We merge on Seq_ID to pull in all the columns from your selex_analysis report
final_order = pd.merge(
    combined_kings[['Seq_ID', 'Cluster_ID', 'King_Type', 'Origin_List']], 
    master, 
    on='Seq_ID', 
    how='left'
)

# 5. Final Sort: Usually helpful to sort by enrichment or total hits
final_order = final_order.sort_values(by='Cumulative_Enrichment', ascending=False)

# 6. Save
final_order.to_csv(snakemake.output.order_list, index=False)
print(f"Synthesis order list created with {len(final_order)} unique sequences.")