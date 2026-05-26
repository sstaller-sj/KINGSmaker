"""
Export per-sample top-N candidates to a multi-sheet Excel workbook.

Usage from repo root:
    python scripts/export_per_sample_top.py          # default N=30
    python scripts/export_per_sample_top.py 50       # N=50
"""

import os
import sys
from collections import defaultdict

import pandas as pd


def load_raw_counts(sample_id):
    """Map sequence -> raw count from counts/{sample_id}_merged_counts.txt"""
    path = f"counts/{sample_id}_merged_counts.txt"
    if not os.path.exists(path):
        return {}
    counts = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            counts[parts[1]] = int(parts[0])
    return counts


def build_sheet_name(sample_id, row):
    lib = row.get('Lib', '?')
    rnd = row.get('Round', '?')
    frac = row.get('Fraction', '?')
    grp = row.get('Group', '?')
    name = f"{sample_id}_Lib{lib}_R{rnd}{frac}_{grp}"
    # Excel sheet names are capped at 31 chars and can't contain : \ / ? * [ ]
    for bad in r':\/?*[]':
        name = name.replace(bad, '_')
    return name[:31]


def main():
    n_per_sample = int(sys.argv[1]) if len(sys.argv) > 1 else 30

    samples = pd.read_csv('input/samples.tsv', sep='\t', dtype={'Sample': str})

    per_sample_sheets = {}                  # sheet_name -> DataFrame (rank-ordered)
    seq_appearances = defaultdict(list)     # sequence -> list of {sample_id, rank, seq_id, lineage, length, fraction}

    for _, row in samples.iterrows():
        sample_id = row['Sample']
        lib = row.get('Lib')
        if pd.isna(lib):
            print(f"  skip {sample_id}: no Lib value in samples.tsv", file=sys.stderr)
            continue

        # Reconstruct the per-sample top-100 filename using the same convention as analyze_selex_enrichment.py
        fname = f"top_100_r{row['Round']}_lib{lib}_{row['Group'].lower()}_{row['Fraction'].lower()}_rpm.csv"
        path = f"results/top_100_lists/{fname}"

        if not os.path.exists(path):
            print(f"  skip {sample_id}: {path} not found", file=sys.stderr)
            continue

        df = pd.read_csv(path)
        if df.empty:
            print(f"  skip {sample_id}: top-100 file is empty", file=sys.stderr)
            continue

        metric_col = f"R{row['Round']}_Lib{lib}_{row['Group']}_{row['Fraction']}_RPM"
        if metric_col not in df.columns:
            print(f"  warn {sample_id}: metric column '{metric_col}' not found, skipping", file=sys.stderr)
            continue

        top_n = df.head(n_per_sample).copy().reset_index(drop=True)

        raw_counts = load_raw_counts(sample_id)
        top_n['Rank'] = top_n.index + 1
        top_n['Raw_Count'] = top_n['Sequence'].map(lambda s: raw_counts.get(s, 0))
        top_n['Length'] = top_n['Sequence'].str.len()
        top_n['RPM'] = top_n[metric_col]

        sheet_df = top_n[['Rank', 'Seq_ID', 'Raw_Count', 'RPM', 'Lineage', 'Length', 'Sequence']]
        sheet_name = build_sheet_name(sample_id, row)
        per_sample_sheets[sheet_name] = sheet_df

        for _, srow in top_n.iterrows():
            seq_appearances[srow['Sequence']].append({
                'sample_id': sample_id,
                'rank': int(srow['Rank']),
                'seq_id': srow['Seq_ID'],
                'lineage': srow['Lineage'],
                'length': int(srow['Length']),
                'fraction': row['Fraction'],
            })

    if not per_sample_sheets:
        print("No data to export. Did the pipeline finish?", file=sys.stderr)
        sys.exit(1)

    # Combined sheet — dedupe by sequence, surface cross-sample hits
    # Drop only Bot-only sequences (no Top signal at all). Tolerate Bot presence if Top signal also exists.
    combined_rows = []
    dropped_bot_only = 0
    for seq, appearances in seq_appearances.items():
        sample_ids = sorted({a['sample_id'] for a in appearances})
        best_rank = min(a['rank'] for a in appearances)
        top_hits = sum(1 for a in appearances if a['fraction'] == 'Top')
        bot_hits = sum(1 for a in appearances if a['fraction'] == 'Bot')
        all_hits = sum(1 for a in appearances if a['fraction'] == 'All')

        # Drop only if the sequence has Bot presence and zero Top presence — FACS-rejected with no positive signal
        if bot_hits > 0 and top_hits == 0:
            dropped_bot_only += 1
            continue

        # Selectivity score: positive when Top dominates, negative when Bot does, 0 when neutral or no sort happened
        selectivity = top_hits - bot_hits

        combined_rows.append({
            'Seq_ID': appearances[0]['seq_id'],
            'Sequence': seq,
            'Lineage': appearances[0]['lineage'],
            'Length': appearances[0]['length'],
            'Total_Hits': len(sample_ids),
            'Top_Hits': top_hits,
            'Bot_Hits': bot_hits,
            'All_Hits': all_hits,
            'Selectivity': selectivity,
            'Best_Rank': best_rank,
            'Samples_Present': ', '.join(sample_ids),
        })

    # Sort: highest selectivity (Top - Bot) first, then most Top hits, then most total hits, then best rank.
    combined_df = pd.DataFrame(combined_rows).sort_values(
        by=['Selectivity', 'Top_Hits', 'Total_Hits', 'Best_Rank'],
        ascending=[False, False, False, True]
    ).reset_index(drop=True)
    print(f"Combined sheet: dropped {dropped_bot_only} Bot-only sequences (no Top presence)", file=sys.stderr)

    out_path = f"results/per_sample_top{n_per_sample}.xlsx"
    os.makedirs('results', exist_ok=True)
    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        combined_df.to_excel(writer, sheet_name='Combined', index=False)
        for sheet_name, df in per_sample_sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    total_rows = sum(len(d) for d in per_sample_sheets.values())
    print(f"Wrote {out_path}: 1 Combined sheet ({len(combined_df)} unique seqs) + {len(per_sample_sheets)} per-sample sheets ({total_rows} total rows)")


if __name__ == '__main__':
    main()
