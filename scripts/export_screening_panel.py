"""
Export a deduplicated flat-CSV screening panel: top-N from each non-Bot sample,
plus any explicitly-requested extra Seq_IDs (e.g. cross-T20 winners).

Usage from repo root:
    python scripts/export_screening_panel.py
    python scripts/export_screening_panel.py --n 5
    python scripts/export_screening_panel.py --n 5 --extra seq_54859f58,seq_3cdba4dd
    python scripts/export_screening_panel.py --include-bot      # include Bot fractions too
"""

import argparse
import os
import sys

import pandas as pd


def load_sample_total(sample_id):
    """Sum of raw counts in counts/{sample_id}_merged_counts.txt (unfiltered FASTQ depth)."""
    path = f"counts/{sample_id}_merged_counts.txt"
    if not os.path.exists(path):
        return 0
    total = 0
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 1:
                try:
                    total += int(parts[0])
                except ValueError:
                    pass
    return total


def lookup_raw_count(sample_id, sequence):
    """Find a sequence's raw count in counts/{sample_id}_merged_counts.txt."""
    path = f"counts/{sample_id}_merged_counts.txt"
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] == sequence:
                return int(parts[0])
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=5, help='Top N sequences per sample (default 5)')
    ap.add_argument('--extra', default='', help='Comma-separated extra Seq_IDs to include (e.g. cross-T20 winners)')
    ap.add_argument('--include-bot', action='store_true', help='Include Bot fraction samples in the panel')
    ap.add_argument('--output', default=None, help='Output CSV path (default: results/screening_panel_top<N>.csv)')
    args = ap.parse_args()

    samples = pd.read_csv('input/samples.tsv', sep='\t', dtype={'Sample': str})
    if 'Lib' not in samples.columns:
        print('ERROR: samples.tsv has no Lib column', file=sys.stderr)
        sys.exit(1)

    extra_ids = {s.strip() for s in args.extra.split(',') if s.strip()}

    panel_rows = []
    seen = set()

    # Cache unfiltered sample totals to avoid re-scanning each count file
    sample_totals = {}

    for _, row in samples.iterrows():
        sample_id = row['Sample']
        lib = row.get('Lib')
        fraction = row['Fraction']

        if pd.isna(lib):
            print(f"  skip {sample_id}: no Lib value", file=sys.stderr)
            continue
        if fraction == 'Bot' and not args.include_bot:
            continue

        fname = f"top_100_r{row['Round']}_lib{lib}_{row['Group'].lower()}_{fraction.lower()}_rpm.csv"
        path = f"results/top_100_lists/{fname}"
        if not os.path.exists(path):
            print(f"  skip {sample_id}: {path} not found", file=sys.stderr)
            continue

        df = pd.read_csv(path)
        metric_col = f"R{row['Round']}_Lib{lib}_{row['Group']}_{fraction}_RPM"
        if metric_col not in df.columns:
            print(f"  warn {sample_id}: metric col '{metric_col}' missing", file=sys.stderr)
            continue

        if sample_id not in sample_totals:
            sample_totals[sample_id] = load_sample_total(sample_id)
        unfiltered_total = sample_totals[sample_id]

        src_label = f"{sample_id}_R{row['Round']}-Lib{lib}-{row['Group']}-{fraction}"

        # Top N from this sample
        for rank_idx, srow in df.head(args.n).iterrows():
            seq_id = srow['Seq_ID']
            if seq_id in seen:
                continue
            seen.add(seq_id)
            seq = srow['Sequence']
            raw_count = lookup_raw_count(sample_id, seq)
            filtered_rpm = srow[metric_col]
            unfiltered_rpm = (raw_count / unfiltered_total * 1e6) if unfiltered_total else 0
            panel_rows.append({
                'Seq_ID': seq_id,
                'Length': len(seq),
                'Lineage': srow['Lineage'],
                'Source_Sample': src_label,
                'Raw_Count': raw_count,
                'Filtered_RPM': round(filtered_rpm, 2),
                'Unfiltered_RPM': round(unfiltered_rpm, 2),
                'Sequence': seq,
            })

        # Pick up requested extras that happen to appear in this sample's top-100 file
        if extra_ids:
            for _, srow in df.iterrows():
                seq_id = srow['Seq_ID']
                if seq_id not in extra_ids or seq_id in seen:
                    continue
                seen.add(seq_id)
                seq = srow['Sequence']
                raw_count = lookup_raw_count(sample_id, seq)
                filtered_rpm = srow[metric_col]
                unfiltered_rpm = (raw_count / unfiltered_total * 1e6) if unfiltered_total else 0
                panel_rows.append({
                    'Seq_ID': seq_id,
                    'Length': len(seq),
                    'Lineage': srow['Lineage'],
                    'Source_Sample': f"(extra) {src_label}",
                    'Raw_Count': raw_count,
                    'Filtered_RPM': round(filtered_rpm, 2),
                    'Unfiltered_RPM': round(unfiltered_rpm, 2),
                    'Sequence': seq,
                })

    if not panel_rows:
        print('No rows produced. Did the pipeline finish?', file=sys.stderr)
        sys.exit(1)

    missing_extras = extra_ids - seen
    if missing_extras:
        print(f"  warn: requested extras not found in any sample's top-100: {sorted(missing_extras)}", file=sys.stderr)

    out_path = args.output or f"results/screening_panel_top{args.n}.csv"
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    pd.DataFrame(panel_rows).to_csv(out_path, index=False)
    print(f"Wrote {out_path}: {len(panel_rows)} unique sequences")


if __name__ == '__main__':
    main()
