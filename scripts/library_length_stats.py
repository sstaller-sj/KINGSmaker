"""
Profile per-library sequence length distributions from count files.
Run from repo root:
    python scripts/library_length_stats.py
"""

import os
import sys
import pandas as pd


def main():
    samples = pd.read_csv('input/samples.tsv', sep='\t', dtype={'Sample': str})

    if 'Lib' not in samples.columns:
        print('ERROR: samples.tsv has no Lib column', file=sys.stderr)
        sys.exit(1)

    for lib in sorted(samples['Lib'].dropna().unique()):
        sample_ids = samples[samples['Lib'] == lib]['Sample'].tolist()

        # Aggregate length stats across all samples in this library
        length_read_counts = {}   # length -> total reads
        length_unique_counts = {} # length -> unique sequences (deduplicated across samples)
        seen_sequences = set()

        for sid in sample_ids:
            path = f'counts/{sid}_merged_counts.txt'
            if not os.path.exists(path):
                continue
            with open(path) as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < 2:
                        continue
                    count = int(parts[0])
                    seq = parts[1]
                    L = len(seq)
                    length_read_counts[L] = length_read_counts.get(L, 0) + count
                    if seq not in seen_sequences:
                        seen_sequences.add(seq)
                        length_unique_counts[L] = length_unique_counts.get(L, 0) + 1

        if not length_read_counts:
            print(f'\nLibrary {lib}: no data')
            continue

        total_reads = sum(length_read_counts.values())
        total_unique = sum(length_unique_counts.values())

        print(f"\n{'=' * 70}")
        print(f"Library {lib}  —  {len(sample_ids)} samples, {total_unique:,} unique seqs, {total_reads:,} reads")
        print(f"{'=' * 70}")

        # Read-weighted percentiles (where the typical read sits)
        sorted_lens = sorted(length_read_counts)
        cumulative = 0
        target = {0.05: None, 0.25: None, 0.50: None, 0.75: None, 0.95: None}
        for L in sorted_lens:
            cumulative += length_read_counts[L]
            frac = cumulative / total_reads
            for p in list(target):
                if target[p] is None and frac >= p:
                    target[p] = L
        print("\nRead-weighted length percentiles (typical read distribution):")
        for p in sorted(target):
            print(f"  {int(p*100):3d}%: {target[p]}")

        # Unique-sequence percentiles (where the library diversity sits)
        sorted_lens_u = sorted(length_unique_counts)
        cumulative_u = 0
        target_u = {0.05: None, 0.25: None, 0.50: None, 0.75: None, 0.95: None}
        for L in sorted_lens_u:
            cumulative_u += length_unique_counts[L]
            frac = cumulative_u / total_unique
            for p in list(target_u):
                if target_u[p] is None and frac >= p:
                    target_u[p] = L
        print("\nUnique-sequence length percentiles (where the library diversity lives):")
        for p in sorted(target_u):
            print(f"  {int(p*100):3d}%: {target_u[p]}")

        # Read-weighted histogram for lengths with >= 0.1% of total reads
        print("\nRead-weighted length histogram (lengths with >=0.1% of reads):")
        threshold = max(total_reads * 0.001, 1)
        shown = [(L, c) for L, c in sorted(length_read_counts.items()) if c >= threshold]
        if shown:
            max_count = max(c for _, c in shown)
            for L, c in shown:
                bar = '#' * max(1, int(c / max_count * 50))
                pct = c / total_reads * 100
                print(f"  {L:4d} | {c:12,d} ({pct:5.1f}%) | {bar}")

        # Suggested filter bounds: 5th-95th of unique-sequence distribution
        # (more conservative than read-weighted because dominant artifacts skew read counts)
        lo, hi = target_u[0.05], target_u[0.95]
        if lo and hi:
            print(f"\nSuggested filter bounds for Library {lib}: ({lo}, {hi})")
            print(f"  (covers 90% of unique sequences; tune ±5 nt as needed)")


if __name__ == '__main__':
    main()
