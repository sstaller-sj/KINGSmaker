#!/usr/bin/env bash
# Generate a report of the KINGSmaker pipeline run state.
# Usage from the repo root:
#     bash scripts/check_outputs.sh > pipeline_check_report.txt
# Then send pipeline_check_report.txt back for review.

set -u

div()  { echo; echo "================================================================"; echo "$1"; echo "================================================================"; }
sub()  { echo; echo "---- $1 ----"; }

# Auto-detect the R3 Top sample IDs from samples.tsv
SINGLET_SAMPLE=$(awk -F'\t' '$4==3 && $5=="Top" && $6=="Singlet" {print $1; exit}' input/samples.tsv)
DOUBLET_SAMPLE=$(awk -F'\t' '$4==3 && $5=="Top" && $6=="Doublet" {print $1; exit}' input/samples.tsv)

# ============================================================
div "1. Per-sample read counts"
# ============================================================
printf "%-40s %12s %12s %10s\n" "file" "total_reads" "unique_seqs" "U/T_ratio"
shopt -s nullglob
for f in counts/*_merged_counts.txt; do
    total=$(awk '{sum += $1} END {print sum+0}' "$f")
    unique=$(wc -l < "$f" | tr -d ' ')
    if [ "$total" -gt 0 ]; then
        ratio=$(awk -v u="$unique" -v t="$total" 'BEGIN {printf "%.3f", u/t}')
    else
        ratio="-"
    fi
    printf "%-40s %12s %12s %10s\n" "$(basename "$f")" "$total" "$unique" "$ratio"
done
shopt -u nullglob

# ============================================================
div "2. Master report structure"
# ============================================================
REPORT=results/selex_enrichment_report.csv
if [ -f "$REPORT" ]; then
    echo "Total unique sequences: $(($(wc -l < "$REPORT") - 1))"
    echo
    echo "Columns:"
    head -1 "$REPORT" | tr ',' '\n' | nl
else
    echo "MISSING: $REPORT"
fi

# ============================================================
div "3. Top-ranked sequence (Stage tag check)"
# ============================================================
if [ -f "$REPORT" ]; then
    head -2 "$REPORT" | awk -F, '
        NR==1 { n=NF; for(i=1;i<=n;i++) header[i]=$i; next }
        { for(i=1;i<=n;i++) printf "  %-30s : %s\n", header[i], $i }
    '
fi

# ============================================================
div "4. Lineage distribution"
# ============================================================
if [ -f "$REPORT" ]; then
    sub "All unique sequences"
    awk -F, 'NR>1 {print $3}' "$REPORT" | sort | uniq -c | sort -rn

    MASTER_SEL=results/top_100_lists/master_top_selection_winners.csv
    if [ -f "$MASTER_SEL" ]; then
        sub "Selection winners pool only"
        awk -F, 'NR>1 {print $3}' "$MASTER_SEL" | sort | uniq -c | sort -rn
    fi
fi

# ============================================================
div "5. R3 T20 top-20 candidates with raw counts"
# ============================================================
echo "Auto-detected R3 T20 Singlet sample: ${SINGLET_SAMPLE:-NONE}"
echo "Auto-detected R3 T20 Doublet sample: ${DOUBLET_SAMPLE:-NONE}"

emit_top20() {
    local file=$1 sample=$2
    if [ ! -f "$file" ]; then echo "MISSING: $file"; return; fi
    if [ -z "$sample" ] || [ ! -f "counts/${sample}_merged_counts.txt" ]; then
        echo "MISSING counts file for sample '$sample'"
        return
    fi
    printf "%-5s %-10s %-18s %s\n" "rank" "raw_count" "lineage" "sequence"
    local rank=0
    awk -F, 'NR>1 && NR<=21 {print $2","$3}' "$file" | while IFS=, read -r seq lineage; do
        rank=$((rank + 1))
        count=$(awk -v s="$seq" '$2==s {print $1; exit}' "counts/${sample}_merged_counts.txt")
        printf "%-5d %-10s %-18s %s\n" "$rank" "${count:-0}" "$lineage" "$seq"
    done
}

SINGLET_FILE=$(ls results/top_100_lists/top_100_r3_*singlet_top_rpm.csv 2>/dev/null | head -1)
DOUBLET_FILE=$(ls results/top_100_lists/top_100_r3_*doublet_top_rpm.csv 2>/dev/null | head -1)

sub "Singlet T20 top-20 ($SINGLET_FILE)"
emit_top20 "$SINGLET_FILE" "$SINGLET_SAMPLE"

sub "Doublet T20 top-20 ($DOUBLET_FILE)"
emit_top20 "$DOUBLET_FILE" "$DOUBLET_SAMPLE"

# ============================================================
div "6. Cross-confirmation: sequences in BOTH T20 gates (top 50 each)"
# ============================================================
S="$SINGLET_FILE"
D="$DOUBLET_FILE"
if [ -n "$S" ] && [ -n "$D" ] && [ -f "$S" ] && [ -f "$D" ]; then
    awk -F, 'NR>1 && NR<=51 {print $2}' "$S" | sort -u > /tmp/s_top50.txt
    awk -F, 'NR>1 && NR<=51 {print $2}' "$D" | sort -u > /tmp/d_top50.txt

    both_n=$(comm -12 /tmp/s_top50.txt /tmp/d_top50.txt | wc -l | tr -d ' ')
    s_n=$(comm -23 /tmp/s_top50.txt /tmp/d_top50.txt | wc -l | tr -d ' ')
    d_n=$(comm -13 /tmp/s_top50.txt /tmp/d_top50.txt | wc -l | tr -d ' ')

    echo "  In BOTH gates: $both_n"
    echo "  Singlet only:  $s_n"
    echo "  Doublet only:  $d_n"

    sub "Sequences in BOTH gates (with raw counts from each)"
    while read -r seq; do
        s_count=$(awk -v s="$seq" '$2==s {print $1; exit}' "counts/${SINGLET_SAMPLE}_merged_counts.txt")
        d_count=$(awk -v s="$seq" '$2==s {print $1; exit}' "counts/${DOUBLET_SAMPLE}_merged_counts.txt")
        printf "S:%-7s D:%-7s %s\n" "${s_count:-0}" "${d_count:-0}" "$seq"
    done < <(comm -12 /tmp/s_top50.txt /tmp/d_top50.txt)
fi

# ============================================================
div "7. Cluster structure"
# ============================================================
for cdata in results/clustering/data/*r3_*singlet_top_rpm*_clustered.csv \
             results/clustering/data/*r3_*doublet_top_rpm*_clustered.csv \
             results/clustering/data/top_100_cumulative_enrichment_clustered.csv \
             results/clustering/data/master_top_selection_winners_clustered.csv; do
    if [ -f "$cdata" ]; then
        sub "$(basename "$cdata" _clustered.csv)"
        cid_col=$(head -1 "$cdata" | tr ',' '\n' | grep -n "^Cluster_ID$" | cut -d: -f1)
        if [ -n "$cid_col" ]; then
            echo "Cluster sizes (-1 = noise/unclustered):"
            awk -F, -v c="$cid_col" 'NR>1 {print $c}' "$cdata" | sort -n | uniq -c
        fi
    fi
done

# ============================================================
div "8. Final synthesis order list"
# ============================================================
ORDER=results/final_synthesis_order_list.csv
if [ -f "$ORDER" ]; then
    echo "Total candidates: $(($(wc -l < "$ORDER") - 1))"
    echo
    sub "Columns"
    head -1 "$ORDER" | tr ',' '\n' | nl
    sub "Top 20 (key columns)"
    awk -F, '
        NR==1 {
            for(i=1;i<=NF;i++) {
                if($i=="Seq_ID") c_id=i
                if($i=="Lineage") c_ln=i
                if($i=="Cluster_ID") c_cl=i
                if($i=="King_Type") c_kt=i
                if($i=="Origin_List") c_ol=i
                if($i=="Cumulative_Enrichment") c_ce=i
                if($i=="Grand_Total_RPM") c_gt=i
            }
            printf "%-13s %-15s %-6s %-20s %-40s %-15s %-15s\n", "Seq_ID", "Lineage", "ClID", "King_Type", "Origin_List", "Cum_Enrich", "Grand_RPM"
            next
        }
        NR<=21 {
            printf "%-13s %-15s %-6s %-20s %-40s %-15s %-15s\n", $c_id, $c_ln, $c_cl, $c_kt, $c_ol, $c_ce, $c_gt
        }
    ' "$ORDER"
else
    echo "MISSING: $ORDER"
fi

# ============================================================
div "9. UMAP and overlap output files"
# ============================================================
echo "UMAP plots:"
ls -la results/clustering/plots/ 2>/dev/null | awk 'NR>1 {printf "  %s  %s\n", $5, $NF}'
echo
echo "Excel overlap workbook:"
ls -la results/metric_overlap_table.xlsx 2>/dev/null | awk '{printf "  %s bytes  %s\n", $5, $NF}'

# Cleanup
rm -f /tmp/s_top50.txt /tmp/d_top50.txt

echo
echo "================================================================"
echo "Report complete."
echo "================================================================"
