#!/usr/bin/env bash
# Part A: verify every pipeline-built day-0 reference against the curated
# reference of its strain. Inputs come from verification/snapshot/<ref>/
# (rsynced from /nfsscratch/amalkova/telo_sra_runs/results/<ref>/_pipeline/).
set -euo pipefail
cd "$(dirname "$0")/.."

SNAP=verification/snapshot
OUT=verification/reports/partA
CUR=${CUR:-verification/curated_refs}
SUMMARY=$OUT/partA_summary.tsv
mkdir -p "$OUT"; rm -f "$SUMMARY"

REFS=(6991_day0 6991_day0_TeloTag 6991_day0_TeloTag_with_selection 6991_day0_reference
      6991_day0_reference_promethion 6991_day0_with_selection 6991_day0_with_selection_repeat
      6991_day0_with_selection_repeat2 7172_day0_with_selection 7302_day0_with_selection)

for ref in "${REFS[@]}"; do
  strain=${ref%%_*}
  lab=$SNAP/$ref/pretelomeric_labels
  working=$(ls "$lab"/yprime_clustering_*/*_working_yprimes.fasta 2>/dev/null | head -1 || true)
  python _pipeline/scripts/verify_day0_reference.py \
    --ours-bed    "$lab/pretelomeric_regions_${ref}_simp.bed" \
    --ours-lib    "$lab/extracted_yprimes_${ref}.fasta" \
    --probe-blast "$lab/pretelomeric_regions_${ref}_probe_blast.txt" \
    --working-lib "${working}" \
    --run-config  "$SNAP/$ref/run_config.yaml" \
    --curated-bed "$CUR/${strain}_features/${strain}_final_features.bed" \
    --curated-lib "$CUR/${strain}_features/repeatmasker_${strain}_all_y_primes.fasta" \
    --strain "$strain" --ref-name "$ref" \
    --out-prefix "$OUT/$ref" --summary-tsv "$SUMMARY"
done
echo; column -t -s $'\t' "$SUMMARY"
