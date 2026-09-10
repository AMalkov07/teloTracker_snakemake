#!/usr/bin/env bash
# Score the curated-library pilot (7302 day0/day4/day5 in results/<s>__curatedYP on Argon):
# snapshot its outputs, day-0 null, Loss classes, chr13L truth set, circle + path reports,
# and a per-read diff against the pipeline-library results with path columns.
set -euo pipefail
cd "$(dirname "$0")/.."
SRC=/home/andrey/argon_scratch/telo_sra_runs/results
SNAP=verification/snapshot_curated; OUT=verification/reports/curated_pilot
CURLIB=verification/curated_refs/7302_features/repeatmasker_7302_all_y_primes.fasta
BED=verification/snapshot/7302_day0_with_selection/pretelomeric_labels/pretelomeric_regions_7302_day0_with_selection_simp.bed
V=_pipeline/scripts/verify_recombination.py
mkdir -p "$OUT"; rm -f "$OUT/partB_report.md"
for s in 7302_day0_with_selection 7302_day4_with_selection 7302_day5_with_selection; do
  [ -d "$SRC/${s}__curatedYP/_pipeline/recombination" ] || { echo "no pilot output for $s"; continue; }
  mkdir -p "$SNAP/$s/recombination"
  rsync -aL --include='*_features.tsv' --include='*_features.tsv.skipped' --include='*_recombination_summary.tsv' --exclude='*' \
     "$SRC/${s}__curatedYP/_pipeline/recombination/" "$SNAP/$s/recombination/"
  for f in _post_telo_trimming.tsv _post_y_prime_probe.tsv _read_summary.tsv; do ln -sf "$(realpath verification/snapshot/$s/$s$f)" "$SNAP/$s/$s$f"; done
  python $V --snapshot "$SNAP" --out "$OUT" loss-vs-truncation --sample $s --day0-bed "$BED"
done
[ -d "$SNAP/7302_day0_with_selection/recombination" ] && python $V --snapshot "$SNAP" --out "$OUT" null --sample 7302_day0_with_selection
python $V --snapshot "$SNAP" --out "$OUT" truth-set --samples 7302_day4_with_selection,7302_day5_with_selection \
   --day0-lib "$CURLIB" --fingerprint-end chr13L --truth-out verification/truth_sets/7302_circular_gains_curatedIDs.tsv
python verification/circle_report.py "$SNAP" "$OUT/circles" 7302_day4_with_selection 7302_day5_with_selection
python verification/path_report.py "$SNAP" "$OUT/path" 7302_day4_with_selection 7302_day5_with_selection 7302_day0_with_selection || true
python verification/diff_runs.py verification/snapshot_v2c_path "$SNAP" "$OUT/diff_v2c_vs_curated" 7302_day4_with_selection 7302_day5_with_selection 7302_day0_with_selection || true
echo "reports in $OUT"
