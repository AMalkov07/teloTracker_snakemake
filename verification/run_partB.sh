#!/usr/bin/env bash
# Part B: recombination verification (B1-B5) on a snapshot of recombination
# outputs. Usage: run_partB.sh <snapshot_dir> <report_dir>
#   v1 (as run on Argon):  run_partB.sh verification/snapshot verification/reports/v1
#   v2 (reprocessed):      run_partB.sh verification/snapshot_v2 verification/reports/v2
set -euo pipefail
cd "$(dirname "$0")/.."
SNAP=${1:-verification/snapshot}
OUT=${2:-verification/reports/v1}
LAB=verification/snapshot            # day-0 beds/libs always come from the original snapshot
V=_pipeline/scripts/verify_recombination.py
mkdir -p "$OUT"; rm -f "$OUT/partB_report.md"
echo "# Part B recombination verification — snapshot: $SNAP" > "$OUT/partB_report.md"

D7302=7302_day0_with_selection; D7172=7172_day0_with_selection
BED7302=$LAB/$D7302/pretelomeric_labels/pretelomeric_regions_${D7302}_simp.bed
BED7172=$LAB/$D7172/pretelomeric_labels/pretelomeric_regions_${D7172}_simp.bed
LIB7302=$LAB/$D7302/pretelomeric_labels/extracted_yprimes_${D7302}.fasta

# B1 null: day-0 self-runs
for d in $D7302 $D7172; do python $V --snapshot "$SNAP" --out "$OUT" null --sample $d; done

# B2 replicates
python $V --snapshot "$SNAP" --out "$OUT" replicates --sample-a 7172_day4_with_selection --sample-b 7172_day4_with_selection_repeat
for pair in "6991_day0_with_selection 6991_day0_with_selection_repeat" "6991_day0_with_selection_repeat 6991_day0_with_selection_repeat2"; do
  set -- $pair; [ -d "$SNAP/$1/recombination" ] && [ -d "$SNAP/$2/recombination" ] && \
    python $V --snapshot "$SNAP" --out "$OUT" replicates --sample-a $1 --sample-b $2 || true
done

# B3 positive control: 7172 chr11L -> chr11R
python $V --snapshot "$SNAP" --out "$OUT" positive-control \
  --samples 7172_day4_with_selection,7172_day6_with_selection,7172_day9_with_selection --chr-end chr11L --expected-source chr11R

# B4 Loss classification
for s in $D7302 7302_day5_with_selection 7302_day6_with_selection_repeat; do
  python $V --snapshot "$SNAP" --out "$OUT" loss-vs-truncation --sample $s --day0-bed "$BED7302"; done
python $V --snapshot "$SNAP" --out "$OUT" loss-vs-truncation --sample $D7172 --day0-bed "$BED7172"

# B5 truth set: chr13L alternating fingerprint in 7302 day5/day6
python $V --snapshot "$SNAP" --out "$OUT" truth-set \
  --samples 7302_day5_with_selection,7302_day6_with_selection_repeat,7302_day4_with_selection,7302_day9_with_selection_repeat \
  --day0-lib "$LIB7302" --fingerprint-end chr13L --truth-out verification/truth_sets/7302_circular_gains.tsv
echo; echo "report: $OUT/partB_report.md"
