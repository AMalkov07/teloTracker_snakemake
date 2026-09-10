#!/usr/bin/env bash
# Reprocess every snapshot sample's *_features.tsv with the v2 attribution
# (no BLAST / RepeatMasker; only the Y' comparison + reconciliation + path are
# re-derived), then re-aggregate summaries.
# Usage: reprocess_snapshot.sh [mode=v2] [snapshot=verification/snapshot] [out=verification/snapshot_v2] [labels=verification/snapshot]
#   labels = snapshot holding the day-0 pretelomeric_labels (bed + Y' library)
set -euo pipefail
cd "$(dirname "$0")/.."
MODE=${1:-v2}
SNAP=${2:-verification/snapshot}
OUT=${3:-verification/snapshot_v2}
LAB=${4:-verification/snapshot}
AF=_pipeline/scripts/analyze_features.py
AG=_pipeline/scripts/aggregate_recombination.py

day0_of() {            # sample -> day-0 reference name
  case "$1" in
    7302_*) echo 7302_day0_with_selection ;;
    7172_*) echo 7172_day0_with_selection ;;
    6991_*) echo "$1" ;;                      # 6991 samples are day-0 self-runs
    *) echo "$1" ;;
  esac
}

for sdir in "$SNAP"/*/recombination; do
  s=$(basename "$(dirname "$sdir")")
  d0=$(day0_of "$s")
  bed=$LAB/$d0/pretelomeric_labels/pretelomeric_regions_${d0}_simp.bed
  lib=$LAB/$d0/pretelomeric_labels/extracted_yprimes_${d0}.fasta
  [ -f "$bed" ] && [ -f "$lib" ] || { echo "skip $s (no day-0 bed/lib for $d0)"; continue; }
  mkdir -p "$OUT/$s/recombination"
  for f in _post_telo_trimming.tsv _post_y_prime_probe.tsv _read_summary.tsv; do
    [ -e "$SNAP/$s/$s$f" ] && ln -sf "$(realpath "$SNAP/$s/$s$f")" "$OUT/$s/$s$f"
  done
  n=0
  for ftsv in "$sdir"/${s}_chr*_features.tsv; do
    ce=$(basename "$ftsv" | sed -E "s/^${s}_(chr[0-9]+[LR])_features.tsv$/\1/")
    out="$OUT/$s/recombination/$(basename "$ftsv")"
    if [ -f "$ftsv.skipped" ]; then cp "$ftsv" "$out"; cp "$ftsv.skipped" "$out.skipped"; continue; fi
    python $AF --reprocess-tsv "$ftsv" --day0-bed "$bed" --y-prime-lib "$lib" --chr-end "$ce" \
      --telo-tsv "$SNAP/$s/${s}_post_telo_trimming.tsv" --probe-tsv "$SNAP/$s/${s}_post_y_prime_probe.tsv" \
      --attribution-mode "$MODE" --output-tsv "$out" > /dev/null
    n=$((n+1))
  done
  python $AG --recombination-dir "$OUT/$s/recombination" --base-name "$s" \
    --output-summary "$OUT/$s/recombination/${s}_recombination_summary.tsv" > /dev/null
  echo "$s: $n chr_ends reprocessed ($MODE)"
done
