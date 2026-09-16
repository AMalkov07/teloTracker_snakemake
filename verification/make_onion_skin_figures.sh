#!/usr/bin/env bash
# Per-chromosome-end "onion skin" schematics of every read that gained Y', from the
# curated-library runs. One PNG per end per sample.
#   make_onion_skin_figures.sh <snapshot> <out_root> <sample> [sample ...]
set -euo pipefail
cd "$(dirname "$0")/.."
SNAP=$1; OUT=$2; shift 2
for s in "$@"; do
  feats=$OUT/$s/feats; mkdir -p "$feats"
  python - "$SNAP" "$s" "$feats" <<'PY'
import sys, os
sys.path.insert(0,'_pipeline/scripts')
from verify_recombination import load_features
snap, s, out = sys.argv[1:4]
df = load_features(snap, s)
if df.empty: print(f'{s}: no features'); sys.exit()
g = df[df.y_prime_recombination_status.isin(["Y' Gain","1st Y' Change","Y' Recombination"])]
n=0
for ce, sub in g.groupby('chr_end'):
    sub.to_csv(f'{out}/{s}_{ce}_features.tsv', sep='\t', index=False); n+=1
print(f'{s}: {len(g)} gain-like reads across {n} ends')
PY
  m=verification/read_id_maps/${s}_read_id_map.tsv
  if [ -f "$m" ]; then
    python _pipeline/scripts/plot_yprime_copies.py "$feats" "$s" "$OUT/$s" --schematic --id-map "$m" > /dev/null
  else
    python _pipeline/scripts/plot_yprime_copies.py "$feats" "$s" "$OUT/$s" --schematic > /dev/null
  fi
  echo "   -> $(ls $OUT/$s/*.png 2>/dev/null | wc -l) figures in $OUT/$s"
done
