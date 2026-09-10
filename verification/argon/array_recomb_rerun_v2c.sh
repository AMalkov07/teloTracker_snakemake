#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd
#$ -N recomb_v2c
#$ -t 1-25
# Re-run ONLY the recombination analysis (step 11+) with the v2 attribution for
# the 23 timepoint/survivor samples + the 7302/7172 day-0 self-runs, in their
# existing run dirs under /nfsscratch/amalkova/telo_sra_runs/<sample>.
#
# Per sample:
#   1. copy the v2 scripts (staged in $STAGE) into the run dir
#   2. patch the run-dir Snakefile (telo/probe inputs, --attribution-mode)
#   3. keep the v1 outputs:  recombination/ -> recombination_v1/  (once)
#   4. delete only the step-11+ outputs that must be regenerated
#   5. snakemake ... --rerun-triggers mtime   (essential: otherwise the changed
#      script/params would re-trigger steps 0-6 as well)
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

RUN_ROOT=/nfsscratch/amalkova/telo_sra_runs
RESULTS=$RUN_ROOT/results
STAGE=$RUN_ROOT/v2_stage
CODE=$STAGE/path                         # analyze_features + yprime_path + aggregate + utils
MANIFEST=$STAGE/recomb_rerun_v2_manifest.txt

read -r SAMPLE DAY0 < <(sed -n "${SGE_TASK_ID}p" "$MANIFEST")
[ -n "${SAMPLE:-}" ] || { echo "ERROR: no sample for task ${SGE_TASK_ID}"; exit 1; }
RUNDIR="$RUN_ROOT/$SAMPLE"
echo "[task ${SGE_TASK_ID}] SAMPLE=${SAMPLE} DAY0=${DAY0} host=$(hostname) start=$(date)"
[ -f "$RUNDIR/_pipeline/Snakefile" ] || { echo "ERROR: no run dir $RUNDIR"; exit 1; }
cd "$RUNDIR"

# 1. v2 scripts
for f in analyze_features.py yprime_path.py aggregate_recombination.py recombination_utils.py; do
  cp "$CODE/$f" "$RUNDIR/_pipeline/scripts/$f"
done
# spacer-only library (variable-size) for this strain's day-0 reference
SPLIB=$STAGE/spacer_variable/$DAY0/pairings_for_spacers/${DAY0}_pairings
[ -d "$SPLIB" ] || { echo "ERROR: missing spacer library $SPLIB"; exit 1; }
SPLIB="$SPLIB" python3 - <<'PY'
import os, re
p = "_pipeline/config.yaml"; s = open(p).read()
s = re.sub(r'^(\s*spacer_lib_dir:).*$', r'\1 "' + os.environ["SPLIB"] + '"', s, flags=re.M)
open(p, "w").write(s); print("spacer_lib_dir ->", os.environ["SPLIB"])
PY
# 2. Snakefile
python3 "$STAGE/patch_snakefile_v2.py" "$RUNDIR/_pipeline/Snakefile"

# 3./4. keep v1, remove only what step 11+ regenerates
P="$RESULTS/$SAMPLE/_pipeline"
# keep the previous run: first one becomes recombination_v1, later ones recombination_<PREV_TAG>
PREV_TAG=${PREV_TAG:-v2b}
if [ ! -d "$P/recombination_v1" ]; then
  cp -a "$P/recombination" "$P/recombination_v1"; echo "backed up recombination -> recombination_v1"
elif [ ! -d "$P/recombination_$PREV_TAG" ] && ls "$P"/recombination/${SAMPLE}_chr*_features.tsv >/dev/null 2>&1; then
  cp -a "$P/recombination" "$P/recombination_$PREV_TAG"; echo "backed up recombination -> recombination_$PREV_TAG"
fi
rm -f "$P"/recombination/${SAMPLE}_chr*_features.tsv "$P"/recombination/${SAMPLE}_chr*_features.tsv.skipped \
      "$P"/recombination/${SAMPLE}_recombination_summary.tsv
rm -rf "$P/recombination_events" "$P/graphs/recombination_tracks"

# 5. run
cd "$RUNDIR"
snakemake -s _pipeline/Snakefile all -c 56 --rerun-triggers mtime -n 2>&1 | grep -E "^(recombination_analyze|recombination_summary|extract_recombination_events|recombination_track_plots|all|total)" || true
snakemake -s _pipeline/Snakefile all -c 56 --rerun-triggers mtime

echo "[task ${SGE_TASK_ID}] DONE ${SAMPLE} end=$(date)"
echo "Outputs: $P/recombination/${SAMPLE}_recombination_summary.tsv  (v1 kept in $P/recombination_v1)"
