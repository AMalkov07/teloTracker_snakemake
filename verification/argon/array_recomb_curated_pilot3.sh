#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd
#$ -N recomb_cur3
#$ -t 1-3
# Pilot: recombination step (11+) with the CURATED strain Y' library instead of
# the pipeline's silhouette-clustered one, for 7302 day4 + day5. Runs in NEW
# run dirs / results dirs (<sample>__curatedYP); the existing v1/v2 outputs in
# results/<sample>/ are not touched (steps 0-6 + alignment outputs are reused
# through a symlink farm).
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

TAG=curatedYP
RUN_ROOT=/nfsscratch/amalkova/telo_sra_runs
RESULTS=$RUN_ROOT/results
STAGE=$RUN_ROOT/v2_stage
SRA=/nfsscratch/amalkova/sra_jobs/outputs
MANIFEST=$STAGE/recomb_curated_pilot2_manifest.txt
CODE=$STAGE/path

read -r SAMPLE DAY0 < <(sed -n "${SGE_TASK_ID}p" "$MANIFEST")
[ -n "${SAMPLE:-}" ] || { echo "ERROR: no sample for task ${SGE_TASK_ID}"; exit 1; }
STRAIN=${DAY0%%_*}
CURLIB=$STAGE/yprime_lib_curated_${STRAIN}.fasta
[ -f "$CURLIB" ] || { echo "ERROR: missing $CURLIB"; exit 1; }
SRC_RUN=$RUN_ROOT/$SAMPLE
RUNDIR=$RUN_ROOT/${SAMPLE}__$TAG
NEWRES=$RESULTS/${SAMPLE}__$TAG
echo "[task ${SGE_TASK_ID}] SAMPLE=$SAMPLE DAY0=$DAY0 lib=$CURLIB host=$(hostname) start=$(date)"

# run dir = copy of the (v2-patched) sample run dir, minus data
mkdir -p "$RUNDIR"
rsync -a --exclude "results/" --exclude "samples_dorado_basecalled/" --exclude ".snakemake/" "$SRC_RUN/" "$RUNDIR/"
for f in analyze_features.py yprime_path.py aggregate_recombination.py recombination_utils.py; do cp "$CODE/$f" "$RUNDIR/_pipeline/scripts/$f"; done
python3 "$CODE/patch_snakefile_v2.py" "$RUNDIR/_pipeline/Snakefile"     # adds --y-prime-id-level when missing
mkdir -p "$RUNDIR/results" "$RUNDIR/samples_dorado_basecalled"
ln -sfn "$RESULTS/$DAY0" "$RUNDIR/results/$DAY0"
# input reads: the SRA download when present, otherwise whatever the source run dir used
# (the day-0 runs were fed an unzipped .fastq under their own name)
rm -f "$RUNDIR/samples_dorado_basecalled/$SAMPLE".fastq.gz "$RUNDIR/samples_dorado_basecalled/$SAMPLE".fastq "$RUNDIR/samples_dorado_basecalled/$SAMPLE".bam
if [ -f "$SRA/$SAMPLE.fastq.gz" ]; then
  ln -sfn "$SRA/$SAMPLE.fastq.gz" "$RUNDIR/samples_dorado_basecalled/$SAMPLE.fastq.gz"
else
  src_in=""
  for ext in bam fastq fastq.gz; do                      # (no ls|head: pipefail + set -e would abort silently)
    [ -f "$SRC_RUN/samples_dorado_basecalled/$SAMPLE.$ext" ] && { src_in="$SRC_RUN/samples_dorado_basecalled/$SAMPLE.$ext"; break; }
  done
  [ -n "$src_in" ] || { echo "ERROR: no input reads for $SAMPLE in $SRA or $SRC_RUN"; exit 1; }
  ln -sfn "$(readlink -f "$src_in")" "$RUNDIR/samples_dorado_basecalled/$(basename "$src_in")"
fi

# results dir = symlink farm onto the existing outputs, with step-11+ outputs removed
PREV_TAG=${PREV_TAG:-v2b}
if [ ! -d "$NEWRES/_pipeline" ]; then
  mkdir -p "$NEWRES/_pipeline"
  cp -as "$RESULTS/$SAMPLE/_pipeline/." "$NEWRES/_pipeline/"
  rm -rf "$NEWRES/_pipeline/recombination_v1" "$NEWRES/_pipeline/recombination_v2a"
elif ls "$NEWRES"/_pipeline/recombination/${SAMPLE}_chr*_features.tsv >/dev/null 2>&1 && [ ! -d "$NEWRES/_pipeline/recombination_$PREV_TAG" ]; then
  cp -a "$NEWRES/_pipeline/recombination" "$NEWRES/_pipeline/recombination_$PREV_TAG"; echo "kept previous pilot run as recombination_$PREV_TAG"
fi
rm -rf "$NEWRES/_pipeline/recombination_events" "$NEWRES/_pipeline/graphs/recombination_tracks"
rm -f "$NEWRES"/_pipeline/recombination/${SAMPLE}_chr*_features.tsv "$NEWRES"/_pipeline/recombination/${SAMPLE}_chr*_features.tsv.skipped \
      "$NEWRES"/_pipeline/recombination/${SAMPLE}_recombination_summary.tsv
ln -sfn "$NEWRES" "$RUNDIR/results/$SAMPLE"

# config: point the recombination step at the curated library
cd "$RUNDIR"
SPLIB=$STAGE/spacer_variable/$DAY0/pairings_for_spacers/${DAY0}_pairings; [ -d "$SPLIB" ] || SPLIB=""
CURLIB="$CURLIB" SPLIB="$SPLIB" python3 - <<'PY'
import os, re
p = "_pipeline/config.yaml"; s = open(p).read(); lib = os.environ["CURLIB"]
s = re.sub(r'^\s*y_prime_lib_override:.*\n', '', s, flags=re.M)
s = s.replace("references:", f'references:\n  y_prime_lib_override: "{lib}"', 1)
level = os.environ.get("ID_LEVEL", "family")            # qsub -v ID_LEVEL=variant keeps Red-Light/Red-Dark apart
s = re.sub(r'^y_prime_id_level:.*\n', '', s, flags=re.M)
s = f'y_prime_id_level: "{level}"\n' + s
splib = os.environ.get("SPLIB", "")                      # spacer-only library, when staged
if splib:
    s = re.sub(r'^(\s*spacer_lib_dir:).*$', r'\1 "' + splib + '"', s, flags=re.M)
open(p, "w").write(s); print("y_prime_lib_override ->", lib, "| y_prime_id_level:", level, "| spacer_lib_dir:", splib or "(unchanged)")
PY
grep -n "y_prime_lib" _pipeline/config.yaml

snakemake -s _pipeline/Snakefile all -c 56 --rerun-triggers mtime -n 2>&1 | grep -E "^(recombination|extract|all|total)" || true
snakemake -s _pipeline/Snakefile all -c 56 --rerun-triggers mtime
echo "[task ${SGE_TASK_ID}] DONE $SAMPLE end=$(date)"
echo "Outputs: $NEWRES/_pipeline/recombination/${SAMPLE}_recombination_summary.tsv"
