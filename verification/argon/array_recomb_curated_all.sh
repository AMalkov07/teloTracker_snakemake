#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd
#$ -N recomb_cur_all
#$ -t 1-33
# Re-run the recombination analysis (step 11+) with the CURATED strain Y' library for every
# sample whose strain has one (8x 6991 day-0 self-runs, 7172 + 7302 day-0, and the 7302/7172
# time points and survivors).  Nothing existing is overwritten: each sample gets its own run
# dir <sample>__curatedYP and results dir results/<sample>__curatedYP, whose steps 0-6 are
# symlinks into the sample's existing outputs (cp -as), so only step 11+ is recomputed.
#
#   Y' library     : curated repeatmasker_<strain>_all_y_primes.fasta   (y_prime_lib_override)
#   Y' ID level    : variant  (ID2_Red-Light vs ID2_Red-Dark kept apart; collapses to family later)
#   spacer library : the spacer-only (variable-size) library for this day-0 reference
#   attribution    : v2, and every read records the matched library ENTRY (y_prime_entries)
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

TAG=curatedYP
ID_LEVEL=${ID_LEVEL:-variant}
RUN_ROOT=/nfsscratch/amalkova/telo_sra_runs
RESULTS=$RUN_ROOT/results
STAGE=$RUN_ROOT/v2_stage
CODE=$STAGE/path
SRA=/nfsscratch/amalkova/sra_jobs/outputs
MANIFEST=$STAGE/recomb_curated_all_manifest.txt

read -r SAMPLE DAY0 < <(sed -n "${SGE_TASK_ID}p" "$MANIFEST")
[ -n "${SAMPLE:-}" ] || { echo "ERROR: no sample for task ${SGE_TASK_ID}"; exit 1; }
STRAIN=${DAY0%%_*}
CURLIB=$STAGE/yprime_lib_curated_${STRAIN}.fasta
SPLIB=$STAGE/spacer_variable/$DAY0/pairings_for_spacers/${DAY0}_pairings
SRC_RUN=$RUN_ROOT/$SAMPLE
RUNDIR=$RUN_ROOT/${SAMPLE}__$TAG
NEWRES=$RESULTS/${SAMPLE}__$TAG
echo "[task ${SGE_TASK_ID}] SAMPLE=$SAMPLE DAY0=$DAY0 STRAIN=$STRAIN id_level=$ID_LEVEL host=$(hostname) start=$(date)"
[ -f "$CURLIB" ] || { echo "ERROR: missing curated library $CURLIB"; exit 1; }
[ -d "$SPLIB" ]  || { echo "ERROR: missing spacer library $SPLIB"; exit 1; }
[ -f "$SRC_RUN/_pipeline/Snakefile" ] || { echo "ERROR: no run dir $SRC_RUN"; exit 1; }

# --- run dir: a copy of the sample's run dir with the current code -----------------
mkdir -p "$RUNDIR"
rsync -a --exclude "results/" --exclude "samples_dorado_basecalled/" --exclude ".snakemake/" "$SRC_RUN/" "$RUNDIR/"
for f in analyze_features.py yprime_path.py aggregate_recombination.py recombination_utils.py; do
  cp "$CODE/$f" "$RUNDIR/_pipeline/scripts/$f"
done
python3 "$CODE/patch_snakefile_v2.py" "$RUNDIR/_pipeline/Snakefile"
mkdir -p "$RUNDIR/results" "$RUNDIR/samples_dorado_basecalled"

# --- results dir: symlink farm over the sample's existing outputs ------------------
# Built atomically via a temp dir and marked complete, so a task that dies part way
# through leaves nothing half-built for the next run to trip over.
PREV_TAG=${PREV_TAG:-prev}
if [ -f "$NEWRES/.farm_ready" ]; then
  if ls "$NEWRES"/_pipeline/recombination/${SAMPLE}_chr*_features.tsv >/dev/null 2>&1 \
     && [ ! -d "$NEWRES/_pipeline/recombination_$PREV_TAG" ]; then
    cp -a "$NEWRES/_pipeline/recombination" "$NEWRES/_pipeline/recombination_$PREV_TAG"
    echo "kept the previous curated run as recombination_$PREV_TAG"
  fi
else
  rm -rf "$NEWRES" "$NEWRES.tmp"
  mkdir -p "$NEWRES.tmp/_pipeline"
  cp -as "$RESULTS/$SAMPLE/_pipeline/." "$NEWRES.tmp/_pipeline/"     # symlinks to steps 0-6; originals never written
  rm -rf "$NEWRES.tmp"/_pipeline/recombination_v1 "$NEWRES.tmp"/_pipeline/recombination_v2a \
         "$NEWRES.tmp"/_pipeline/recombination_v2b "$NEWRES.tmp"/_pipeline/recombination_$PREV_TAG
  mkdir -p "$NEWRES"; mv "$NEWRES.tmp/_pipeline" "$NEWRES/_pipeline"; rmdir "$NEWRES.tmp"
  touch "$NEWRES/.farm_ready"
fi
rm -rf "$NEWRES/_pipeline/recombination_events" "$NEWRES/_pipeline/graphs/recombination_tracks"
rm -f "$NEWRES"/_pipeline/recombination/${SAMPLE}_chr*_features.tsv \
      "$NEWRES"/_pipeline/recombination/${SAMPLE}_chr*_features.tsv.skipped \
      "$NEWRES"/_pipeline/recombination/${SAMPLE}_recombination_summary.tsv

ln -sfn "$RESULTS/$DAY0" "$RUNDIR/results/$DAY0"     # day-0 reference / labels (read)
ln -sfn "$NEWRES" "$RUNDIR/results/$SAMPLE"          # this sample's outputs (write; wins if SAMPLE==DAY0)

# --- input reads: the SRA download, else whatever the source run dir used ----------
rm -f "$RUNDIR/samples_dorado_basecalled/$SAMPLE".{bam,fastq,fastq.gz}
if [ -f "$SRA/$SAMPLE.fastq.gz" ]; then
  ln -sfn "$SRA/$SAMPLE.fastq.gz" "$RUNDIR/samples_dorado_basecalled/$SAMPLE.fastq.gz"
else
  src_in=""
  for ext in bam fastq fastq.gz; do
    [ -f "$SRC_RUN/samples_dorado_basecalled/$SAMPLE.$ext" ] && { src_in="$SRC_RUN/samples_dorado_basecalled/$SAMPLE.$ext"; break; }
  done
  [ -n "$src_in" ] || { echo "ERROR: no input reads for $SAMPLE"; exit 1; }
  ln -sfn "$(readlink -f "$src_in")" "$RUNDIR/samples_dorado_basecalled/$(basename "$src_in")"
fi

# --- config: curated library, variant-level IDs, spacer-only library ---------------
cd "$RUNDIR"
CURLIB="$CURLIB" SPLIB="$SPLIB" ID_LEVEL="$ID_LEVEL" python3 - <<'PY'
import os, re
p = "_pipeline/config.yaml"; s = open(p).read()
s = re.sub(r'^\s*y_prime_lib_override:.*\n', '', s, flags=re.M)
s = s.replace("references:", 'references:\n  y_prime_lib_override: "' + os.environ["CURLIB"] + '"', 1)
s = re.sub(r'^(\s*spacer_lib_dir:).*$', r'\1 "' + os.environ["SPLIB"] + '"', s, flags=re.M)
s = re.sub(r'^y_prime_id_level:.*\n', '', s, flags=re.M)
s = 'y_prime_id_level: "' + os.environ["ID_LEVEL"] + '"\n' + s
open(p, "w").write(s)
print("y_prime_lib_override ->", os.environ["CURLIB"])
print("spacer_lib_dir       ->", os.environ["SPLIB"])
print("y_prime_id_level     ->", os.environ["ID_LEVEL"])
PY

snakemake -s _pipeline/Snakefile all -c 56 --rerun-triggers mtime -n 2>&1 | grep -E "^(recombination_analyze|recombination_summary|extract_recombination_events|recombination_track_plots|all|total)" || true
snakemake -s _pipeline/Snakefile all -c 56 --rerun-triggers mtime
echo "[task ${SGE_TASK_ID}] DONE $SAMPLE end=$(date)"
echo "Outputs: $NEWRES/_pipeline/recombination/${SAMPLE}_recombination_summary.tsv"
