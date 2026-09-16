#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd
#$ -N elemYPfix
#$ -t 1-8
# As array_yprime_element_lib.sh, but against element libraries rebuilt from BOUNDARY-CORRECTED
# BEDs: chr16L_Y_Prime_1 is trimmed by the 77 bp it over-extends into non-Y' sequence (measured
# per reference by verification/fix_yprime_boundary.py, 10-17 partners agreeing on every one).
# The ONLY difference from the elemYP run is those 77 bp, so any change in the read-to-element
# matching is attributable to that boundary alone.
#
# Re-run the recombination step for the eight 6991 day-0 populations against an
# ELEMENT-LEVEL Y' library: every reference Y' element is its own entry with its own ID.
# The point is not the recombination call (meaningless at this resolution) but the
# per-copy assignment it records -- for each Y' copy on a read, the single reference
# element RepeatMasker matched best. Every grouping scheme is then scored offline as a
# relabelling of that one fixed assignment, so the schemes are compared on identical data.
#
# Nothing existing is overwritten: outputs go to results/<sample>__elemYP, whose steps 0-6
# are symlinks (cp -as) into the sample's existing outputs.
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

TAG=elemYPfix
RUN_ROOT=/nfsscratch/amalkova/telo_sra_runs
RESULTS=$RUN_ROOT/results
STAGE=$RUN_ROOT/v2_stage
CODE=$STAGE/path
MANIFEST=$STAGE/elem_lib_manifest.txt

read -r SAMPLE DAY0 < <(sed -n "${SGE_TASK_ID}p" "$MANIFEST")
[ -n "${SAMPLE:-}" ] || { echo "ERROR: no sample for task ${SGE_TASK_ID}"; exit 1; }
SRC_RUN=$RUN_ROOT/$SAMPLE
RUNDIR=$RUN_ROOT/${SAMPLE}__$TAG
NEWRES=$RESULTS/${SAMPLE}__$TAG
D0P=$RESULTS/$DAY0/_pipeline
ASM=$D0P/assembly_$DAY0/assembly_${DAY0}_dorado_reference.fasta
BED=$D0P/pretelomeric_labels/pretelomeric_regions_${DAY0}_simp.bed
SPLIB=$STAGE/spacer_variable/$DAY0/pairings_for_spacers/${DAY0}_pairings
ELEMDIR=$STAGE/element_libs_fixed
ELEMLIB=$ELEMDIR/elem_${DAY0}.fasta
echo "[task ${SGE_TASK_ID}] SAMPLE=$SAMPLE DAY0=$DAY0 host=$(hostname) start=$(date)"
[ -f "$ASM" ] || { echo "ERROR: missing assembly $ASM"; exit 1; }
[ -f "$BED" ] || { echo "ERROR: missing bed $BED"; exit 1; }
[ -d "$SPLIB" ] || { echo "ERROR: missing spacer library $SPLIB"; exit 1; }

# --- element library for this day-0 assembly (built once, per task-safe) -----------
# The corrected libraries are built locally from the boundary-fixed BEDs and staged; rebuilding
# here would silently regenerate them from the UNcorrected BED and undo the fix.
[ -s "$ELEMLIB" ] || { echo "ERROR: corrected element library not staged: $ELEMLIB"; exit 1; }
echo "element library: $ELEMLIB ($(grep -c '^>' "$ELEMLIB") entries)"

# --- run dir with the current code -------------------------------------------------
mkdir -p "$RUNDIR"
rsync -a --exclude "results/" --exclude "samples_dorado_basecalled/" --exclude ".snakemake/" "$SRC_RUN/" "$RUNDIR/"
for f in analyze_features.py yprime_path.py aggregate_recombination.py recombination_utils.py; do
  cp "$CODE/$f" "$RUNDIR/_pipeline/scripts/$f"
done
python3 "$CODE/patch_snakefile_v2.py" "$RUNDIR/_pipeline/Snakefile"
mkdir -p "$RUNDIR/results" "$RUNDIR/samples_dorado_basecalled"

# --- results dir: symlink farm over the sample's existing outputs -------------------
if [ ! -f "$NEWRES/.farm_ready" ]; then
  rm -rf "$NEWRES" "$NEWRES.tmp"
  mkdir -p "$NEWRES.tmp/_pipeline"
  cp -as "$RESULTS/$SAMPLE/_pipeline/." "$NEWRES.tmp/_pipeline/"
  rm -rf "$NEWRES.tmp"/_pipeline/recombination_v1 "$NEWRES.tmp"/_pipeline/recombination_v2a \
         "$NEWRES.tmp"/_pipeline/recombination_v2b
  mkdir -p "$NEWRES"; mv "$NEWRES.tmp/_pipeline" "$NEWRES/_pipeline"; rmdir "$NEWRES.tmp"
  touch "$NEWRES/.farm_ready"
fi
rm -rf "$NEWRES/_pipeline/recombination_events" "$NEWRES/_pipeline/graphs/recombination_tracks"
rm -f "$NEWRES"/_pipeline/recombination/${SAMPLE}_chr*_features.tsv \
      "$NEWRES"/_pipeline/recombination/${SAMPLE}_chr*_features.tsv.skipped \
      "$NEWRES"/_pipeline/recombination/${SAMPLE}_recombination_summary.tsv

ln -sfn "$RESULTS/$DAY0" "$RUNDIR/results/$DAY0"
ln -sfn "$NEWRES" "$RUNDIR/results/$SAMPLE"
rm -f "$RUNDIR/samples_dorado_basecalled/$SAMPLE".{bam,fastq,fastq.gz}
src_in=""
for ext in bam fastq fastq.gz; do
  [ -f "$SRC_RUN/samples_dorado_basecalled/$SAMPLE.$ext" ] && { src_in="$SRC_RUN/samples_dorado_basecalled/$SAMPLE.$ext"; break; }
done
[ -n "$src_in" ] || { echo "ERROR: no input reads for $SAMPLE"; exit 1; }
ln -sfn "$(readlink -f "$src_in")" "$RUNDIR/samples_dorado_basecalled/$(basename "$src_in")"

cd "$RUNDIR"
ELEMLIB="$ELEMLIB" SPLIB="$SPLIB" python3 - <<'PY'
import os, re
p = "_pipeline/config.yaml"; s = open(p).read()
s = re.sub(r'^\s*y_prime_lib_override:.*\n', '', s, flags=re.M)
s = s.replace("references:", 'references:\n  y_prime_lib_override: "' + os.environ["ELEMLIB"] + '"', 1)
s = re.sub(r'^(\s*spacer_lib_dir:).*$', r'\1 "' + os.environ["SPLIB"] + '"', s, flags=re.M)
s = re.sub(r'^y_prime_id_level:.*\n', '', s, flags=re.M)
s = 'y_prime_id_level: "variant"\n' + s
open(p, "w").write(s)
print("y_prime_lib_override ->", os.environ["ELEMLIB"])
PY

snakemake -s _pipeline/Snakefile all -c 56 --rerun-triggers mtime
echo "[task ${SGE_TASK_ID}] DONE $SAMPLE end=$(date)"
