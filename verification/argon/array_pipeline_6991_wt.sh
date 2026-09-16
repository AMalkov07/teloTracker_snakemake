#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 56
#$ -j y
#$ -cwd
#$ -N wt6991
#$ -t 1-8
# Full pipeline (steps 0-6) plus the recombination step for the eight 6991 (WT) day-4 / day-5
# populations, against a 6991 day-0 reference, using the CURATED 6991 Y' library at variant
# level and the spacer-only library -- the same settings as the curated re-runs, so the output
# is directly comparable with Supplementary Data 5.
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

DAY0=${DAY0:-6991_day0_with_selection_repeat}   # cleanest 6991 reference (Part A: counts exact, no major feature diffs)
STRAIN=6991
RUN_ROOT=/nfsscratch/amalkova/telo_sra_runs
RESULTS=$RUN_ROOT/results
STAGE=$RUN_ROOT/v2_stage
CODE=$STAGE/path
CODE_SRC=$RUN_ROOT/7172_day0_with_selection          # a complete run-dir copy (run_pipeline.py + _pipeline)
SRA=/nfsscratch/amalkova/sra_jobs/outputs
MANIFEST=$STAGE/wt6991_manifest.txt

SAMPLE=$(sed -n "${SGE_TASK_ID}p" "$MANIFEST")
[ -n "${SAMPLE:-}" ] || { echo "ERROR: no sample for task ${SGE_TASK_ID}"; exit 1; }
CURLIB=$STAGE/yprime_lib_curated_${STRAIN}.fasta
SPLIB=$STAGE/spacer_variable/$DAY0/pairings_for_spacers/${DAY0}_pairings
RUNDIR=$RUN_ROOT/$SAMPLE
echo "[task ${SGE_TASK_ID}] SAMPLE=$SAMPLE DAY0=$DAY0 host=$(hostname) start=$(date)"
[ -f "$CURLIB" ] || { echo "ERROR: missing $CURLIB"; exit 1; }
[ -d "$SPLIB" ]  || { echo "ERROR: missing $SPLIB"; exit 1; }
[ -s "$SRA/$SAMPLE.fastq.gz" ] || { echo "ERROR: no reads at $SRA/$SAMPLE.fastq.gz"; exit 1; }

# run dir with the current code
mkdir -p "$RUNDIR"
rsync -a --exclude "results/" --exclude "samples_dorado_basecalled/" --exclude ".snakemake/" "$CODE_SRC/" "$RUNDIR/"
for f in analyze_features.py yprime_path.py aggregate_recombination.py recombination_utils.py; do
  cp "$CODE/$f" "$RUNDIR/_pipeline/scripts/$f"
done
python3 "$CODE/patch_snakefile_v2.py" "$RUNDIR/_pipeline/Snakefile"
mkdir -p "$RUNDIR/results" "$RUNDIR/samples_dorado_basecalled" "$RESULTS/$SAMPLE"
ln -sfn "$RESULTS/$DAY0" "$RUNDIR/results/$DAY0"         # day-0 reference + labels (read)
ln -sfn "$RESULTS/$SAMPLE" "$RUNDIR/results/$SAMPLE"     # this sample's outputs (write)
ln -sfn "$SRA/$SAMPLE.fastq.gz" "$RUNDIR/samples_dorado_basecalled/$SAMPLE.fastq.gz"

cd "$RUNDIR"
# config: base = this sample, strain/day0 = the 6991 day-0 reference
SAMPLE="$SAMPLE" DAY0="$DAY0" python3 -c "import os,sys; sys.path.insert(0,'.'); import run_pipeline as r; r.write_snakemake_config({'base_name':os.environ['SAMPLE'],'strain':os.environ['DAY0'],'day0_base_name':os.environ['DAY0']})"
DAY0="$DAY0" CURLIB="$CURLIB" SPLIB="$SPLIB" python3 - <<'PY'
import os, re
p="_pipeline/config.yaml"; s=open(p).read(); d=os.environ["DAY0"]
ref=f"results/{d}/_pipeline/assembly_{d}/assembly_{d}_dorado_reference.fasta"
bed=f"results/{d}/_pipeline/pretelomeric_labels/pretelomeric_regions_{d}_simp.bed"
if "day0_ref:" not in s:
    s=s.replace("references:", 'references:\n  day0_ref: "'+ref+'"\n  day0_bed: "'+bed+'"', 1)
s=re.sub(r'^\s*y_prime_lib_override:.*\n','',s,flags=re.M)
s=s.replace("references:", 'references:\n  y_prime_lib_override: "'+os.environ["CURLIB"]+'"', 1)
s=re.sub(r'^(\s*spacer_lib_dir:).*$', r'\1 "'+os.environ["SPLIB"]+'"', s, flags=re.M)
s=re.sub(r'^y_prime_id_level:.*\n','',s,flags=re.M)
s='y_prime_id_level: "variant"\n'+s
open(p,"w").write(s)
print("day0_ref/bed, curated library, spacer-only library and variant IDs set")
PY
grep -E "base_name|strain|day0_ref|y_prime_lib_override|spacer_lib_dir|y_prime_id_level" _pipeline/config.yaml

snakemake -s _pipeline/Snakefile all -c 56
P=$RESULTS/$SAMPLE/_pipeline
src=$P/${SAMPLE}_raw.fastq; [ -f "$src" ] || src=$P/${SAMPLE}.fastq
[ -f "$src" ] && python3 "$STAGE/make_read_id_map.py" "$src" "$P/${SAMPLE}_read_id_map.tsv" || echo "WARN: no raw fastq for the read-ID map"
echo "[task ${SGE_TASK_ID}] DONE $SAMPLE end=$(date)"
