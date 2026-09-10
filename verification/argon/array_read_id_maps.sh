#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 2
#$ -j y
#$ -cwd
#$ -N readid_map
#$ -t 1-27
# Write results/<sample>/_pipeline/<sample>_read_id_map.tsv (pipeline read_id -> original ONT UUID)
# from each sample's raw FASTQ.
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate consensus
RUN_ROOT=/nfsscratch/amalkova/telo_sra_runs; RESULTS=$RUN_ROOT/results; STAGE=$RUN_ROOT/v2_stage
SAMPLE=$(sed -n "${SGE_TASK_ID}p" "$STAGE/read_id_map_manifest.txt")
[ -n "$SAMPLE" ] || { echo "no sample for task $SGE_TASK_ID"; exit 1; }
P=$RESULTS/$SAMPLE/_pipeline
src=$P/${SAMPLE}_raw.fastq; [ -f "$src" ] || src=$P/${SAMPLE}.fastq
[ -f "$src" ] || { echo "no fastq for $SAMPLE"; exit 1; }
python "$STAGE/make_read_id_map.py" "$src" "$P/${SAMPLE}_read_id_map.tsv"
echo "[task $SGE_TASK_ID] DONE $SAMPLE"
