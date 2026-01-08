#$ -q UI-GPU
#$ -pe smp 56
#$ -l gpu=true
#$ -j y
# -o /dev/null
#$ -cwd

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# ============================================================================
# Configuration - Edit these values for your analysis
# ============================================================================

# Required inputs
BASE_NAME="dorado_7093_day0_repeat2_PromethION_no_tag_yes_rejection"
STRAIN_ID="7093"

# Paths relative to Snakemake pipeline outputs
INPUT_TSV="results/${BASE_NAME}/${BASE_NAME}_post_y_prime_probe.tsv"
READS_FASTQ="results/${BASE_NAME}/${BASE_NAME}.fastq"

# Output configuration
OUTPUT_DIR="results/${BASE_NAME}/assembly_${STRAIN_ID}"
PREFIX="assembly_${STRAIN_ID}"

# Thread configuration
THREADS=56

###### Below does NOT need to be adjusted ######

# Reference configuration
REFERENCE="/Shared/malkova_lab/Ivan/nanopore_sequencing/reference_files/6991_only_to_anchors.fasta"
ADAPTER_FILE="/Shared/malkova_lab/Ivan/nanopore_sequencing/offical_nanopore_adapter_seq+trunc.txt"

# Dorado configuration
DORADO_MODE="docker"  # "docker" or "local"
DORADO_IMAGE="/Shared/malkova_lab/dorado_files/docker_images_of_dorado/dorado-11-26-2025.sif"


# ============================================================================
# Setup
# ============================================================================

echo "Starting telomere extension and polishing pipeline"
echo "Date: $(date)"
echo "Working directory: $(pwd)"

conda activate consensus

# ============================================================================
# Step 0: Run Snakemake pipeline to generate required inputs
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 0: Running Snakemake pipeline to generate inputs"
echo "========================================================================"

# Run Snakemake until y_prime_analysis rule completes
snakemake -s Snakefile --until y_prime_analysis -c 56

# Check if Snakemake completed successfully
if [ $? -ne 0 ]; then
    echo "ERROR: Snakemake pipeline failed. Exiting."
    exit 1
fi

echo "Snakemake pipeline completed successfully"
echo "Required inputs generated:"
echo "  - TSV: ${INPUT_TSV}"
echo "  - FASTQ: ${READS_FASTQ}"

# Verify required files exist
if [ ! -f "${INPUT_TSV}" ]; then
    echo "ERROR: Required file not found: ${INPUT_TSV}"
    exit 1
fi

if [ ! -f "${READS_FASTQ}" ]; then
    echo "ERROR: Required file not found: ${READS_FASTQ}"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Define output file paths
SELECTED_READS_TEXT_FILE="${OUTPUT_DIR}/${PREFIX}_selected_read_chr_arm_pairs.txt"
SELECTED_READS_IDS_ONLY_FILE="${OUTPUT_DIR}/${PREFIX}_selected_read_ids_only.txt"
SELECTED_FASTQ="${OUTPUT_DIR}/${PREFIX}_selected_reads.fastq"
TRIMMED_FASTQ="${OUTPUT_DIR}/${PREFIX}_trimmed.fastq"
EXTENDED_REF="${OUTPUT_DIR}/${PREFIX}_extended_to_telomere_reference.fasta"
INITIAL_BAM_FILE="${OUTPUT_DIR}/${PREFIX}_initial_mapped.bam"
EXTENDED_BAM_FILE="${OUTPUT_DIR}/${PREFIX}_extended_mapped.bam"
FLYE_DIR="${OUTPUT_DIR}/flye_polish"
FLYE_REF="${OUTPUT_DIR}/${PREFIX}_flye_reference.fasta"
DORADO_ALIGNMENT="${OUTPUT_DIR}/${PREFIX}_dorado_aligned.bam"
DORADO_DIR="${OUTPUT_DIR}/dorado_polish"
DORADO_REF="${OUTPUT_DIR}/${PREFIX}_dorado_reference.fasta"
FINAL_ALIGNMENT="${OUTPUT_DIR}/${PREFIX}_final_aligment_to_dorado_ref.bam"

# ============================================================================
# Step 1: Trim adapters with porechop_abi
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 1: Trimming adapters with porechop_abi"
echo "========================================================================"

porechop_abi -i "${READS_FASTQ}" -o "${TRIMMED_FASTQ}" \
    -cap "${ADAPTER_FILE}" -t $THREADS --no_split -ddb

echo "Adapter-trimmed FASTQ saved to: ${TRIMMED_FASTQ}"

# ============================================================================
# Step 2: Select 75th percentile reads
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 2: Selecting 75th percentile reads from TSV"
echo "========================================================================"

python run_subtelomere_reference_pipeline.py select_reads \
    --input-tsv "${INPUT_TSV}" \
    --output-dir "${OUTPUT_DIR}" \
    --output-file "${SELECTED_READS_TEXT_FILE}"

echo "Selected reads saved to: ${SELECTED_READS_TEXT_FILE}"

# ============================================================================
# Step 3: Extract selected reads from FASTQ
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 3: Extracting selected reads from trimmed FASTQ"
echo "========================================================================"

python run_subtelomere_reference_pipeline.py extract_reads \
    --reads-fastq "${TRIMMED_FASTQ}" \
    --read-ids-file "${SELECTED_READS_IDS_ONLY_FILE}" \
    --output-fastq "${SELECTED_FASTQ}"

echo "Selected FASTQ saved to: ${SELECTED_FASTQ}"

# ============================================================================
# Step 4: Map reads to reference
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 4: Mapping reads to reference with minimap2"
echo "========================================================================"

minimap2 -ax map-ont -L -t "${THREADS}" "${REFERENCE}" "${SELECTED_FASTQ}" \
    | samtools sort -@ "${THREADS}" -o "${INITIAL_BAM_FILE}"

samtools index -@ "${THREADS}" "${INITIAL_BAM_FILE}"

echo "BAM file created: ${INITIAL_BAM_FILE}"

# ============================================================================
# Step 5: Extend reference using soft-clipped bases
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 5: Extending reference using soft-clipped bases"
echo "========================================================================"

python run_subtelomere_reference_pipeline.py extend_reference \
    --bamfile "${INITIAL_BAM_FILE}" \
    --reference "${REFERENCE}" \
    --read-ids-file "${SELECTED_READS_IDS_ONLY_FILE}" \
    --output-fasta "${EXTENDED_REF}"

echo "Extended reference saved to: ${EXTENDED_REF}"

# ============================================================================
# Step 6: Polish with Flye
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 6: Polishing with Flye"
echo "========================================================================"

mkdir -p "${FLYE_DIR}"

flye --polish-target "${EXTENDED_REF}" \
     --nano-hq "${SELECTED_FASTQ}" \
     --out-dir "${FLYE_DIR}" \
     --threads "${THREADS}" \
     -i 3

cp $FLYE_DIR/polished_3.fasta $FLYE_REF

echo "Flye polished reference: ${FLYE_REF}"

# ============================================================================
# Step 7: Aligning with Dorado
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 7: Aligning with Dorado"
echo "========================================================================"

if [ "${DORADO_MODE}" == "docker" ]; then

    echo "Running dorado via docker..."

    singularity exec "${DORADO_IMAGE}" \
    dorado aligner "$FLYE_REF" "$TRIMMED_FASTQ" | samtools sort -@ "$THREADS" > \
    "$DORADO_ALIGNMENT"

else

    echo "Running dorado locally..."

    dorado aligner "$FLYE_REF" "$TRIMMED_FASTQ" | samtools sort -@ "$THREADS" > \
    "$DORADO_ALIGNMENT"

fi

samtools index -@ "$THREADS" "$DORADO_ALIGNMENT"

echo "Dorado aligned: ${DORADO_ALIGNMENT}"

# ============================================================================
# Step 8: Polishing with Dorado
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 8: Polishing with Dorado"
echo "========================================================================"

echo "Polishing: ${DORADO_ALIGNMENT}"

if [[ "$DORADO_MODE" == "docker" ]]; then

    echo "Running dorado via docker..."

    singularity exec --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True "${DORADO_IMAGE}" \
    dorado polish "${DORADO_ALIGNMENT}" "$FLYE_REF" \
    --threads $THREADS --batchsize 4 --draft-batchsize 500K \
    --bam-chunk 20000 --bam-subchunk 2000 \
    --window-len 7000 --window-overlap 3000 \
    --ignore-read-groups --vcf --min-mapq 20 --min-depth 5 \
    --output-dir "${DORADO_DIR}"

elif [[ "$DORADO_MODE" == "local" ]]; then

    echo "Running dorado locally..."

    dorado polish "${DORADO_ALIGNMENT}" "$FLYE_REF" \
    --threads $THREADS --batchsize 4 --draft-batchsize 500K \
    --bam-chunk 20000 --bam-subchunk 2000 \
    --window-len 7000 --window-overlap 3000 \
    --ignore-read-groups --vcf --min-mapq 20 --min-depth 5 \
    --output-dir "${DORADO_DIR}"

else
    echo "ERROR: unknown mode '$DORADO_MODE'"
    exit 1
fi

cp $DORADO_DIR/consensus.fasta $DORADO_REF

echo "Dorado polished reference saved to: ${DORADO_REF}"


# ============================================================================
# Step 9: Final ALignment
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 9: Final ALignment"
echo "========================================================================"

minimap2 -ax map-ont -L -t "${THREADS}" "${DORADO_REF}" "${SELECTED_FASTQ}" \
    | samtools sort -@ "${THREADS}" -o "${FINAL_ALIGNMENT}"

samtools index -@ "${THREADS}" "${FINAL_ALIGNMENT}"

echo "BAM file created: ${FINAL_ALIGNMENT}"

# ============================================================================
# Pipeline Complete
# ============================================================================

echo ""
echo "========================================================================"
echo "Pipeline completed successfully!"
echo "========================================================================"
echo "Final outputs:"
echo "  - Selected reads: ${SELECTED_READS}"
echo "  - Extended reference: ${EXTENDED_REF}"
echo "  - Flye polished reference: ${FLYE_REF}"
echo "  - Dorado polished reference: ${DORADO_REF}"
echo "  - Final alignment to ${DORADO_REF}:"
echo "          ${FINAL_ALIGNMENT}:"
echo ""
echo "End time: $(date)"

