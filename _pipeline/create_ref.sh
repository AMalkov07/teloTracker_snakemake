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
BASE_NAME="dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection"
STRAIN_ID="7575"
ANCHOR_SET="telomerase_shutoff_anchors"

# Paths relative to Snakemake pipeline outputs
INPUT_TSV="results/${BASE_NAME}/_pipeline/${BASE_NAME}_post_y_prime_probe.tsv"
READS_FASTQ="results/${BASE_NAME}/_pipeline/${BASE_NAME}.fastq"
READS_FASTA="results/${BASE_NAME}/_pipeline/${BASE_NAME}.fasta"
ALL_MATCHES_TSV="results/${BASE_NAME}/_pipeline/blast/all_matches_${BASE_NAME}_blasted_${ANCHOR_SET}.tsv"

# Output configuration
OUTPUT_DIR="results/${BASE_NAME}/_pipeline/assembly_${STRAIN_ID}"
PREFIX="assembly_${STRAIN_ID}"

# Thread configuration
THREADS=56

# Polishing scope: "genome" polishes the whole extended reference with every read.
# "ends" polishes each arm (its extension + POLISH_END_FLANK bp of the base contig) with only
# the reads anchored to that arm, and leaves the base contig interiors unpolished. Use "ends"
# for end-enriched libraries: interiors at a few x are corrupted by genome-wide polishing.
POLISH_SCOPE="genome"
POLISH_END_FLANK=5000

###### Below does NOT need to be adjusted ######

# Reference configuration (paths relative to the repo root)
REFERENCE="_pipeline/references/6991_only_to_anchors.fasta"
ADAPTER_FILE="_pipeline/references/offical_nanopore_adapter_seq+trunc.txt"
SCRIPTS_DIR="_pipeline/scripts"

# Dorado configuration
DORADO_MODE="docker"  # "docker" or "local"
DORADO_IMAGE=""  # Singularity .sif with dorado (docker mode); run_pipeline.py sets it from dorado_image
DORADO_MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.2.0"  # Basecaller model used for original reads


# ============================================================================
# Setup
# ============================================================================

echo "Starting telomere extension and polishing pipeline"
echo "Date: $(date)"
echo "Working directory: $(pwd)"

# Polishing needs dorado; stop now rather than after the read processing and assembly.
if [ "${DORADO_MODE}" = "docker" ] && [ ! -f "${DORADO_IMAGE}" ]; then
    echo "ERROR: DORADO_MODE is docker but the dorado image is not set or not found: '${DORADO_IMAGE}'"
    echo "       Set dorado_image in pipeline_config.yaml to a Singularity .sif with dorado."
    exit 1
fi

conda activate consensus

# ============================================================================
# Step 0: Run Snakemake pipeline to generate required inputs
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 0: Running Snakemake pipeline to generate inputs"
echo "========================================================================"

# Run Snakemake until y_prime_analysis rule completes
snakemake -s _pipeline/Snakefile through_y_prime_analysis -c "${THREADS}"

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
# Note: The Python script creates the ids-only file by appending _only before .txt
SELECTED_READS_IDS_ONLY_FILE="${OUTPUT_DIR}/${PREFIX}_selected_read_chr_arm_pairs_only.txt"
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

# The scaffold is chosen by consensus among the reads bracketing the 75th
# percentile, not by trusting that single read. Needs random access to the reads
# (.fai) and the anchor positions, to compare the region each candidate would
# actually graft onto the reference.
if [ ! -f "${READS_FASTA}.fai" ]; then
    echo "Indexing ${READS_FASTA} for consensus scaffold selection"
    samtools faidx "${READS_FASTA}"
fi

python "${SCRIPTS_DIR}/run_subtelomere_reference_pipeline.py" select_reads \
    --input-tsv "${INPUT_TSV}" \
    --all-matches-tsv "${ALL_MATCHES_TSV}" \
    --reads-fasta "${READS_FASTA}" \
    --output-dir "${OUTPUT_DIR}" \
    --output-file "${SELECTED_READS_TEXT_FILE}" \
    --report-tsv "${OUTPUT_DIR}/${PREFIX}_scaffold_selection_report.tsv" \
    --threads "${THREADS}"

echo "Selected reads saved to: ${SELECTED_READS_TEXT_FILE}"

# ============================================================================
# Step 3: Extract selected reads from FASTQ
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 3: Extracting selected reads from trimmed FASTQ"
echo "========================================================================"

python "${SCRIPTS_DIR}/run_subtelomere_reference_pipeline.py" extract_reads \
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

# -Y soft-clips supplementary alignments so they keep the full read sequence.
# extend_reference_multi selects a scaffold's clip by contig+end rather than by
# the primary flag, because minimap2 scores by alignment LENGTH and these reads
# are mostly shared subtelomeric repeat -- the correct placement is often
# demoted to supplementary (7858 9L: chr9 @ 99.8% lost to chr15 @ 90.0%).
# Without -Y those records are hard-clipped and their clip bases are gone.
minimap2 -ax map-ont -Y -L -t "${THREADS}" "${REFERENCE}" "${SELECTED_FASTQ}" \
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

python "${SCRIPTS_DIR}/run_subtelomere_reference_pipeline.py" extend_reference \
    --bamfile "${INITIAL_BAM_FILE}" \
    --reference "${REFERENCE}" \
    --read-ids-file "${SELECTED_READS_IDS_ONLY_FILE}" \
    --chr-arm-pairs-file "${SELECTED_READS_TEXT_FILE}" \
    --output-fasta "${EXTENDED_REF}"

# Rebuild the .fai now, since the assembly dir is reused across runs and htslib
# trusts a pre-existing index rather than regenerating it. A stale index from a
# previous run points at wrong byte offsets, and any htslib random-access reader
# (dorado polish) then fetches garbage. dorado aligner builds its own minimap2
# index so it is unaffected, which is why this stayed latent until extension
# stopped failing earlier.
rm -f "${EXTENDED_REF}.fai"
samtools faidx "${EXTENDED_REF}"

echo "Extended reference saved to: ${EXTENDED_REF}"

# Write a BAM's header with the basecaller model on its @RG line; dorado polish
# reads it to choose its polishing model.
add_model_to_header() {  # in.bam out.bam
    local hdr="${2%.bam}.header.sam"
    samtools view -H "$1" > "$hdr"
    if grep -q "^@RG" "$hdr"; then
        sed -i "s/^@RG\t/@RG\tDS:basecall_model=${DORADO_MODEL}\t/" "$hdr"
    else
        echo -e "@RG\tID:unknown\tDS:basecall_model=${DORADO_MODEL}" >> "$hdr"
    fi
    samtools reheader "$hdr" "$1" > "$2"
    samtools index "$2"
    rm "$hdr"
}

run_dorado() {  # dorado subcommand and arguments
    if [[ "$DORADO_MODE" == "docker" ]]; then
        singularity exec --nv --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True "${DORADO_IMAGE}" dorado "$@"
    elif [[ "$DORADO_MODE" == "local" ]]; then
        dorado "$@"
    else
        echo "ERROR: unknown mode '$DORADO_MODE'"; return 1
    fi
}

if [[ "$POLISH_SCOPE" == "ends" ]]; then

# ============================================================================
# Steps 6-8 (ends): polish each arm with its own anchored reads
# ============================================================================

echo ""
echo "========================================================================"
echo "Steps 6-8: Polishing each arm with the reads anchored to it (POLISH_SCOPE=ends)"
echo "========================================================================"

ENDS_DIR="${OUTPUT_DIR}/end_polish"
rm -rf "${ENDS_DIR}"

python "${SCRIPTS_DIR}/polish_ends.py" prepare \
    --extended "${EXTENDED_REF}" --base "${REFERENCE}" --tsv "${INPUT_TSV}" \
    --reads "${TRIMMED_FASTQ}" --flank "${POLISH_END_FLANK}" --outdir "${ENDS_DIR}"

while IFS=$'\t' read -r END_ID CONTIG SIDE EXT_BP FLANK_BP SEG_BP N_READS DO_POLISH; do
    [[ "$DO_POLISH" == "yes" ]] || continue
    D="${ENDS_DIR}/${END_ID}"
    echo "--- ${END_ID}: ${N_READS} anchored reads, ${SEG_BP} bp segment"
    if flye --polish-target "$D/segment.fasta" --nano-hq "$D/reads.fastq" \
            --out-dir "$D/flye" --threads "${THREADS}" -i 3 > "$D/flye.stdout" 2>&1; then
        cp "$D/flye/polished_3.fasta" "$D/flye.fasta"
    else
        echo "WARNING: ${END_ID} Flye polish failed (see $D/flye.stdout); arm left unpolished"
        continue
    fi
    if run_dorado aligner "$D/flye.fasta" "$D/reads.fastq" | samtools sort -@ "$THREADS" > "$D/aligned.bam" \
       && add_model_to_header "$D/aligned.bam" "$D/aligned_with_model.bam" \
       && run_dorado polish "$D/aligned_with_model.bam" "$D/flye.fasta" \
              --threads $THREADS --batchsize 4 --draft-batchsize 500K \
              --bam-chunk 20000 --bam-subchunk 2000 \
              --window-len 7000 --window-overlap 3000 \
              --ignore-read-groups --vcf --min-mapq 20 --min-depth 5 \
              --output-dir "$D/dorado" > "$D/dorado.stdout" 2>&1; then
        cp "$D/dorado/consensus.fasta" "$D/dorado.fasta"
    else
        echo "WARNING: ${END_ID} dorado polish failed (see $D/dorado.stdout); arm keeps its Flye polish"
        cp "$D/flye.fasta" "$D/dorado.fasta"
    fi
done < <(tail -n +2 "${ENDS_DIR}/manifest.tsv")

python "${SCRIPTS_DIR}/polish_ends.py" splice --extended "${EXTENDED_REF}" \
    --outdir "${ENDS_DIR}" --polished-name flye.fasta --output "${FLYE_REF}"
python "${SCRIPTS_DIR}/polish_ends.py" splice --extended "${EXTENDED_REF}" \
    --outdir "${ENDS_DIR}" --polished-name dorado.fasta --output "${DORADO_REF}"

rm -f "${FLYE_REF}.fai" "${DORADO_REF}.fai"
samtools faidx "${FLYE_REF}"
samtools faidx "${DORADO_REF}"

echo "Flye polished reference (arms only): ${FLYE_REF}"
echo "Dorado polished reference (arms only) saved to: ${DORADO_REF}"

elif [[ "$POLISH_SCOPE" == "genome" ]]; then

# ============================================================================
# Step 6: Polish with Flye
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 6: Polishing with Flye"
echo "========================================================================"

mkdir -p "${FLYE_DIR}"

flye --polish-target "${EXTENDED_REF}" \
     --nano-hq "${TRIMMED_FASTQ}" \
     --out-dir "${FLYE_DIR}" \
     --threads "${THREADS}" \
     -i 3

cp $FLYE_DIR/polished_3.fasta $FLYE_REF

# Same staleness hazard as EXTENDED_REF above: dorado polish reads FLYE_REF by
# faidx random access, so a leftover index from a prior run is what aborted the
# reruns ("Could not fetch sequence for chr9_extended").
rm -f "${FLYE_REF}.fai"
samtools faidx "${FLYE_REF}"

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

# First, we need to add the basecaller model to the BAM header
# dorado polish requires this information to select the correct polishing model
DORADO_ALIGNMENT_WITH_MODEL="${OUTPUT_DIR}/${PREFIX}_dorado_aligned_with_model.bam"

echo "Adding basecaller model to BAM header..."
echo "Model: ${DORADO_MODEL}"

# Extract header, add DS tag with basecaller model, and reheader the BAM
samtools view -H "${DORADO_ALIGNMENT}" > "${OUTPUT_DIR}/header.sam"

# Check if there's already a @RG line, if not add one with the model info
if grep -q "^@RG" "${OUTPUT_DIR}/header.sam"; then
    # Add DS tag to existing @RG line(s)
    sed -i "s/^@RG\t/@RG\tDS:basecall_model=${DORADO_MODEL}\t/" "${OUTPUT_DIR}/header.sam"
else
    # Add a new @RG line with the model info
    echo -e "@RG\tID:unknown\tDS:basecall_model=${DORADO_MODEL}" >> "${OUTPUT_DIR}/header.sam"
fi

samtools reheader "${OUTPUT_DIR}/header.sam" "${DORADO_ALIGNMENT}" > "${DORADO_ALIGNMENT_WITH_MODEL}"
samtools index "${DORADO_ALIGNMENT_WITH_MODEL}"
rm "${OUTPUT_DIR}/header.sam"

echo "Created BAM with model info: ${DORADO_ALIGNMENT_WITH_MODEL}"

mkdir -p "${DORADO_DIR}"

if [[ "$DORADO_MODE" == "docker" ]]; then

    echo "Running dorado via docker..."

    singularity exec --nv --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True "${DORADO_IMAGE}" \
    dorado polish "${DORADO_ALIGNMENT_WITH_MODEL}" "$FLYE_REF" \
    --threads $THREADS --batchsize 4 --draft-batchsize 500K \
    --bam-chunk 20000 --bam-subchunk 2000 \
    --window-len 7000 --window-overlap 3000 \
    --ignore-read-groups --vcf --min-mapq 20 --min-depth 5 \
    --output-dir "${DORADO_DIR}"

elif [[ "$DORADO_MODE" == "local" ]]; then

    echo "Running dorado locally..."

    dorado polish "${DORADO_ALIGNMENT_WITH_MODEL}" "$FLYE_REF" \
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

# Final reference; label_regions and recombination read it by faidx random
# access, so give it a fresh index for the same reason as above.
rm -f "${DORADO_REF}.fai"
samtools faidx "${DORADO_REF}"

echo "Dorado polished reference saved to: ${DORADO_REF}"

# Genome-wide polishing at low depth deletes sequence rather than correcting it; a large
# length change means the library does not cover the interiors (use POLISH_SCOPE="ends").
awk -v e="$(awk '{s+=$2} END{print s}' "${EXTENDED_REF}.fai")" '{s+=$2} END{
    d = (s - e) / e * 100
    printf "Total length: %d bp extended -> %d bp polished (%+.2f%%)\n", e, s, d
    if (d < -2 || d > 2) print "WARNING: polishing changed the genome length by more than 2%; " \
        "interiors are probably under-covered. Consider POLISH_SCOPE=\"ends\"."
}' "${DORADO_REF}.fai"

# The genome total hides losses on single chromosomes (the 7172 / 7302 rebuilds lost 188 kb
# and 83 kb of interior sequence, under 2 % of the genome), so compare every contig.
# Flag any that changed by more than 1 % or 10 kb; the table is kept next to the reference.
POLISH_CHECK_TSV="${OUTPUT_DIR}/polish_length_check.tsv"
awk -F'\t' -v OFS='\t' -v out="${POLISH_CHECK_TSV}" '
    NR == FNR { ext[$1] = $2; next }
    {
        n_pol++
        if (!($1 in ext)) { unmatched++; next }
        d = $2 - ext[$1]; pct = d / ext[$1] * 100
        flag = (pct < -1 || pct > 1 || d < -10000 || d > 10000) ? "CHANGED" : "ok"
        if (!header++) print "contig", "extended_bp", "polished_bp", "change_bp", "change_pct", "flag" > out
        printf "%s\t%d\t%d\t%+d\t%+.2f\t%s\n", $1, ext[$1], $2, d, pct, flag > out
        if (flag == "CHANGED") bad[++nbad] = sprintf("  %-20s %12d -> %12d bp (%+d bp, %+.2f%%)", $1, ext[$1], $2, d, pct)
    }
    END {
        if (unmatched == n_pol) {
            print "NOTE: polished contig names do not match the extended reference; per-contig check skipped."
            exit
        }
        if (unmatched) printf "NOTE: %d polished contig(s) have no counterpart in the extended reference.\n", unmatched
        if (nbad) {
            print "========================================================================"
            printf "WARNING: polishing changed the length of %d contig(s) by more than 1%% or 10 kb:\n", nbad
            for (i = 1; i <= nbad; i++) print bad[i]
            print "Interiors are probably under-covered; consider POLISH_SCOPE=\"ends\"."
            print "========================================================================"
        } else {
            print "Per-contig length check: no contig changed by more than 1% or 10 kb."
        }
        print "Per-contig lengths: " out
    }' "${EXTENDED_REF}.fai" "${DORADO_REF}.fai"

else
    echo "ERROR: unknown POLISH_SCOPE '${POLISH_SCOPE}' (use \"genome\" or \"ends\")"
    exit 1
fi


# ============================================================================
# Step 9: Final ALignment
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 9: Final ALignment"
echo "========================================================================"

if [ "${DORADO_MODE}" == "docker" ]; then

    echo "Running dorado via docker..."

    singularity exec "${DORADO_IMAGE}" \
    dorado aligner "$DORADO_REF" "$TRIMMED_FASTQ" | samtools sort -@ "$THREADS" > \
    "$FINAL_ALIGNMENT"

else

    echo "Running dorado locally..."

    dorado aligner "$DORADO_REF" "$TRIMMED_FASTQ" | samtools sort -@ "$THREADS" > \
    "$FINAL_ALIGNMENT"

fi

samtools index -@ "$THREADS" "$FINAL_ALIGNMENT"

echo "Dorado aligned: ${FINAL_ALIGNMENT}"

# ============================================================================
# Pipeline Complete
# ============================================================================

echo ""
echo "========================================================================"
echo "Pipeline completed successfully!"
echo "========================================================================"
echo "Final outputs:"
echo "  - Selected reads: ${SELECTED_FASTQ}"
echo "  - Extended reference: ${EXTENDED_REF}"
echo "  - Flye polished reference: ${FLYE_REF}"
echo "  - Dorado polished reference: ${DORADO_REF}"
echo "  - Final alignment to ${DORADO_REF}:"
echo "          ${FINAL_ALIGNMENT}:"
echo ""
echo "End time: $(date)"