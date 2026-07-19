#!/bin/bash
#$ -q TELOMERE2,UI
#$ -pe smp 56
#$ -j y
#$ -cwd

conda activate consensus

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# ============================================================================
# Configuration - Edit these values for your reference
# ============================================================================

BASE_NAME="dorado_7302_day0_PromethION_no_tag_yes_rejection"
STRAIN_ID="7302"

# Input: Reference FASTA file
REFERENCE_DIR="results/${BASE_NAME}/_pipeline/assembly_${STRAIN_ID}"
#REFERENCE_DIR="labeling_test"
REFERENCE_FASTA="${REFERENCE_DIR}/assembly_${STRAIN_ID}_dorado_reference.fasta"
#REFERENCE_FASTA="${REFERENCE_DIR}/${STRAIN_ID}.fasta"

# Reference sequences for labeling (from strain 6991)
REFERENCES_DIR="_pipeline/references"
ANCHORS_FASTA="${REFERENCES_DIR}/test_anchors.fasta"
YPRIMES_FASTA="${REFERENCES_DIR}/repeatmasker_6991_all_y_primes.fasta"
XPRIMES_FASTA="${REFERENCES_DIR}/6991_xprimes.fasta"  # X prime sequences for detection
PROBE_FASTA="${REFERENCES_DIR}/probe.fasta"  # Y prime probe for verification

# Scripts directory
SCRIPTS_DIR="_pipeline/scripts"

# Output configuration
OUTPUT_DIR="results/${BASE_NAME}/_pipeline/pretelomeric_labels"
PREFIX="pretelomeric_regions_${STRAIN_ID}"

# Thread configuration
THREADS=56

# Y prime clustering configuration (consumed by cluster_yprimes_paper_method.py)
#   YPRIME_LINKAGE:   'complete' | 'average' | 'single'
#     complete = merge clusters only when every cross-cluster pair is >= threshold (strict)
#     average  = merge when the mean cross-cluster pair is >= threshold (default)
#     single   = merge if any cross-cluster pair is >= threshold (most permissive)
#   YPRIME_STOP_MODE: 'silhouette' | 'threshold'
#     silhouette = data-driven, picks the k with the highest silhouette score (default)
#     threshold  = fixed 97% identity cutoff
YPRIME_LINKAGE="average"
YPRIME_STOP_MODE="silhouette"

# BLAST parameters
MIN_PIDENT=75.0        # Minimum percent identity (lowered for cross-strain comparison)
MIN_LENGTH=100         # Minimum alignment length
EVALUE=1e-5           # E-value threshold (relaxed for cross-strain)

# Boundary adjustment parameters (to maximize ITS regions between features)
# Trims telomeric bases from feature boundaries:
#   - R arm: trims T and G from feature ends until hitting A or C
#   - L arm: trims A and C from feature ends until hitting T or G
ADJUST_BOUNDARIES=true  # Set to true to trim telomeric bases from feature boundaries
BOUNDARY_WINDOW=50      # Minimum feature size to maintain after trimming (bp)

# Debugging options
DEBUG_BOUNDARIES=false  # Set to true to enable detailed boundary adjustment debugging

# ============================================================================
# Setup
# ============================================================================

echo "Starting pre-telomeric region labeling pipeline"
echo "Date: $(date)"
echo "Working directory: $(pwd)"

# Verify input files exist
if [ ! -f "${REFERENCE_FASTA}" ]; then
    echo "ERROR: Reference FASTA not found: ${REFERENCE_FASTA}"
    echo "Please run create_ref.sh first to generate the reference"
    exit 1
fi

if [ ! -f "${ANCHORS_FASTA}" ]; then
    echo "ERROR: Anchors FASTA not found: ${ANCHORS_FASTA}"
    exit 1
fi

if [ ! -f "${YPRIMES_FASTA}" ]; then
    echo "ERROR: Y primes FASTA not found: ${YPRIMES_FASTA}"
    exit 1
fi

if [ ! -f "${PROBE_FASTA}" ]; then
    echo "WARNING: Probe FASTA not found: ${PROBE_FASTA}"
    echo "Y prime verification will be skipped"
    PROBE_FASTA=""
fi

if [ ! -f "${XPRIMES_FASTA}" ]; then
    echo "WARNING: X primes FASTA not found: ${XPRIMES_FASTA}"
    echo "X prime detection will be skipped"
    XPRIMES_FASTA=""
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo ""
echo "========================================================================"
echo "Configuration Summary"
echo "========================================================================"
echo "Reference FASTA: ${REFERENCE_FASTA}"
echo "Anchors FASTA: ${ANCHORS_FASTA}"
echo "Y primes FASTA: ${YPRIMES_FASTA}"
echo "X primes FASTA: ${XPRIMES_FASTA:-'Not provided (X prime detection skipped)'}"
echo "Probe FASTA: ${PROBE_FASTA:-'Not provided (verification skipped)'}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "Min percent identity: ${MIN_PIDENT}"
echo "Min alignment length: ${MIN_LENGTH}"
echo "E-value threshold: ${EVALUE}"
echo "Trim telomeric bases from boundaries: ${ADJUST_BOUNDARIES}"
if [ "${ADJUST_BOUNDARIES}" = "true" ]; then
    echo "Min feature size after trimming: ${BOUNDARY_WINDOW}bp"
fi
echo "Debug boundary adjustments: ${DEBUG_BOUNDARIES}"
echo ""

# ============================================================================
# Run Labeling Pipeline
# ============================================================================

echo ""
echo "========================================================================"
echo "Running Pre-telomeric Region Labeling"
echo "========================================================================"

# Build command with optional arguments
PROBE_ARG=""
if [ -n "${PROBE_FASTA}" ]; then
    PROBE_ARG="--probe ${PROBE_FASTA}"
fi

XPRIME_ARG=""
if [ -n "${XPRIMES_FASTA}" ]; then
    XPRIME_ARG="--xprimes ${XPRIMES_FASTA}"
fi

BOUNDARY_ARG=""
if [ "${ADJUST_BOUNDARIES}" = "true" ]; then
    BOUNDARY_ARG="--adjust-boundaries --boundary-window ${BOUNDARY_WINDOW}"
fi

DEBUG_ARG=""
if [ "${DEBUG_BOUNDARIES}" = "true" ]; then
    DEBUG_ARG="--debug-boundaries"
fi

python "${SCRIPTS_DIR}/label_pretelomeric_regions.py" \
    --reference "${REFERENCE_FASTA}" \
    --anchors "${ANCHORS_FASTA}" \
    --yprimes "${YPRIMES_FASTA}" \
    --output-dir "${OUTPUT_DIR}" \
    --prefix "${PREFIX}" \
    --threads "${THREADS}" \
    --min-pident "${MIN_PIDENT}" \
    --min-length "${MIN_LENGTH}" \
    --evalue "${EVALUE}" \
    ${XPRIME_ARG} \
    ${PROBE_ARG} \
    ${BOUNDARY_ARG} \
    ${DEBUG_ARG}

# ============================================================================
# Extract Y Prime Sequences to FASTA
# ============================================================================

echo ""
echo "========================================================================"
echo "Extracting Y Prime Sequences to FASTA"
echo "========================================================================"

LABELED_TSV="${OUTPUT_DIR}/${PREFIX}.tsv"
YPRIME_OUTPUT_FASTA="${OUTPUT_DIR}/extracted_yprimes_${STRAIN_ID}.fasta"

if [ -f "${LABELED_TSV}" ]; then
    python "${SCRIPTS_DIR}/extract_yprime_fasta.py" \
        --labeled-tsv "${LABELED_TSV}" \
        --reference "${REFERENCE_FASTA}" \
        --output "${YPRIME_OUTPUT_FASTA}" \
        --strain "${STRAIN_ID}"
    echo "Y prime sequences extracted to: ${YPRIME_OUTPUT_FASTA}"

    # Cluster Y primes by pairwise identity (automated group assignment)
    echo ""
    echo "Clustering Y prime sequences by pairwise identity..."
    CLUSTERED_YPRIME_FASTA="${OUTPUT_DIR}/extracted_yprimes_${STRAIN_ID}_clustered.fasta"
    CLUSTER_OUTPUT_DIR="${OUTPUT_DIR}/yprime_clustering_${STRAIN_ID}"

    python "${SCRIPTS_DIR}/cluster_yprimes_paper_method.py" \
        "${YPRIME_OUTPUT_FASTA}" \
        --output-fasta "${CLUSTERED_YPRIME_FASTA}" \
        --output-dir   "${CLUSTER_OUTPUT_DIR}" \
        --linkage      "${YPRIME_LINKAGE}" \
        --stop-mode    "${YPRIME_STOP_MODE}" \
        --threads      "${THREADS}"

    # Replace extracted FASTA with clustered version
    mv "${CLUSTERED_YPRIME_FASTA}" "${YPRIME_OUTPUT_FASTA}"
    echo "Y prime sequences clustered: ${YPRIME_OUTPUT_FASTA}"
    echo "Clustering results: ${CLUSTER_OUTPUT_DIR}/"
else
    echo "WARNING: Labeled TSV not found, skipping Y prime extraction"
fi

echo ""
echo "========================================================================"
echo "Pipeline completed successfully!"
echo "========================================================================"
echo "Output files:"
echo "  - GFF3 annotation: ${OUTPUT_DIR}/${PREFIX}.gff3"
echo "  - BED annotation: ${OUTPUT_DIR}/${PREFIX}.bed"
echo "  - Simplified BED: ${OUTPUT_DIR}/${PREFIX}_simp.bed"
echo "  - TSV table: ${OUTPUT_DIR}/${PREFIX}.tsv"
echo "  - Structure viz: ${OUTPUT_DIR}/${PREFIX}_structure.txt"
echo "  - Anchor BLAST: ${OUTPUT_DIR}/${PREFIX}_anchor_blast.txt"
echo "  - Y prime BLAST: ${OUTPUT_DIR}/${PREFIX}_yprime_blast.txt"
if [ -n "${XPRIMES_FASTA}" ]; then
    echo "  - X prime BLAST: ${OUTPUT_DIR}/${PREFIX}_xprime_blast.txt"
fi
echo "  - Y prime FASTA: ${YPRIME_OUTPUT_FASTA}"
echo ""
echo "Visualization tips:"
echo "  - Load the GFF3 file in a genome browser (e.g., IGV, JBrowse)"
echo "  - Use the BED file with bedtools for region extraction"
echo "  - Analyze the TSV file for quantitative analysis"
echo "  - View the structure visualization for per-chromosome-end summaries"
echo ""

# ============================================================================
# Create Reference Files for Pairing Analysis
# ============================================================================

echo ""
echo "========================================================================"
echo "Creating Reference Files for Pairing Analysis"
echo "========================================================================"

# Define paths for pairing creation
SIMPLIFIED_BED="${OUTPUT_DIR}/${PREFIX}_simp.bed"
EXTRACTED_YPRIMES="${OUTPUT_DIR}/extracted_yprimes_${STRAIN_ID}.fasta"
SPACER_SEQUENCES="${OUTPUT_DIR}/${STRAIN_ID}_50kb_chopped_up_spacer_sequences.fasta"
X_ELEMENT_SEQUENCES="${OUTPUT_DIR}/${STRAIN_ID}_whole_x_regions_sequences.fasta"

# Verify required files exist
if [ ! -f "${SIMPLIFIED_BED}" ]; then
    echo "ERROR: Simplified BED file not found: ${SIMPLIFIED_BED}"
    echo "Skipping pairing creation"
else
    if [ ! -f "${EXTRACTED_YPRIMES}" ]; then
        echo "ERROR: Extracted Y primes not found: ${EXTRACTED_YPRIMES}"
        echo "Skipping pairing creation"
    else
        # Step 1: Create chopped spacer sequences
        echo ""
        echo "Step 1: Creating chopped spacer sequences..."
        python _pipeline/scripts/make_chopped_spacer_sequences.py \
            "${STRAIN_ID}" \
            "${SIMPLIFIED_BED}" \
            "${REFERENCE_FASTA}" \
            "${OUTPUT_DIR}/" \
            --fixed-50kb

        if [ -f "${SPACER_SEQUENCES}" ]; then
            echo "Created: ${SPACER_SEQUENCES}"
        else
            echo "WARNING: Spacer sequences file not created"
        fi

        # Step 2: Create X element sequences
        echo ""
        echo "Step 2: Creating X element sequences..."
        python _pipeline/scripts/make_x_element_sequences.py \
            "${STRAIN_ID}" \
            "${SIMPLIFIED_BED}" \
            "${REFERENCE_FASTA}" \
            "${OUTPUT_DIR}/"

        if [ -f "${X_ELEMENT_SEQUENCES}" ]; then
            echo "Created: ${X_ELEMENT_SEQUENCES}"
        else
            echo "WARNING: X element sequences file not created"
        fi

        # Step 3: Create spacer pairings
        echo ""
        echo "Step 3: Creating spacer pairings for RepeatMasker..."
        if [ -f "${SPACER_SEQUENCES}" ]; then
            mkdir -p "${OUTPUT_DIR}/pairings_for_spacers"
            python _pipeline/scripts/make_databases_of_pairings_for_spacers.py \
                "${STRAIN_ID}" \
                "${EXTRACTED_YPRIMES}" \
                "${SPACER_SEQUENCES}" \
                "${OUTPUT_DIR}/pairings_for_spacers/"
            echo "Created: ${OUTPUT_DIR}/pairings_for_spacers/${STRAIN_ID}_pairings/"
        else
            echo "WARNING: Skipping spacer pairings - spacer sequences not available"
        fi

        # Step 4: Cluster X elements so chr_ends with indistinguishable
        # x-element sequences share a single cluster ID at analysis time.
        echo ""
        echo "Step 4: Clustering X elements..."
        if [ -f "${X_ELEMENT_SEQUENCES}" ]; then
            CLUSTERED_X_FASTA="${OUTPUT_DIR}/clustered_x_elements_${STRAIN_ID}.fasta"
            X_CLUSTER_WORKDIR="${OUTPUT_DIR}/x_element_clustering_${STRAIN_ID}"
            mkdir -p "${X_CLUSTER_WORKDIR}"
            python _pipeline/scripts/cluster_x_elements.py \
                "${X_ELEMENT_SEQUENCES}" \
                --output-fasta "${CLUSTERED_X_FASTA}" \
                --output-dir "${X_CLUSTER_WORKDIR}" \
                --threads "${THREADS}"
            echo "Clustered x-element library: ${CLUSTERED_X_FASTA}"
            echo "Cluster work dir:            ${X_CLUSTER_WORKDIR}"
        else
            echo "WARNING: Skipping X element clustering - X element sequences not available"
        fi

        # Extracted Y primes already live in ${EXTRACTED_YPRIMES} (the per-run
        # pretelomeric_labels/ dir), which is where the recombination step
        # picks them up. No copy to references/ — keeping the Y' file
        # alongside the run it came from makes the sample/file pairing
        # obvious.
        echo ""
        echo "Extracted Y primes available at: ${EXTRACTED_YPRIMES}"
    fi
fi

# ============================================================================
# Reference structure diagrams (auto-generated overview)
# ============================================================================

echo ""
echo "========================================================================"
echo "Generating reference structure diagrams"
echo "========================================================================"

REF_STRUCT_DIR="${OUTPUT_DIR}/reference_structure_figures"
SIMP_BED_FOR_PLOT="${OUTPUT_DIR}/${PREFIX}_simp.bed"

if [ -f "${SIMP_BED_FOR_PLOT}" ]; then
    python _pipeline/scripts/plot_reference_structure.py \
        --bed        "${SIMP_BED_FOR_PLOT}" \
        --yprime-lib "${EXTRACTED_YPRIMES}" \
        --reference  "${REFERENCE_FASTA}" \
        --output-dir "${REF_STRUCT_DIR}" \
        || echo "WARNING: reference structure plot generation failed (non-fatal)"
else
    echo "Skipping reference structure plots: ${SIMP_BED_FOR_PLOT} not found"
fi

# ============================================================================
# Final Summary
# ============================================================================

echo ""
echo "========================================================================"
echo "Pipeline completed successfully!"
echo "========================================================================"
echo ""
echo "Reference files created for strain ${STRAIN_ID}:"
echo "  - Extracted Y primes: ${EXTRACTED_YPRIMES}"
echo "  - Spacer sequences: ${SPACER_SEQUENCES}"
echo "  - X element sequences: ${X_ELEMENT_SEQUENCES}"
echo "  - Spacer pairings: ${OUTPUT_DIR}/pairings_for_spacers/${STRAIN_ID}_pairings/"
echo "  - Clustered X elements: ${OUTPUT_DIR}/clustered_x_elements_${STRAIN_ID}.fasta"
echo "  - Reference structure diagrams: ${REF_STRUCT_DIR}/"
echo ""
echo "You can now run the telomere analysis pipeline with:"
echo "  qsub telomere_analysis_complete.sh dorado_${STRAIN_ID}_dayX_PromethION_no_tag_yes_rejection"
echo ""
echo "End time: $(date)"
