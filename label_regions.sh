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
#REFERENCE_DIR="results/${BASE_NAME}/assembly_${STRAIN_ID}"
REFERENCE_DIR="labeling_test"
#REFERENCE_FASTA="${REFERENCE_DIR}/assembly_${STRAIN_ID}_dorado_reference.fasta"
REFERENCE_FASTA="${REFERENCE_DIR}/${STRAIN_ID}.fasta"

# Reference sequences for labeling (from strain 6991)
REFERENCES_DIR="/Shared/malkova_lab/Ivan/nanopore_sequencing/reference_files"
ANCHORS_FASTA="${REFERENCES_DIR}/test_anchors.fasta"
YPRIMES_FASTA="${REFERENCES_DIR}/repeatmasker_6991_all_y_primes.fasta"
XPRIMES_FASTA="${REFERENCES_DIR}/6991_xprimes.fasta"  # X prime sequences for detection
PROBE_FASTA="${REFERENCES_DIR}/probe.fasta"  # Y prime probe for verification

# Scripts directory
SCRIPTS_DIR="/Shared/malkova_lab/Ivan/nanopore_sequencing/reference_files/scripts"

# Output configuration
OUTPUT_DIR="results/${BASE_NAME}/pretelomeric_labels"
PREFIX="pretelomeric_regions_${STRAIN_ID}"

# Thread configuration
THREADS=56

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
SPACER_SEQUENCES="references/${STRAIN_ID}_50kb_chopped_up_spacer_sequences.fasta"
X_ELEMENT_SEQUENCES="references/${STRAIN_ID}_whole_x_regions_sequences.fasta"

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
        python scripts/make_chopped_spacer_sequences.py \
            "${STRAIN_ID}" \
            "${SIMPLIFIED_BED}" \
            "${REFERENCE_FASTA}" \
            "references/" \
            --fixed-50kb

        if [ -f "${SPACER_SEQUENCES}" ]; then
            echo "Created: ${SPACER_SEQUENCES}"
        else
            echo "WARNING: Spacer sequences file not created"
        fi

        # Step 2: Create X element sequences
        echo ""
        echo "Step 2: Creating X element sequences..."
        python scripts/make_x_element_sequences.py \
            "${STRAIN_ID}" \
            "${SIMPLIFIED_BED}" \
            "${REFERENCE_FASTA}" \
            "references/"

        if [ -f "${X_ELEMENT_SEQUENCES}" ]; then
            echo "Created: ${X_ELEMENT_SEQUENCES}"
        else
            echo "WARNING: X element sequences file not created"
        fi

        # Step 3: Create spacer pairings
        echo ""
        echo "Step 3: Creating spacer pairings for RepeatMasker..."
        if [ -f "${SPACER_SEQUENCES}" ]; then
            python scripts/make_databases_of_pairings_for_spacers.py \
                "${STRAIN_ID}" \
                "${EXTRACTED_YPRIMES}" \
                "${SPACER_SEQUENCES}" \
                "references/pairings_for_spacers/"
            echo "Created: references/pairings_for_spacers/${STRAIN_ID}_pairings/"
        else
            echo "WARNING: Skipping spacer pairings - spacer sequences not available"
        fi

        # Step 4: Create X element pairings
        echo ""
        echo "Step 4: Creating X element pairings for RepeatMasker..."
        if [ -f "${X_ELEMENT_SEQUENCES}" ]; then
            python scripts/make_databases_of_pairings_for_x_elements.py \
                "${STRAIN_ID}" \
                "${EXTRACTED_YPRIMES}" \
                "${X_ELEMENT_SEQUENCES}" \
                "references/pairings_for_x_element_ends/"
            echo "Created: references/pairings_for_x_element_ends/${STRAIN_ID}_pairings/"
        else
            echo "WARNING: Skipping X element pairings - X element sequences not available"
        fi

        # Copy extracted Y primes to references directory for easy access
        echo ""
        echo "Copying extracted Y primes to references directory..."
        cp "${EXTRACTED_YPRIMES}" "references/extracted_yprimes_${STRAIN_ID}.fasta"
        echo "Copied: references/extracted_yprimes_${STRAIN_ID}.fasta"
    fi
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
echo "  - Extracted Y primes: references/extracted_yprimes_${STRAIN_ID}.fasta"
echo "  - Spacer sequences: references/${STRAIN_ID}_50kb_chopped_up_spacer_sequences.fasta"
echo "  - X element sequences: references/${STRAIN_ID}_whole_x_regions_sequences.fasta"
echo "  - Spacer pairings: references/pairings_for_spacers/${STRAIN_ID}_pairings/"
echo "  - X element pairings: references/pairings_for_x_element_ends/${STRAIN_ID}_pairings/"
echo ""
echo "You can now run the telomere analysis pipeline with:"
echo "  qsub telomere_analysis_complete.sh dorado_${STRAIN_ID}_dayX_PromethION_no_tag_yes_rejection"
echo ""
echo "End time: $(date)"
