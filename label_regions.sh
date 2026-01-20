#!/bin/bash
#$ -q UI-GPU
#$ -pe smp 56
#$ -l gpu=true
#$ -j y
#$ -cwd

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# ============================================================================
# Configuration - Edit these values to match your create_ref.sh setup
# ============================================================================

# Must match the values from create_ref.sh
BASE_NAME="dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection"
STRAIN_ID="7575"

# Input: Reference created by create_ref.sh (Dorado polished reference)
REFERENCE_DIR="results/${BASE_NAME}/assembly_${STRAIN_ID}"
REFERENCE_FASTA="${REFERENCE_DIR}/assembly_${STRAIN_ID}_dorado_reference.fasta"

# Reference sequences for labeling (from strain 6991)
ANCHORS_FASTA="references/test_anchors.fasta"
YPRIMES_FASTA="references/repeatmasker_6991_all_y_primes.fasta"
XPRIMES_FASTA="references/6991_xprimes.fasta"  # X prime sequences for detection
PROBE_FASTA="references/probe.fasta"  # Y prime probe for verification

# Output configuration
OUTPUT_DIR="${REFERENCE_DIR}/pretelomeric_labels"
PREFIX="pretelomeric_regions_${STRAIN_ID}"

# Thread configuration
THREADS=56

# BLAST parameters
MIN_PIDENT=75.0        # Minimum percent identity (lowered for cross-strain comparison)
MIN_LENGTH=100         # Minimum alignment length
EVALUE=1e-5           # E-value threshold (relaxed for cross-strain)

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

python scripts/label_pretelomeric_regions.py \
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
    ${PROBE_ARG}

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
    python scripts/extract_yprime_fasta.py \
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
echo "End time: $(date)"
