#!/bin/bash
# Test script for Y prime detection with fragment merging

set -e

echo "================================================================================"
echo "TESTING: Y Prime Detection with Fragment Merging"
echo "================================================================================"
echo ""

# Configuration
BASE_NAME="dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection"
STRAIN_ID="7575"

# Input files
REFERENCE="results/${BASE_NAME}/assembly_${STRAIN_ID}/assembly_${STRAIN_ID}_dorado_reference.fasta"
ANCHORS="references/test_anchors.fasta"
YPRIMES="references/repeatmasker_6991_all_y_primes.fasta"

# Output directory
OUTPUT_DIR="results/${BASE_NAME}/assembly_${STRAIN_ID}/test_yprime_detection"
PREFIX="test_yprimes_${STRAIN_ID}"

# BLAST parameters
THREADS=8
MIN_PIDENT=75.0
MIN_LENGTH=100
EVALUE=1e-5

echo "Configuration:"
echo "  Reference: $REFERENCE"
echo "  Anchors: $ANCHORS"
echo "  Y primes: $YPRIMES"
echo "  Output: $OUTPUT_DIR"
echo "  Threads: $THREADS"
echo ""

# Verify files exist
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference file not found: $REFERENCE"
    exit 1
fi

if [ ! -f "$ANCHORS" ]; then
    echo "ERROR: Anchors file not found: $ANCHORS"
    exit 1
fi

if [ ! -f "$YPRIMES" ]; then
    echo "ERROR: Y primes file not found: $YPRIMES"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Running Y prime detection with fragment merging..."
echo ""

# Run the script
python scripts/label_pretelomeric_regions.py \
    --reference "$REFERENCE" \
    --anchors "$ANCHORS" \
    --yprimes "$YPRIMES" \
    --output-dir "$OUTPUT_DIR" \
    --prefix "$PREFIX" \
    --threads "$THREADS" \
    --min-pident "$MIN_PIDENT" \
    --min-length "$MIN_LENGTH" \
    --evalue "$EVALUE"

echo ""
echo "================================================================================"
echo "Test completed! Check results in:"
echo "  $OUTPUT_DIR"
echo "================================================================================"
echo ""
echo "Output files:"
ls -lh "$OUTPUT_DIR"
echo ""
echo "To view Y prime results:"
echo "  grep y_prime $OUTPUT_DIR/${PREFIX}.tsv"
echo ""
echo "To count Y primes per chr end:"
echo "  awk '\$7==\"y_prime\" {print \$1}' $OUTPUT_DIR/${PREFIX}.tsv | sort | uniq -c"
