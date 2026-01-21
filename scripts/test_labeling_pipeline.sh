#!/bin/bash
# Test script to validate the labeling pipeline setup

set -e

echo "================================================================================"
echo "Pre-Telomeric Region Labeling Pipeline - Validation Test"
echo "================================================================================"
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

ERRORS=0
WARNINGS=0

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check if file exists
check_file() {
    if [ -f "$1" ]; then
        echo -e "${GREEN}✓${NC} Found: $1"
        return 0
    else
        echo -e "${RED}✗${NC} Missing: $1"
        ((ERRORS++))
        return 1
    fi
}

# Function to check if directory exists
check_dir() {
    if [ -d "$1" ]; then
        echo -e "${GREEN}✓${NC} Found: $1"
        return 0
    else
        echo -e "${YELLOW}!${NC} Directory not found: $1"
        ((WARNINGS++))
        return 1
    fi
}

echo "1. Checking Required Commands"
echo "------------------------------"

if command_exists blastn; then
    BLAST_VERSION=$(blastn -version | head -1)
    echo -e "${GREEN}✓${NC} blastn is installed: $BLAST_VERSION"
else
    echo -e "${RED}✗${NC} blastn not found (required for BLAST)"
    echo "  Install with: conda install -c bioconda blast"
    ((ERRORS++))
fi

if command_exists python3; then
    PYTHON_VERSION=$(python3 --version)
    echo -e "${GREEN}✓${NC} python3 is installed: $PYTHON_VERSION"
else
    echo -e "${RED}✗${NC} python3 not found"
    ((ERRORS++))
fi

if python3 -c "import pandas" 2>/dev/null; then
    PANDAS_VERSION=$(python3 -c "import pandas; print(pandas.__version__)")
    echo -e "${GREEN}✓${NC} pandas is installed: version $PANDAS_VERSION"
else
    echo -e "${RED}✗${NC} pandas not found"
    echo "  Install with: pip install pandas"
    ((ERRORS++))
fi

echo ""
echo "2. Checking Pipeline Scripts"
echo "-----------------------------"

check_file "label_regions.sh"
check_file "scripts/label_pretelomeric_regions.py"
check_file "scripts/visualize_labeled_regions.py"
check_file "scripts/extract_labeled_sequences.py"

# Check if scripts are executable
for script in label_regions.sh scripts/label_pretelomeric_regions.py \
              scripts/visualize_labeled_regions.py scripts/extract_labeled_sequences.py; do
    if [ -x "$script" ]; then
        echo -e "${GREEN}✓${NC} Executable: $script"
    else
        echo -e "${YELLOW}!${NC} Not executable: $script (run: chmod +x $script)"
        ((WARNINGS++))
    fi
done

echo ""
echo "3. Checking Reference Files"
echo "----------------------------"

check_file "references/test_anchors.fasta"
check_file "references/repeatmasker_6991_all_y_primes.fasta"

# Count sequences in reference files
if [ -f "references/test_anchors.fasta" ]; then
    ANCHOR_COUNT=$(grep -c "^>" references/test_anchors.fasta)
    echo "  → $ANCHOR_COUNT anchor sequences found"
fi

if [ -f "references/repeatmasker_6991_all_y_primes.fasta" ]; then
    YPRIME_COUNT=$(grep -c "^>" references/repeatmasker_6991_all_y_primes.fasta)
    echo "  → $YPRIME_COUNT Y prime sequences found"
fi

echo ""
echo "4. Checking Documentation"
echo "-------------------------"

check_file "LABELING_PIPELINE_README.md"
check_file "WORKFLOW_EXAMPLE.md"
check_file "QUICK_START_LABELING.md"
check_file "LABELING_PIPELINE_SUMMARY.txt"

echo ""
echo "5. Checking create_ref.sh Output"
echo "---------------------------------"

# Try to find any assembly directories
ASSEMBLY_DIRS=$(find results -name "assembly_*" -type d 2>/dev/null | head -5)

if [ -z "$ASSEMBLY_DIRS" ]; then
    echo -e "${YELLOW}!${NC} No assembly directories found in results/"
    echo "  This is expected if you haven't run create_ref.sh yet"
    echo "  Run create_ref.sh first, then run label_regions.sh"
    ((WARNINGS++))
else
    echo "Found assembly directories:"
    echo "$ASSEMBLY_DIRS" | while read dir; do
        echo "  - $dir"

        # Check for dorado reference
        DORADO_REF=$(find "$dir" -name "*_dorado_reference.fasta" 2>/dev/null)
        if [ -n "$DORADO_REF" ]; then
            echo -e "    ${GREEN}✓${NC} Has dorado reference: $(basename $DORADO_REF)"
        else
            echo -e "    ${YELLOW}!${NC} No dorado reference found"
        fi

        # Check for existing labels
        LABELS_DIR="$dir/pretelomeric_labels"
        if [ -d "$LABELS_DIR" ]; then
            echo -e "    ${GREEN}✓${NC} Labels already generated: $LABELS_DIR"

            # Count label files
            if [ -f "$LABELS_DIR"/*.tsv ]; then
                TSV_FILE=$(ls "$LABELS_DIR"/*.tsv 2>/dev/null | head -1)
                if [ -n "$TSV_FILE" ]; then
                    LABEL_COUNT=$(tail -n +2 "$TSV_FILE" | wc -l)
                    echo "      → $LABEL_COUNT labeled regions found"
                fi
            fi
        fi
    done
fi

echo ""
echo "6. Configuration Check"
echo "----------------------"

# Extract configuration from label_regions.sh
if [ -f "label_regions.sh" ]; then
    BASE_NAME=$(grep "^BASE_NAME=" label_regions.sh | cut -d'"' -f2)
    STRAIN_ID=$(grep "^STRAIN_ID=" label_regions.sh | cut -d'"' -f2)

    echo "Current configuration in label_regions.sh:"
    echo "  BASE_NAME: $BASE_NAME"
    echo "  STRAIN_ID: $STRAIN_ID"
    echo ""

    # Check if this matches an existing directory
    EXPECTED_DIR="results/${BASE_NAME}/assembly_${STRAIN_ID}"
    if [ -d "$EXPECTED_DIR" ]; then
        echo -e "${GREEN}✓${NC} Configuration matches existing directory: $EXPECTED_DIR"

        # Check for reference file
        EXPECTED_REF="${EXPECTED_DIR}/assembly_${STRAIN_ID}_dorado_reference.fasta"
        if [ -f "$EXPECTED_REF" ]; then
            echo -e "${GREEN}✓${NC} Reference file exists: $(basename $EXPECTED_REF)"

            # Get reference stats
            REF_SIZE=$(du -h "$EXPECTED_REF" | cut -f1)
            REF_SEQS=$(grep -c "^>" "$EXPECTED_REF")
            echo "    Size: $REF_SIZE, Sequences: $REF_SEQS"
            echo ""
            echo -e "${GREEN}✓✓✓ Ready to run label_regions.sh! ✓✓✓${NC}"
        else
            echo -e "${YELLOW}!${NC} Reference file not found: $(basename $EXPECTED_REF)"
            echo "  Run create_ref.sh first to generate the reference"
        fi
    else
        echo -e "${YELLOW}!${NC} Expected directory not found: $EXPECTED_DIR"
        echo "  Either:"
        echo "    1. Run create_ref.sh first, or"
        echo "    2. Update BASE_NAME and STRAIN_ID in label_regions.sh"
    fi
fi

echo ""
echo "================================================================================"
echo "Validation Summary"
echo "================================================================================"

if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo -e "${GREEN}✓ All checks passed!${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. Run: bash label_regions.sh"
    echo "  2. Check outputs in: results/${BASE_NAME}/assembly_${STRAIN_ID}/pretelomeric_labels/"
    echo "  3. Visualize with: python scripts/visualize_labeled_regions.py --input-tsv <tsv_file>"
elif [ $ERRORS -eq 0 ]; then
    echo -e "${YELLOW}⚠ Validation completed with $WARNINGS warning(s)${NC}"
    echo ""
    echo "Warnings are not critical, but you should review them."
else
    echo -e "${RED}✗ Validation failed with $ERRORS error(s) and $WARNINGS warning(s)${NC}"
    echo ""
    echo "Please fix the errors above before running the pipeline."
    exit 1
fi

echo ""
echo "For more information, see:"
echo "  - QUICK_START_LABELING.md (quick reference)"
echo "  - LABELING_PIPELINE_README.md (full documentation)"
echo "  - WORKFLOW_EXAMPLE.md (examples and workflows)"
echo ""
