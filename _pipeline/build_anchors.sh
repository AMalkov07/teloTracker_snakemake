#$ -q UI
#$ -pe smp 56
#$ -j y
#$ -cwd

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# ============================================================================
# Configuration - Edit these values for your analysis
# ============================================================================

# Required inputs
BASE_NAME="6212_TLC1_and_tlc1_delta_day0_with_selection"
STRAIN_ID="6212"
ANCHOR_SET="6212_anchors"

# Thread configuration
THREADS=56

# Assembly configuration
GENOME_SIZE="12m"
TARGET_BASES="1200000000"  # ~100x of 12 Mb (Wick 2026: excess depth degrades assembly)
STOP_AFTER_ASSEMBLY="true" # assess assembly QC before deriving anchors from it

###### Below does NOT need to be adjusted ######

# Paths relative to the repo root
READS_FASTQ="results/${BASE_NAME}/_pipeline/${BASE_NAME}.fastq"
OUTPUT_DIR="results/${BASE_NAME}/_pipeline/anchor_build_${STRAIN_ID}"
ADAPTER_FILE="_pipeline/references/offical_nanopore_adapter_seq+trunc.txt"
S288C_REF="_pipeline/references/s288c/s288c_genome.fasta"
MITO_CONTIG="NC_001224.1"
MITO_REF="_pipeline/references/s288c/s288c_mito.fasta"
SCRIPTS_DIR="_pipeline/scripts"

# ============================================================================
# Setup
# ============================================================================

echo "Building a strain-specific anchor set"
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "Strain: ${STRAIN_ID}   Sample: ${BASE_NAME}"

# run_pipeline.py prepends an equivalent prologue; guard so the script also works
# when qsub'ed directly.
if ! declare -F conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate consensus
export LANG=C LC_ALL=C

mkdir -p "${OUTPUT_DIR}/qc"

DENOVO_TRIMMED="${OUTPUT_DIR}/${STRAIN_ID}_denovo_trimmed.fastq"
DENOVO_NOMITO="${OUTPUT_DIR}/${STRAIN_ID}_denovo_nomito.fastq"
DENOVO_SUBSAMPLED="${OUTPUT_DIR}/${STRAIN_ID}_denovo_subsampled.fastq"
FLYE_ASM_DIR="${OUTPUT_DIR}/flye_denovo"
ASSEMBLY="${FLYE_ASM_DIR}/assembly.fasta"

# ============================================================================
# Step 0: Produce filtered reads (anchor-independent)
# ============================================================================
# NOTE: create_ref.sh cannot be reused here. Its Step 0 runs
# `snakemake through_y_prime_analysis`, which depends on rule blast_anchors and
# therefore on the very anchors we are trying to build. Target the filtered-FASTQ
# rule directly instead -- prep_raw_fastq -> filter_reads touches no anchor.

echo ""
echo "========================================================================"
echo "Step 0: Generating filtered reads (no anchors required)"
echo "========================================================================"

snakemake -s _pipeline/Snakefile "${READS_FASTQ}" -c "${THREADS}"

if [ ! -f "${READS_FASTQ}" ]; then
    echo "ERROR: Required file not found: ${READS_FASTQ}"
    exit 1
fi
echo "Filtered reads: ${READS_FASTQ}"

# ============================================================================
# Step 1: Trim adapters
# ============================================================================
# Deliberately NOT written to assembly_<strain>/<prefix>_trimmed.fastq -- that path
# belongs to create_ref.sh and pre-creating it invites the stale-file failures already
# documented at create_ref.sh:207-212.

echo ""
echo "========================================================================"
echo "Step 1: Trimming adapters with porechop_abi"
echo "========================================================================"

if [ -s "${DENOVO_TRIMMED}" ]; then
    echo "Already present, skipping: ${DENOVO_TRIMMED}"
else
    porechop_abi -i "${READS_FASTQ}" -o "${DENOVO_TRIMMED}" \
        -cap "${ADAPTER_FILE}" -t "${THREADS}" --no_split -ddb
fi
echo "Adapter-trimmed FASTQ: ${DENOVO_TRIMMED}"

# ============================================================================
# Step 1.5: Deplete mitochondrial reads
# ============================================================================
# THIS IS THE STEP THAT MAKES THE ASSEMBLY POSSIBLE AT ALL.
#
# Measured on 6212 (30,000 random raw reads mapped to S288C):
#   mitochondrial  97,413,984 bp  (81.9% of all aligned bases)
#   nuclear        21,556,817 bp  (18.1%)
# Among reads >=2 kb the nuclear share is only 16.5%. Flye was therefore handed
# 1.2 Gb that was ~5/6 mitochondrial: ~16x real nuclear depth, while the mito
# contigs sat at 9,000-11,000x. That skews Flye's global coverage model -- genuine
# single-copy nuclear sequence looks like noise next to a 10,000x organelle -- and
# the assembly collapsed to 546 kb of the 12.07 Mb genome.
#
# Depleting mito leaves ~433 Mb / ~36x honest nuclear coverage.
#
# Threshold note: yeast has NUMTs (mito fragments inserted in the nuclear genome),
# so a read is only discarded when a mito alignment covers >=50% of its length.
# Discarding on any mito hit would throw away real subtelomeric reads.

echo ""
echo "========================================================================"
echo "Step 1.5: Depleting mitochondrial reads"
echo "========================================================================"

if [ -s "${DENOVO_NOMITO}" ]; then
    echo "Already present, skipping: ${DENOVO_NOMITO}"
else
    if [ ! -s "${MITO_REF}" ]; then
        samtools faidx "${S288C_REF}" "${MITO_CONTIG}" > "${MITO_REF}"
    fi

    MITO_IDS="${OUTPUT_DIR}/qc/mito_read_ids.txt"
    minimap2 -x map-ont -t "${THREADS}" --secondary=no "${MITO_REF}" "${DENOVO_TRIMMED}" 2>/dev/null \
        | awk '{ cov[$1] += $4 - $3; qlen[$1] = $2 }
               END { for (r in cov) if (cov[r] >= 0.5 * qlen[r]) print r }' \
        > "${MITO_IDS}"

    echo "Reads identified as mitochondrial: $(wc -l < "${MITO_IDS}")"
    seqkit grep -v -f "${MITO_IDS}" "${DENOVO_TRIMMED}" > "${DENOVO_NOMITO}"
fi
seqkit stats -T "${DENOVO_NOMITO}" | tee "${OUTPUT_DIR}/qc/nomito_stats.tsv"

# ============================================================================
# Step 2: Subsample UNIFORMLY AT RANDOM to the target depth
# ============================================================================
# With mito depleted this is usually a no-op (available < target), which is the
# intended behaviour -- it exists only to bound memory on a genuinely deep library.
#
# Do NOT bound memory by taking the longest reads instead. That was tried here and
# the result is misleading: a longest-first cut to 720 Mb (threshold 28,989 bp,
# 15,714 / 259,353 reads) gave a LARGER assembly (3.3 Mb) than random sampling
# (546 kb), which looks like an argument for length selection but is not. Mito
# content is strongly length-dependent in this library -- 96% of 2-10 kb reads but
# only 68% of >20 kb reads -- so the length cut was acting as a crude, accidental
# mito filter. Step 1.5 does that job directly and correctly. Wick 2026
# (rrwick.github.io/2026/02/05/read_qc_testing.html) advises against length-selecting
# QC for the usual reasons; the accidental benefit seen here does not overturn that.
#
# Length/quality filtering already happened upstream in filter_reads.py (>=2000 bp,
# qs>=10), which is the equivalent of Wick's `chopper -q 10 -l 1000` first pass.

echo ""
echo "========================================================================"
echo "Step 2: Random subsample to ~${TARGET_BASES} bases"
echo "========================================================================"

if [ -s "${DENOVO_SUBSAMPLED}" ]; then
    echo "Already present, skipping: ${DENOVO_SUBSAMPLED}"
else
    LEN_TABLE="${OUTPUT_DIR}/qc/read_lengths.tsv"
    if [ ! -s "${LEN_TABLE}" ]; then
        seqkit fx2tab -nil "${DENOVO_NOMITO}" | cut -f2 > "${LEN_TABLE}"
    fi

    TOTAL_BASES=$(awk '{s+=$1} END{print s}' "${LEN_TABLE}")
    PROPORTION=$(awk -v t="${TARGET_BASES}" -v s="${TOTAL_BASES}" \
                     'BEGIN{p = t/s; if (p > 1) p = 1; printf "%.4f", p}')
    echo "Available: ${TOTAL_BASES} bp   target: ${TARGET_BASES} bp   proportion: ${PROPORTION}"

    if [ "${PROPORTION}" == "1.0000" ]; then
        echo "Available is at or below target; using all reads."
        ln -sf "$(basename "${DENOVO_NOMITO}")" "${DENOVO_SUBSAMPLED}"
    else
        seqkit sample -p "${PROPORTION}" -s 42 "${DENOVO_NOMITO}" > "${DENOVO_SUBSAMPLED}"
    fi
fi
seqkit stats -T "${DENOVO_SUBSAMPLED}" | tee "${OUTPUT_DIR}/qc/subsampled_stats.tsv"

# ============================================================================
# Step 3: De novo assembly with Flye
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 3: De novo assembly with Flye"
echo "========================================================================"

if [ -s "${ASSEMBLY}" ]; then
    echo "Already present, skipping assembly: ${ASSEMBLY}"
else
    flye --nano-hq "${DENOVO_SUBSAMPLED}" \
         --out-dir "${FLYE_ASM_DIR}" \
         --genome-size "${GENOME_SIZE}" \
         --asm-coverage 50 \
         --threads "${THREADS}" \
         --no-alt-contigs
fi

if [ ! -s "${ASSEMBLY}" ]; then
    echo "ERROR: Flye produced no assembly at ${ASSEMBLY}"
    exit 1
fi
echo "Assembly: ${ASSEMBLY}"

# ============================================================================
# Step 4: Assembly QC
# ============================================================================

echo ""
echo "========================================================================"
echo "Step 4: Assembly QC"
echo "========================================================================"

python3 "${SCRIPTS_DIR}/assembly_qc.py" \
    --assembly "${ASSEMBLY}" \
    --assembly-info "${FLYE_ASM_DIR}/assembly_info.txt" \
    --s288c "${S288C_REF}" \
    --threads "${THREADS}" \
    --out-prefix "${OUTPUT_DIR}/qc/${STRAIN_ID}"

if [ "${STOP_AFTER_ASSEMBLY}" == "true" ]; then
    echo ""
    echo "========================================================================"
    echo "STOPPING after assembly, as configured."
    echo ""
    echo "Review ${OUTPUT_DIR}/qc/${STRAIN_ID}_assembly_qc.txt before continuing."
    echo "The telomere-completeness table is the one that matters: ends without a"
    echo "telomeric tract will place their anchor further inward than the 6991 set."
    echo "Set STOP_AFTER_ASSEMBLY=false to proceed to anchor derivation."
    echo "========================================================================"
    exit 0
fi

echo ""
echo "Remaining stages (contig assignment, anchor derivation) not yet enabled."
echo "DONE: $(date)"
