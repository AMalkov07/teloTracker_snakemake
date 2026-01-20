# Pre-Telomeric Region Labeling Pipeline

This pipeline labels pre-telomeric regions (anchor, X prime, and Y prime) in a newly assembled reference genome using known sequences from a reference strain.

## Overview

After running `create_ref.sh` to create a high-quality reference genome, this pipeline:

1. Uses BLAST to find similar anchor and Y prime sequences from strain 6991 in your new strain
2. Identifies anchor regions (closest to centromere)
3. Identifies Y prime regions (middle subtelomeric region)
4. Infers X prime regions (between anchor and Y prime, or after anchor if no Y prime)
5. Outputs annotations in GFF3, BED, and TSV formats for visualization and analysis

## Prerequisites

- Completed `create_ref.sh` pipeline (generates the reference FASTA)
- BLAST+ installed (`blastn` command available)
- Python 3 with pandas
- Reference sequences:
  - `references/test_anchors.fasta` (32 anchor sequences)
  - `references/repeatmasker_6991_all_y_primes.fasta` (23 Y prime sequences)

## Quick Start

### 1. Run the pipeline

After `create_ref.sh` completes successfully:

```bash
# Submit as a job (recommended)
qsub label_regions.sh

# Or run directly
bash label_regions.sh
```

### 2. Check outputs

The pipeline creates several output files in:
```
results/${BASE_NAME}/assembly_${STRAIN_ID}/pretelomeric_labels/
```

Output files:
- `pretelomeric_regions_${STRAIN_ID}.gff3` - Genome annotation in GFF3 format
- `pretelomeric_regions_${STRAIN_ID}.bed` - BED format for bedtools/IGV
- `pretelomeric_regions_${STRAIN_ID}.tsv` - Tabular format for analysis
- `pretelomeric_regions_${STRAIN_ID}_anchor_blast.txt` - Raw BLAST results for anchors
- `pretelomeric_regions_${STRAIN_ID}_yprime_blast.txt` - Raw BLAST results for Y primes

## Configuration

Edit `label_regions.sh` to customize:

```bash
# Match these to your create_ref.sh settings
BASE_NAME="dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection"
STRAIN_ID="7575"

# BLAST parameters
MIN_PIDENT=75.0        # Minimum percent identity (75% for cross-strain)
MIN_LENGTH=100         # Minimum alignment length (100bp)
EVALUE=1e-5           # E-value threshold
```

### Parameter Guidelines

- **MIN_PIDENT**: Lower values (70-80%) allow more divergent matches between strains
- **MIN_LENGTH**: Minimum 100bp recommended to filter out spurious hits
- **EVALUE**: 1e-5 is relaxed for cross-strain comparison; use 1e-10 for same-strain

## Output Format Details

### GFF3 Format

Standard genome annotation format with attributes:
- `ID`: Unique identifier for each feature
- `Name`: Human-readable name (e.g., chr1L_anchor)
- `chr_end`: Chromosome end identifier
- `source_seq`: Original sequence name from reference
- `pident`: Percent identity from BLAST
- `bitscore`: BLAST alignment score
- `query_coverage`: Percentage of query sequence aligned

### BED Format

6-column BED format:
```
chr  start  end  name  score  strand
```

### TSV Format

Tabular data with columns:
- `chr_end`: Chromosome end (e.g., chr1L, chr1R)
- `chr`: Chromosome/contig name in reference
- `start`, `end`: Region coordinates
- `length`: Region length in bp
- `strand`: + or -
- `region_type`: anchor, x_prime, or y_prime
- `source`: Origin of annotation
- `pident`, `bitscore`, `query_coverage`: BLAST quality metrics

## Region Identification Logic

### Anchor Regions
- Direct BLAST alignment from `test_anchors.fasta`
- Typically closest to the centromere
- One per chromosome end

### Y Prime Regions
- Direct BLAST alignment from `repeatmasker_6991_all_y_primes.fasta`
- May have 0, 1, or multiple Y primes per chromosome end
- Classified as Long/Short and Solo/Tandem

### X Prime Regions
- **Inferred** (not directly aligned)
- If Y prime present: region between anchor and Y prime
- If no Y prime: ~3kb region after anchor (toward telomere)
- Minimum 100bp length required

## Visualization

### IGV (Integrative Genomics Viewer)

1. Load your reference FASTA:
   ```
   results/${BASE_NAME}/assembly_${STRAIN_ID}/assembly_${STRAIN_ID}_dorado_reference.fasta
   ```

2. Load the GFF3 file:
   ```
   results/${BASE_NAME}/assembly_${STRAIN_ID}/pretelomeric_labels/pretelomeric_regions_${STRAIN_ID}.gff3
   ```

3. Navigate to chromosome ends to see labeled regions

### Bedtools

Extract sequences for specific regions:

```bash
# Extract all anchor sequences
bedtools getfasta -fi reference.fasta -bed pretelomeric_regions_7575.bed | \
  grep -A1 "anchor" > anchor_sequences.fasta

# Get statistics on region sizes
awk '$4 ~ /anchor/ {print $3-$2}' pretelomeric_regions_7575.bed | \
  datamash mean 1 median 1 min 1 max 1
```

## Troubleshooting

### No BLAST hits found

- Check that BLAST+ is installed: `blastn -version`
- Verify reference files exist in `references/` directory
- Try relaxing parameters (lower MIN_PIDENT, higher EVALUE)

### Reference FASTA not found

- Ensure `create_ref.sh` completed successfully
- Check `BASE_NAME` and `STRAIN_ID` match between scripts
- Verify the Dorado polished reference exists

### Python errors

- Ensure pandas is installed: `pip install pandas`
- Check Python 3 is available: `python3 --version`

## Advanced Usage

### Running with different references

```bash
python scripts/label_pretelomeric_regions.py \
    --reference /path/to/your/reference.fasta \
    --anchors references/test_anchors.fasta \
    --yprimes references/repeatmasker_6991_all_y_primes.fasta \
    --output-dir /path/to/output \
    --prefix my_labels \
    --threads 8 \
    --min-pident 70.0 \
    --min-length 50 \
    --evalue 1e-3
```

### Filtering results by quality

```python
import pandas as pd

# Load TSV results
df = pd.read_csv('pretelomeric_regions_7575.tsv', sep='\t')

# Filter high-confidence anchors (>90% identity, >80% coverage)
high_conf = df[
    (df['region_type'] == 'anchor') &
    (df['pident'] > 90) &
    (df['query_coverage'] > 80)
]

print(f"High-confidence anchors: {len(high_conf)}")
```

## Integration with Existing Pipeline

This pipeline integrates with your existing TeloTracker workflow:

```
Snakefile (y_prime_analysis)
    ↓
create_ref.sh (reference assembly and polishing)
    ↓
label_regions.sh (THIS PIPELINE - region labeling)
    ↓
[Your downstream analysis]
```

The labeled regions can be used for:
- Identifying recombination breakpoints
- Analyzing subtelomeric structure variation
- Guiding repeat element classification
- Comparing telomere organization across strains

## Citation

If you use this pipeline, please cite:
- BLAST: Altschul et al. (1990)
- Your TeloTracker pipeline publication (when available)

## Contact

For issues or questions about this pipeline, please open an issue in the repository.
