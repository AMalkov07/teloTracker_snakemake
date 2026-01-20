# Complete Workflow Example: From Reads to Labeled Reference

This document provides a step-by-step example of the complete TeloTracker pipeline, from raw reads to a fully labeled reference genome with annotated pre-telomeric regions.

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────┐
│ Step 1: Snakemake Pipeline                                  │
│ - Basecalling (if needed)                                   │
│ - Adapter detection                                         │
│ - Y prime analysis                                          │
│ Output: results/${BASE_NAME}/${BASE_NAME}_post_y_prime_probe.tsv │
│         results/${BASE_NAME}/${BASE_NAME}.fastq             │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 2: Reference Creation (create_ref.sh)                  │
│ - Select 75th percentile reads                              │
│ - Adapter trimming with porechop_abi                        │
│ - Initial alignment to reference                            │
│ - Reference extension using soft-clipped bases              │
│ - Flye polishing                                            │
│ - Dorado alignment and polishing                            │
│ Output: assembly_${STRAIN_ID}_dorado_reference.fasta        │
│         assembly_${STRAIN_ID}_final_alignment_to_dorado_ref.bam │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 3: Region Labeling (label_regions.sh) ← NEW!           │
│ - BLAST anchors from reference strain                       │
│ - BLAST Y primes from reference strain                      │
│ - Infer X prime regions                                     │
│ Output: pretelomeric_regions_${STRAIN_ID}.gff3/.bed/.tsv    │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 4: Analysis and Visualization                          │
│ - Visualize labels in genome browser                        │
│ - Extract labeled sequences                                 │
│ - Downstream analysis                                       │
└─────────────────────────────────────────────────────────────┘
```

## Complete Example

### Step 1: Run Snakemake Pipeline

This step is already integrated into `create_ref.sh`, so you typically don't run it separately:

```bash
# This is done automatically in create_ref.sh at Step 0
snakemake -s Snakefile --until y_prime_analysis -c 56
```

### Step 2: Create Reference Genome

Edit `create_ref.sh` configuration:

```bash
# Required inputs
BASE_NAME="dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection"
STRAIN_ID="7575"
```

Run the pipeline:

```bash
# Submit as job
qsub create_ref.sh

# Or run directly
bash create_ref.sh
```

Expected outputs in `results/${BASE_NAME}/assembly_${STRAIN_ID}/`:
- `assembly_7575_dorado_reference.fasta` - Final polished reference
- `assembly_7575_final_alignment_to_dorado_ref.bam` - Alignment of all reads to reference
- Various intermediate files from each step

### Step 3: Label Pre-Telomeric Regions

After `create_ref.sh` completes successfully:

```bash
# Submit as job
qsub label_regions.sh

# Or run directly
bash label_regions.sh
```

Expected outputs in `results/${BASE_NAME}/assembly_${STRAIN_ID}/pretelomeric_labels/`:
- `pretelomeric_regions_7575.gff3` - Genome annotation
- `pretelomeric_regions_7575.bed` - BED format
- `pretelomeric_regions_7575.tsv` - Tabular data
- `pretelomeric_regions_7575_anchor_blast.txt` - BLAST results for anchors
- `pretelomeric_regions_7575_yprime_blast.txt` - BLAST results for Y primes

### Step 4: Visualize and Analyze Results

#### View Summary Statistics

```bash
python scripts/visualize_labeled_regions.py \
    --input-tsv results/dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection/assembly_7575/pretelomeric_labels/pretelomeric_regions_7575.tsv
```

Example output:
```
================================================================================
SUMMARY STATISTICS
================================================================================

Region counts by type:
  anchor         :  32
  x_prime        :  32
  y_prime        :  18

Total chromosome ends with labels: 32

Region size statistics (bp):
  anchor:
    Count:   32
    Mean:    1847
    Median:  1823
    Min:     1245
    Max:     2567
  ...
```

#### Export Summary Table

```bash
python scripts/visualize_labeled_regions.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --export-summary summary_by_chr_end.tsv
```

#### Extract Sequences

Extract all labeled regions:

```bash
python scripts/extract_labeled_sequences.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --reference results/.../assembly_7575_dorado_reference.fasta \
    --output-dir extracted_sequences \
    --prefix strain_7575
```

Extract only anchors:

```bash
python scripts/extract_labeled_sequences.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --reference results/.../assembly_7575_dorado_reference.fasta \
    --output-dir extracted_sequences \
    --prefix strain_7575 \
    --region-types anchor
```

Extract specific chromosome ends:

```bash
python scripts/extract_labeled_sequences.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --reference results/.../assembly_7575_dorado_reference.fasta \
    --output-dir extracted_sequences \
    --prefix strain_7575_chr1 \
    --chr-ends chr1L chr1R
```

#### Visualize in IGV

1. Open IGV
2. Load reference: `Genomes > Load Genome from File`
   - Select: `assembly_7575_dorado_reference.fasta`
3. Load annotations: `File > Load from File`
   - Select: `pretelomeric_regions_7575.gff3`
4. Navigate to chromosome ends to see labeled regions

#### Use bedtools for Analysis

Get region sizes:

```bash
awk '{print $3-$2, $4}' pretelomeric_regions_7575.bed | \
    sort -k2 | \
    datamash -g 2 mean 1 median 1 min 1 max 1 count 1
```

Extract sequences using bedtools:

```bash
bedtools getfasta \
    -fi assembly_7575_dorado_reference.fasta \
    -bed pretelomeric_regions_7575.bed \
    -name -s > all_labeled_regions.fasta
```

Find overlaps with other features:

```bash
# If you have other annotations (e.g., genes.bed)
bedtools intersect \
    -a pretelomeric_regions_7575.bed \
    -b genes.bed \
    -wa -wb > overlapping_features.txt
```

## Downstream Analysis Examples

### Compare Region Organization Across Strains

If you've run the pipeline for multiple strains:

```bash
# Strain 6991
python scripts/visualize_labeled_regions.py \
    --input-tsv results/.../pretelomeric_regions_6991.tsv \
    --export-summary summary_6991.tsv

# Strain 7575
python scripts/visualize_labeled_regions.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --export-summary summary_7575.tsv

# Compare in R or Python
```

### Analyze Y Prime Presence/Absence

```python
import pandas as pd

df = pd.read_csv('pretelomeric_regions_7575.tsv', sep='\t')

# Get all chromosome ends
chr_ends_all = set(df['chr_end'])

# Get chr ends with Y primes
chr_ends_with_yprime = set(df[df['region_type'] == 'y_prime']['chr_end'])

# Chr ends without Y primes
chr_ends_no_yprime = chr_ends_all - chr_ends_with_yprime

print(f"Chromosome ends with Y primes: {len(chr_ends_with_yprime)}")
print(f"Chromosome ends without Y primes: {len(chr_ends_no_yprime)}")
print(f"Y prime presence: {sorted(chr_ends_with_yprime)}")
print(f"Y prime absence: {sorted(chr_ends_no_yprime)}")
```

### Measure X Prime Lengths

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('pretelomeric_regions_7575.tsv', sep='\t')
xprimes = df[df['region_type'] == 'x_prime']

# Plot distribution
plt.figure(figsize=(10, 6))
plt.hist(xprimes['length'], bins=30, edgecolor='black')
plt.xlabel('X Prime Length (bp)')
plt.ylabel('Frequency')
plt.title('Distribution of X Prime Lengths in Strain 7575')
plt.savefig('xprime_length_distribution.png', dpi=300)
plt.close()

# Summary statistics
print(f"X Prime Length Statistics:")
print(f"  Mean: {xprimes['length'].mean():.0f} bp")
print(f"  Median: {xprimes['length'].median():.0f} bp")
print(f"  Range: {xprimes['length'].min():.0f} - {xprimes['length'].max():.0f} bp")
```

### Identify Recombination Breakpoints

Combine with your existing Y prime analysis:

```python
import pandas as pd

# Load labeled regions
labels = pd.read_csv('pretelomeric_regions_7575.tsv', sep='\t')

# Load Y prime probe results (from Snakemake pipeline)
yprimes = pd.read_csv('results/.../post_y_prime_probe.tsv', sep='\t')

# Filter for reads with Y prime changes
changed = yprimes[yprimes['delta_y_prime_sign'].isin(['+', '-'])]

# For each changed read, identify which labeled region it maps to
# This would require overlap analysis with the alignment BAM file
# See existing scripts like y_prime_analysis.py for patterns
```

## Troubleshooting Common Issues

### Issue: No Y primes detected

**Possible causes:**
- Y primes genuinely absent on some chromosome ends (expected)
- BLAST parameters too stringent for cross-strain comparison

**Solutions:**
- Check raw BLAST results: `less pretelomeric_regions_7575_yprime_blast.txt`
- Relax BLAST parameters in `label_regions.sh`:
  ```bash
  MIN_PIDENT=70.0  # Lower from 75.0
  EVALUE=1e-3      # Higher from 1e-5
  ```

### Issue: Anchors missing for some chromosome ends

**Possible causes:**
- Reference quality issues
- Chromosome end not properly assembled
- Anchor sequences diverged significantly

**Solutions:**
- Check which chr ends are missing anchors:
  ```bash
  python scripts/visualize_labeled_regions.py \
      --input-tsv pretelomeric_regions_7575.tsv | \
      grep "without anchor"
  ```
- Inspect those regions manually in IGV
- Consider lowering MIN_PIDENT for anchors

### Issue: X primes too short or too long

**Expected behavior:**
- X prime lengths are inferred, not directly measured
- Typical range: 500bp - 5kb

**Solutions:**
- If too short (<500bp): Likely anchor and Y prime are very close
- If too long (>5kb): May indicate missing Y prime
- Manually verify regions in IGV
- Consider adjusting X prime inference logic in the Python script

## Advanced Customization

### Modify X Prime Inference

Edit `scripts/label_pretelomeric_regions.py`, function `infer_x_prime_regions()`:

```python
# Change default X prime length when no Y prime found
x_end = x_start + 5000  # Change from 3000 to 5000
```

### Add Custom Region Types

To add a new region type (e.g., spacer elements):

1. Add BLAST search in `label_pretelomeric_regions.py`
2. Parse results in `assign_regions_to_chr_ends()`
3. Update output writers to include new type

### Integrate with Existing RepeatMasker Pipeline

You already have RepeatMasker utilities in `scripts/repeatmasker_utils.py`. To integrate:

```bash
# Run RepeatMasker on extracted sequences
for region in anchor x_prime y_prime; do
    RepeatMasker -species saccharomyces \
        extracted_sequences/strain_7575_${region}.fasta
done

# Parse results with existing tools
python scripts/make_y_prime_repeatmasker_tsv.py ...
```

## Next Steps

After labeling regions, you can:

1. **Sequence Analysis**: BLAST labeled sequences against databases
2. **Comparative Genomics**: Compare region organization between strains
3. **Recombination Analysis**: Map recombination breakpoints to labeled regions
4. **Functional Studies**: Test effects of specific regions experimentally
5. **Publication**: Use labeled reference for figures and supplementary data

## File Organization

Recommended directory structure after running all pipelines:

```
teloTracker_snakemake/
├── create_ref.sh                    # Reference creation pipeline
├── label_regions.sh                 # Region labeling pipeline (NEW)
├── LABELING_PIPELINE_README.md      # Documentation (NEW)
├── WORKFLOW_EXAMPLE.md              # This file (NEW)
├── references/
│   ├── test_anchors.fasta
│   └── repeatmasker_6991_all_y_primes.fasta
├── scripts/
│   ├── label_pretelomeric_regions.py      # Main labeling script (NEW)
│   ├── visualize_labeled_regions.py       # Visualization script (NEW)
│   ├── extract_labeled_sequences.py       # Sequence extraction (NEW)
│   └── [existing scripts...]
└── results/
    └── dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection/
        └── assembly_7575/
            ├── assembly_7575_dorado_reference.fasta
            ├── assembly_7575_final_alignment_to_dorado_ref.bam
            └── pretelomeric_labels/                        # (NEW)
                ├── pretelomeric_regions_7575.gff3
                ├── pretelomeric_regions_7575.bed
                ├── pretelomeric_regions_7575.tsv
                ├── pretelomeric_regions_7575_anchor_blast.txt
                └── pretelomeric_regions_7575_yprime_blast.txt
```

---

**Questions or issues?** Check the LABELING_PIPELINE_README.md or open an issue in the repository.
