# Pre-Telomeric Region Labeling Pipeline

A comprehensive pipeline for identifying and annotating anchor, X prime, and Y prime regions in yeast subtelomeric sequences.

## 🚀 Quick Start

```bash
# 1. Validate setup
bash scripts/test_labeling_pipeline.sh

# 2. Run the pipeline (after create_ref.sh completes)
bash label_regions.sh

# 3. View results
python scripts/visualize_labeled_regions.py \
    --input-tsv results/*/assembly_*/pretelomeric_labels/pretelomeric_regions_*.tsv
```

**See [QUICK_START_LABELING.md](QUICK_START_LABELING.md) for detailed quick start guide.**

## 📋 Overview

This pipeline extends your TeloTracker workflow by adding automated labeling of pre-telomeric regions:

- **Anchor regions**: Sequences closest to centromere (~1-3 kb)
- **Y prime regions**: Variable subtelomeric elements (0-2 per chromosome end)
- **X prime regions**: Inferred regions between anchor and telomere

### How It Works

1. Uses BLAST to align known sequences from reference strain (6991) to your new strain
2. Filters alignments by quality (percent identity, coverage)
3. Assigns regions to specific chromosome ends
4. Infers X prime positions based on anchor and Y prime locations
5. Outputs annotations in multiple formats (GFF3, BED, TSV)

## 📁 Files Created

### Main Scripts
- **[label_regions.sh](label_regions.sh)** - Main pipeline script (run after create_ref.sh)
- **[scripts/label_pretelomeric_regions.py](scripts/label_pretelomeric_regions.py)** - Core labeling logic
- **[scripts/visualize_labeled_regions.py](scripts/visualize_labeled_regions.py)** - Summary and visualization
- **[scripts/extract_labeled_sequences.py](scripts/extract_labeled_sequences.py)** - Sequence extraction
- **[scripts/test_labeling_pipeline.sh](scripts/test_labeling_pipeline.sh)** - Validation test

### Documentation
- **[QUICK_START_LABELING.md](QUICK_START_LABELING.md)** - 5-minute quick start
- **[LABELING_PIPELINE_README.md](LABELING_PIPELINE_README.md)** - Comprehensive documentation
- **[WORKFLOW_EXAMPLE.md](WORKFLOW_EXAMPLE.md)** - Complete workflow with examples
- **[LABELING_PIPELINE_SUMMARY.txt](LABELING_PIPELINE_SUMMARY.txt)** - Pipeline overview
- **README_LABELING.md** (this file) - Main entry point

## 🔧 Prerequisites

- ✅ Completed `create_ref.sh` pipeline
- ✅ BLAST+ installed (`blastn` command)
- ✅ Python 3 with pandas
- ✅ Reference files in `references/` directory:
  - `test_anchors.fasta` (32 sequences)
  - `repeatmasker_6991_all_y_primes.fasta` (23 sequences)

Run the validation test to check your setup:
```bash
bash scripts/test_labeling_pipeline.sh
```

## 📊 Pipeline Workflow

```
create_ref.sh (Dorado polished reference)
    ↓
label_regions.sh
    ├── BLAST anchors
    ├── BLAST Y primes
    ├── Filter and assign regions
    └── Infer X primes
    ↓
Output Files:
    ├── pretelomeric_regions_${STRAIN_ID}.gff3  (genome browser)
    ├── pretelomeric_regions_${STRAIN_ID}.bed   (bedtools)
    ├── pretelomeric_regions_${STRAIN_ID}.tsv   (analysis)
    └── *_blast.txt                             (raw BLAST)
```

## 💻 Usage Examples

### Basic Usage

```bash
# Run the pipeline
bash label_regions.sh
```

### Visualize Results

```bash
# Summary statistics
python scripts/visualize_labeled_regions.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv

# Export summary table
python scripts/visualize_labeled_regions.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --export-summary summary_table.tsv
```

### Extract Sequences

```bash
# Extract all regions
python scripts/extract_labeled_sequences.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --reference results/.../assembly_7575_dorado_reference.fasta \
    --output-dir extracted_sequences \
    --prefix strain_7575

# Extract only anchors
python scripts/extract_labeled_sequences.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --reference results/.../assembly_7575_dorado_reference.fasta \
    --output-dir extracted_sequences \
    --prefix strain_7575 \
    --region-types anchor

# Extract specific chromosome ends
python scripts/extract_labeled_sequences.py \
    --input-tsv results/.../pretelomeric_regions_7575.tsv \
    --reference results/.../assembly_7575_dorado_reference.fasta \
    --output-dir extracted_sequences \
    --prefix strain_7575 \
    --chr-ends chr1L chr1R
```

### Load in IGV

1. `Genomes > Load Genome from File` → Select reference FASTA
2. `File > Load from File` → Select GFF3 file
3. Navigate to chromosome ends

## 📖 Documentation

| Document | Purpose |
|----------|---------|
| **[QUICK_START_LABELING.md](QUICK_START_LABELING.md)** | Get started in 5 minutes |
| **[LABELING_PIPELINE_README.md](LABELING_PIPELINE_README.md)** | Full documentation, configuration, troubleshooting |
| **[WORKFLOW_EXAMPLE.md](WORKFLOW_EXAMPLE.md)** | Complete examples and downstream analysis |
| **[LABELING_PIPELINE_SUMMARY.txt](LABELING_PIPELINE_SUMMARY.txt)** | Technical overview |

## ⚙️ Configuration

Edit `label_regions.sh` to match your `create_ref.sh` settings:

```bash
# Must match create_ref.sh
BASE_NAME="dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection"
STRAIN_ID="7575"

# BLAST parameters (adjust for your comparison)
MIN_PIDENT=75.0     # Minimum percent identity (75% for cross-strain)
MIN_LENGTH=100      # Minimum alignment length
EVALUE=1e-5        # E-value threshold
```

### Parameter Guidelines

| Scenario | MIN_PIDENT | EVALUE | Notes |
|----------|-----------|---------|-------|
| Same strain | 95.0 | 1e-20 | High stringency |
| Closely related | 85.0 | 1e-10 | Moderate stringency |
| Cross-strain (default) | 75.0 | 1e-5 | Relaxed for divergence |
| Distantly related | 70.0 | 1e-3 | Very relaxed |

## 📈 Expected Results

For a typical yeast genome (16 chromosomes = 32 chromosome ends):

| Region Type | Expected Count | Notes |
|-------------|---------------|-------|
| Anchors | 30-32 | Should find most/all |
| Y primes | 10-20 | Varies by strain |
| X primes | 30-32 | Inferred from anchors |

## 🔍 Output Formats

### GFF3 Format
Standard genome annotation with attributes:
- ID, Name, chr_end, source_seq
- Quality metrics: pident, bitscore, query_coverage
- Load in IGV, JBrowse, Artemis

### BED Format
6-column BED:
- chr, start, end, name, score, strand
- Compatible with bedtools suite

### TSV Format
Tabular data with columns:
- chr_end, chr, start, end, length, strand
- region_type, source
- pident, bitscore, query_coverage

## 🔬 Integration with Existing Pipeline

This pipeline complements your existing tools:

```
TeloTracker Pipeline:
├── Snakefile (y_prime_analysis)
├── create_ref.sh (assembly & polishing)
├── label_regions.sh (NEW - region labeling)
└── Your analysis scripts:
    ├── y_prime_analysis.py
    ├── get_stats_of_recombination.py
    ├── repeatmasker_utils.py
    └── plot_switch_distances_*.py
```

You can now:
- Map recombination events to specific labeled regions
- Analyze subtelomeric structure variation
- Compare region organization across strains
- Guide repeat element classification

## 🐛 Troubleshooting

| Issue | Solution |
|-------|----------|
| "Reference FASTA not found" | Run `create_ref.sh` first |
| "blastn: command not found" | Install BLAST+: `conda install blast` |
| No Y primes found | Normal for some chr ends; check BLAST output |
| Low anchor count | Try lowering MIN_PIDENT to 70 |
| Module not found error | Install pandas: `pip install pandas` |

Run the validation test for detailed diagnostics:
```bash
bash scripts/test_labeling_pipeline.sh
```

## 📚 Examples from Documentation

### Analyze Region Sizes (Python)

```python
import pandas as pd

df = pd.read_csv('pretelomeric_regions_7575.tsv', sep='\t')

# Region size statistics
for region_type in ['anchor', 'x_prime', 'y_prime']:
    subset = df[df['region_type'] == region_type]
    print(f"\n{region_type}:")
    print(f"  Mean: {subset['length'].mean():.0f} bp")
    print(f"  Median: {subset['length'].median():.0f} bp")
    print(f"  Range: {subset['length'].min()}-{subset['length'].max()} bp")
```

### Extract with bedtools

```bash
# Get all anchor sequences
bedtools getfasta \
    -fi reference.fasta \
    -bed pretelomeric_regions_7575.bed \
    -name -s | \
    grep -A1 "anchor" > anchor_sequences.fasta
```

### Compare Strains

```bash
# Generate summaries for multiple strains
for strain in 6991 7575 7093; do
    python scripts/visualize_labeled_regions.py \
        --input-tsv results/.../pretelomeric_regions_${strain}.tsv \
        --export-summary summary_${strain}.tsv
done

# Compare in R or Python
```

## 🎯 Next Steps

After running the labeling pipeline:

1. **Visualize** in IGV to validate region assignments
2. **Extract** sequences for specific regions of interest
3. **Analyze** with your existing scripts (y_prime_analysis.py, etc.)
4. **Compare** across multiple strains
5. **Map** recombination breakpoints to labeled regions
6. **Publish** results with labeled reference

## 📝 Citation

If you use this pipeline, please cite:
- BLAST: Altschul et al. (1990)
- Your TeloTracker publication (when available)

## 🤝 Contributing

Found a bug or have a suggestion? Please open an issue in the repository.

## 📄 License

This pipeline is part of the TeloTracker suite. Check repository for license information.

---

**Created:** January 2026
**Version:** 1.0
**Maintainer:** Your Lab/Team Name

For questions, see the documentation or open an issue.
