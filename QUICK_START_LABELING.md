# Quick Start: Labeling Pre-Telomeric Regions

Get your reference genome labeled with anchor, X prime, and Y prime regions in 5 minutes.

## Prerequisites

✅ Completed `create_ref.sh` successfully
✅ BLAST+ installed (`blastn` command available)
✅ Python 3 with pandas installed

## Step 1: Run the Labeling Pipeline (1 command)

```bash
bash label_regions.sh
```

Or submit as a job:

```bash
qsub label_regions.sh
```

**What it does:**
- BLAST anchors from strain 6991 against your new reference
- BLAST Y primes from strain 6991 against your new reference
- Infers X prime regions between anchors and Y primes
- Generates GFF3, BED, and TSV annotation files

**Time:** ~5-15 minutes depending on reference size

## Step 2: View Results

### Quick summary:

```bash
cd results/dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection/assembly_7575/pretelomeric_labels/

# View summary
python ../../../../scripts/visualize_labeled_regions.py \
    --input-tsv pretelomeric_regions_7575.tsv
```

### Output files:

- **pretelomeric_regions_7575.gff3** - Load in IGV or other genome browsers
- **pretelomeric_regions_7575.bed** - Use with bedtools
- **pretelomeric_regions_7575.tsv** - Open in Excel or analyze with Python/R

## Step 3: Load in IGV (Optional)

1. Open IGV
2. `Genomes > Load Genome from File` → Select your reference FASTA
3. `File > Load from File` → Select the GFF3 file
4. Navigate to chromosome ends to see labeled regions

## Common Adjustments

### If you get too few hits:

Edit `label_regions.sh` and relax parameters:

```bash
MIN_PIDENT=70.0    # Lower from 75.0 (allow more divergence)
EVALUE=1e-3        # Higher from 1e-5 (less stringent)
```

### If you want to extract sequences:

```bash
python scripts/extract_labeled_sequences.py \
    --input-tsv pretelomeric_regions_7575.tsv \
    --reference assembly_7575_dorado_reference.fasta \
    --output-dir extracted_sequences \
    --prefix my_strain
```

## Expected Results

For a typical yeast genome (16 chromosomes = 32 chr ends):
- **Anchors:** ~30-32 (should find most/all)
- **Y primes:** ~10-20 (varies by strain)
- **X primes:** ~30-32 (inferred from anchors)

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Reference FASTA not found" | Run `create_ref.sh` first |
| "blastn: command not found" | Install BLAST+: `conda install blast` |
| No Y primes found | Normal for some chr ends; check BLAST results |
| Script fails with "No module named 'pandas'" | Install pandas: `pip install pandas` |

## Next Steps

- View detailed docs: [LABELING_PIPELINE_README.md](LABELING_PIPELINE_README.md)
- See complete workflow: [WORKFLOW_EXAMPLE.md](WORKFLOW_EXAMPLE.md)
- Analyze results with your existing scripts in `scripts/`

## File Locations

```
results/
└── ${BASE_NAME}/
    └── assembly_${STRAIN_ID}/
        ├── assembly_${STRAIN_ID}_dorado_reference.fasta  ← Your reference
        └── pretelomeric_labels/                           ← NEW labels here
            ├── pretelomeric_regions_${STRAIN_ID}.gff3
            ├── pretelomeric_regions_${STRAIN_ID}.bed
            └── pretelomeric_regions_${STRAIN_ID}.tsv
```

---

**Questions?** See full documentation in [LABELING_PIPELINE_README.md](LABELING_PIPELINE_README.md)
