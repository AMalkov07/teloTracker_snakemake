# Anchor-Only Extraction Test Guide

## Overview

The script has been modified to test **ONLY anchor extraction**, with Y prime and X prime detection commented out.

## What Was Changed

### Commented Out:
1. ✅ Y prime BLAST search
2. ✅ Y prime parsing
3. ✅ X prime inference
4. ✅ Y prime/X prime output in summary

### Still Active:
1. ✅ Anchor BLAST search
2. ✅ Anchor parsing and filtering
3. ✅ Strict chromosome matching (synteny enforcement)
4. ✅ Quality threshold checking (95% identity)
5. ✅ Quality report generation
6. ✅ GFF3/BED/TSV output (anchors only)

## How to Run the Test

### Option 1: Use the Test Script (Recommended)

```bash
bash test_anchor_only.sh
```

This script:
- Uses 8 threads (faster than default)
- Uses relaxed BLAST parameters (75% identity, 1e-5 E-value)
- Creates output in `results/.../test_anchor_only/`
- Shows clear test mode messages

### Option 2: Direct Python Call

```bash
python scripts/label_pretelomeric_regions.py \
    --reference results/dorado_fast5_7575_day0_PromethION_no_tag_yes_rejection/assembly_7575/assembly_7575_dorado_reference.fasta \
    --anchors references/test_anchors.fasta \
    --yprimes references/repeatmasker_6991_all_y_primes.fasta \
    --output-dir test_output \
    --prefix test_anchors \
    --threads 8 \
    --min-pident 75.0 \
    --min-length 100 \
    --evalue 1e-5
```

Note: `--yprimes` is still required by argparse but won't be used

## Expected Console Output

```
================================================================================
Pre-telomeric Region Labeling Pipeline
================================================================================
Reference: results/.../assembly_7575_dorado_reference.fasta
Anchors: references/test_anchors.fasta
Y primes: references/repeatmasker_6991_all_y_primes.fasta
Output directory: test_output

Step 1: Running BLAST for anchor regions...
Running BLAST: blastn -query references/test_anchors.fasta -subject ...
Anchor BLAST results: test_output/test_anchors_anchor_blast.txt

Step 2: Parsing and filtering BLAST results (ANCHORS ONLY - TESTING MODE)...
Found 32 anchor hits

Step 3: Assigning ANCHOR regions to chromosome ends (with synteny enforcement)...
         - Enforcing strict chromosome matching
         - Requiring ≥95% identity for high-quality anchors
Detected regions for 32 chromosome ends

Summary by chromosome end (ANCHORS ONLY - TESTING MODE):
  chr1L: anchor=1
  chr1R: anchor=1
  chr2L: anchor=1
  chr2R: anchor=1
  chr3L: anchor=1
  chr3R: anchor=1
  chr4L: anchor=1
  chr4R: anchor=1
  chr5L: anchor=1
  chr5R: anchor=1
  chr6L: anchor=1
  chr6R: anchor=1
  chr7L: anchor=1
  chr7R: anchor=1
  chr8L: anchor=1
  chr8R: anchor=1
  chr9L: anchor=1
  chr9R: anchor=1
  chr10L: anchor=1
  chr10R: anchor=1
  chr11L: anchor=1
  chr11R: anchor=1
  chr12L: anchor=1
  chr12R: anchor=1
  chr13L: anchor=1
  chr13R: anchor=1
  chr14L: anchor=1
  chr14R: anchor=1
  chr15L: anchor=1
  chr15R: anchor=1
  chr16L: anchor=1
  chr16R: anchor=1

================================================================================
QUALITY REPORT
================================================================================

✅ No quality issues detected!
All anchors found with high identity (≥95%) on expected chromosomes.

================================================================================

Step 6: Writing output files...
GFF3 file: test_output/test_anchors.gff3
BED file: test_output/test_anchors.bed
TSV file: test_output/test_anchors.tsv
Quality report: test_output/test_anchors_quality_report.txt

================================================================================
Pipeline completed successfully!
================================================================================
```

## Expected Output Files

In `results/.../test_anchor_only/` (or your specified output directory):

### 1. `test_anchors_7575_anchor_blast.txt`
**Raw BLAST output** (14 columns, tab-separated)

Example lines:
```
chr1L_anchor	chr1_extended	96.5	1823	45	3	1	1823	5234	7056	0.0	2456.7	1823	230218
chr1R_anchor	chr1_extended	95.8	1756	52	2	1	1756	228450	230205	0.0	2389.3	1756	230218
chr2L_anchor	chr2_extended	97.2	1645	38	2	1	1645	3456	5100	0.0	2234.1	1645	813184
```

Columns:
- Query ID (e.g., chr1L_anchor)
- Subject ID (e.g., chr1_extended)
- % identity (e.g., 96.5)
- Alignment length (e.g., 1823)
- Mismatches, gaps
- Query start/end, Subject start/end
- E-value, Bitscore
- Query/Subject lengths

### 2. `test_anchors_7575.tsv`
**Tabular annotation data** (ONLY anchors, no Y primes or X primes)

Example:
```tsv
chr_end	chr	start	end	length	strand	region_type	source	pident	bitscore	query_coverage
chr1L	chr1_extended	5234	7056	1823	+	anchor	chr1L_anchor	96.5	2456.7	100.0
chr1R	chr1_extended	228450	230205	1756	+	anchor	chr1R_anchor	95.8	2389.3	100.0
chr2L	chr2_extended	3456	5100	1645	+	anchor	chr2L_anchor	97.2	2234.1	100.0
```

**What to look for:**
- **32 rows** (one per chromosome end: chr1L-chr16R)
- **region_type column**: Only "anchor" (no y_prime or x_prime)
- **pident column**: Should mostly be ≥95% (high quality)
- **chr column**: Should match expected (chr1_extended for chr1L/R, etc.)

### 3. `test_anchors_7575.gff3`
**Genome annotation format** (ONLY anchors)

Example:
```gff
##gff-version 3
chr1_extended	TeloTracker	anchor	5234	7056	.	+	.	ID=anchor_1;Name=chr1L_anchor;chr_end=chr1L;source_seq=chr1L_anchor;pident=96.50;bitscore=2456.70;query_coverage=100.00
chr1_extended	TeloTracker	anchor	228450	230205	.	+	.	ID=anchor_2;Name=chr1R_anchor;chr_end=chr1R;source_seq=chr1R_anchor;pident=95.80;bitscore=2389.30;query_coverage=100.00
```

**What to look for:**
- Only "anchor" features (column 3)
- No "y_prime" or "x_prime" features
- Positions within valid chromosome ranges

### 4. `test_anchors_7575.bed`
**BED format** (6 columns, ONLY anchors)

Example:
```bed
chr1_extended	5234	7056	chr1L_anchor	2457	+
chr1_extended	228450	230205	chr1R_anchor	2389	+
chr2_extended	3456	5100	chr2L_anchor	2234	+
```

**What to look for:**
- 32 lines (one per anchor)
- All names end with "_anchor"
- No y_prime or x_prime entries

### 5. `test_anchors_7575_quality_report.txt`
**Quality control report**

Best case (all good):
```
================================================================================
PRE-TELOMERIC REGION LABELING - QUALITY REPORT
================================================================================

No quality issues detected!
All anchors found with high identity (≥95%) on expected chromosomes.
```

With warnings:
```
================================================================================
PRE-TELOMERIC REGION LABELING - QUALITY REPORT
================================================================================

ANCHOR LOW QUALITY MATCHES (<95% identity, INCLUDED WITH WARNING):
--------------------------------------------------------------------------------
Query: chr3R_anchor
  Chromosome: chr3_extended
  Position: 5,234 - 7,123
  Identity: 89.20% | Coverage: 94.50% | Bitscore: 1234.50

MISSING ANCHORS:
--------------------------------------------------------------------------------
Chromosome ends without anchors: chr14L, chr14R
```

## What to Check in Results

### 1. Count the Anchors
```bash
# Should show 32 anchors (one per chromosome end)
grep -c "anchor" test_output/test_anchors_7575.tsv
```

Expected: **32**

### 2. Check Identity Values
```bash
# View all identity values
awk 'NR>1 {print $1, $9}' test_output/test_anchors_7575.tsv | sort
```

Expected:
- Most should be **95-100%**
- A few might be **85-95%** (strain divergence)
- Very few should be **<85%**

### 3. Check Chromosome Matching
```bash
# Verify chr ends match chromosomes
awk 'NR>1 {print $1, $2}' test_output/test_anchors_7575.tsv
```

Expected:
```
chr1L chr1_extended
chr1R chr1_extended
chr2L chr2_extended
chr2R chr2_extended
...
```

Each `chrNL` and `chrNR` should map to `chrN_extended`

### 4. Check for Mismatches
```bash
# Should be empty (no mismatches)
grep "CHROMOSOME MISMATCHES" test_output/test_anchors_7575_quality_report.txt
```

Expected: **Empty** (no mismatches found)

If you see mismatches, this indicates:
- Chromosomal rearrangement
- Assembly issue
- Mislabeled sequences

### 5. Check for Missing Anchors
```bash
grep "MISSING ANCHORS" test_output/test_anchors_7575_quality_report.txt
```

Expected: **Empty or minimal** (0-2 missing)

If many anchors missing:
- Reference assembly incomplete
- Try lowering `--min-pident` to 70.0
- Check those chr ends in IGV

## Interpreting Results

### Scenario 1: Perfect Result ✅
```
32 anchors found
All with ≥95% identity
All on correct chromosomes (chr1L → chr1_extended, etc.)
No quality warnings
```

**Meaning:** Excellent! Strains are very similar, high-quality assembly.

### Scenario 2: Good Result with Minor Warnings ⚠️
```
30-32 anchors found
Most with ≥95% identity, 2-3 with 90-95%
All on correct chromosomes
```

**Meaning:** Good! Some sequence divergence expected. Manually verify low-quality ones.

### Scenario 3: Concerning Result 🔴
```
<28 anchors found
Several with <90% identity
Some chromosome mismatches detected
```

**Meaning:**
- Check assembly quality
- Verify reference file integrity
- May indicate significant rearrangement or assembly errors
- Investigate in IGV

## Troubleshooting

### Problem: No BLAST output
**Check:**
```bash
ls -lh test_output/test_anchors_7575_anchor_blast.txt
```

If empty (0 bytes):
- BLAST found no hits
- Lower `MIN_PIDENT` to 70.0
- Check reference file format

### Problem: Too many missing anchors
**Try:**
```bash
# Edit test_anchor_only.sh
MIN_PIDENT=70.0  # Lower from 75.0

# Re-run
bash test_anchor_only.sh
```

### Problem: Unexpected chromosome names
**Check:**
```bash
# View reference chromosome names
grep "^>" results/.../assembly_7575_dorado_reference.fasta
```

Expected: `chr1_extended`, `chr2_extended`, etc.

If different, the chromosome matching logic may need adjustment.

## Next Steps After Testing

1. **Review results** - Check all quality metrics look good
2. **Verify in IGV** - Load GFF3 and visually confirm anchor positions
3. **If satisfied** - Uncomment Y prime and X prime sections
4. **Run full pipeline** - Use `label_regions.sh` for complete labeling

## Uncommenting Full Pipeline

When ready to enable Y prime and X prime:

1. Edit `scripts/label_pretelomeric_regions.py`
2. Remove the `##` comment markers from lines:
   - ~625-631: Y prime BLAST
   - ~638-641: Y prime parsing
   - ~662-665: X prime inference
   - ~671-674: Summary output
3. Delete the empty yprime_df creation line (643-644)
4. Run full pipeline with `bash label_regions.sh`

Or restore from git if you committed before testing.

## Summary

**What You're Testing:**
- ✅ BLAST finds anchor sequences
- ✅ Synteny enforcement works (chr1L → chr1_extended)
- ✅ Quality thresholds work (95% identity warning)
- ✅ Output files generated correctly
- ✅ Quality report accurate

**What You're NOT Testing:**
- ❌ Y prime detection
- ❌ X prime inference
- ❌ Multi-region output

This focused test ensures the core anchor detection logic works correctly before adding complexity!
