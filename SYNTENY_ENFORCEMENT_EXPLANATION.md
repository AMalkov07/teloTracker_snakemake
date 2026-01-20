# Synteny Enforcement in Anchor Detection

## Overview

The `label_pretelomeric_regions.py` script has been updated to **enforce strict chromosome matching** and provide **clear quality warnings** for anchor and Y prime detection.

## What Changed

### 1. Strict Chromosome Matching (Synteny Enforcement)

**Previous Behavior:**
- BLAST searched the entire genome for anchor matches
- `chr1L_anchor` could match **any chromosome** in the new reference
- Results were labeled by the query name, not actual position
- Could lead to misleading labels if chromosomes were rearranged

**New Behavior:**
- `chr1L_anchor` **must match** a chromosome named `chr1*` (e.g., `chr1_extended`)
- Mismatches are **detected and excluded** from results
- **Detailed quality report** shows what was excluded and why

### 2. Quality Threshold for Anchors (95% Identity)

**Expectation:**
- Anchors should be highly conserved between closely related strains
- ≥95% sequence identity expected for proper matches

**Implementation:**
- Anchors with <95% identity trigger a **⚠️ WARNING**
- These anchors are still **included** in results but flagged for manual review
- Suggests potential sequence divergence or assembly issues

### 3. Quality Reporting

**Console Output:**
Clear, color-coded warnings during pipeline execution:
- ⚠️ **CRITICAL**: Chromosome mismatches (excluded)
- ⚠️ **WARNING**: Low quality anchors <95% (included with warning)
- ⚠️ **WARNING**: Missing anchors (not found)
- ⚠️ **INFO**: Y prime issues (expected variation)
- ✅ **SUCCESS**: No issues detected

**File Output:**
- `pretelomeric_regions_7575_quality_report.txt` - Full report for review
- Contains all mismatches, low-quality hits, and missing anchors

## How Chromosome Matching Works

### Function: `check_chromosome_match(expected_chr, subject_chr)`

**Purpose:** Match chromosome names accounting for variations

**Examples:**

| Query | Expected | Subject | Match? | Reason |
|-------|----------|---------|--------|--------|
| chr1L_anchor | chr1 | chr1_extended | ✅ Yes | chr1_extended starts with chr1 |
| chr1L_anchor | chr1 | chr1_7575 | ✅ Yes | chr1_7575 starts with chr1 |
| chr1L_anchor | chr1 | chr10_extended | ❌ No | Would match chr1 with chr10 (prevented!) |
| chr2R_anchor | chr2 | chr2_extended | ✅ Yes | Correct match |
| chr5L_anchor | chr5 | chr3_extended | ❌ No | Wrong chromosome |

**Algorithm:**
```python
def check_chromosome_match(expected_chr, subject_chr):
    # Handles naming like chr1_extended, chr1_7575, etc.
    if subject_lower.startswith(expected_lower):
        # Ensure chr1 doesn't match chr10, chr11, etc.
        next_char = subject_lower[len(expected_lower)]
        return not next_char.isdigit()  # Next char must be non-digit
    return False
```

**Why this works for your reference:**
Your reference uses names like:
- `chr1_extended`, `chr2_extended`, ..., `chr16_extended`
- The underscore after the number ensures clean matching
- `chr1` matches `chr1_extended` ✅
- `chr1` does NOT match `chr10_extended` ✅

## Quality Report Sections

### 1. CRITICAL: Anchor Chromosome Mismatches

**What it means:**
An anchor expected to be on one chromosome was found on a different chromosome.

**Example:**
```
⚠️  CRITICAL: ANCHOR CHROMOSOME MISMATCHES DETECTED ⚠️
The following anchors matched the WRONG chromosome:
--------------------------------------------------------------------------------
  chr5L_anchor          → Expected: chr5        | Found: chr3_extended
                         Identity: 88.5% | Coverage: 92.3% | Bitscore: 1456.2
```

**Possible causes:**
- Chromosomal rearrangement between strains
- Anchor duplication/translocation
- Assembly error
- Mislabeled reference sequences

**Action:**
- These hits are **EXCLUDED** from results
- Manually inspect in IGV to verify
- Consider if biological rearrangement is expected

### 2. WARNING: Low Quality Anchor Matches (<95%)

**What it means:**
An anchor matched the correct chromosome but with lower than expected identity.

**Example:**
```
⚠️  WARNING: LOW QUALITY ANCHOR MATCHES (<95% identity) ⚠️
The following anchors have LOWER than expected sequence identity:
--------------------------------------------------------------------------------
  chr3R_anchor          → chr3_extended    | Position:  5,234 - 7,123
                         Identity: 89.2% | Coverage: 94.5% | Bitscore: 1234.5
```

**Possible causes:**
- Natural sequence divergence between strains
- Indels or SNPs in anchor region
- Partial assembly of anchor region
- Reference quality issues

**Action:**
- These anchors are **INCLUDED** in results
- Review manually in IGV
- Consider if divergence is expected for your strains
- May indicate interesting evolutionary changes

### 3. WARNING: Missing Anchors

**What it means:**
No anchor was found for certain chromosome ends.

**Example:**
```
⚠️  WARNING: MISSING ANCHORS ⚠️
No anchors found for 4 chromosome ends:
--------------------------------------------------------------------------------
  chr7L   chr7R   chr14L  chr14R
```

**Possible causes:**
- These regions not assembled in reference
- Anchors diverged >25% (beyond detection threshold)
- True biological loss of anchor sequences
- BLAST parameters too stringent

**Action:**
- Check if these chromosome ends are properly assembled
- Try lowering `--min-pident` threshold (e.g., to 70%)
- Check raw BLAST output for near-misses
- May need manual curation

### 4. INFO: Y Prime Mismatches/Low Quality

**What it means:**
Y primes are more variable and can be translocated, so these are informational.

**Action:**
- Y prime mismatches are excluded (like anchors)
- Low quality Y primes (<80%) are included
- This variation is expected biological diversity

## Expected Results for Your Reference

For `assembly_7575_dorado_reference.fasta`:

**Chromosomes:** chr1_extended through chr16_extended (16 chromosomes)

**Expected:**
- 32 chromosome ends total (16 × 2)
- ~30-32 anchors found (most or all chr ends)
- 0-20 Y primes (varies by strain)
- ~30-32 X primes (inferred from anchors)

**Good Result:**
```
✅ No quality issues detected!
All anchors found with high identity (≥95%) on expected chromosomes.
```

**Typical Result:**
```
⚠️  WARNING: MISSING ANCHORS ⚠️
No anchors found for 2 chromosome ends: chr14L, chr14R

⚠️  WARNING: LOW QUALITY ANCHOR MATCHES (<95% identity) ⚠️
3 anchor(s) with 90-95% identity (sequence divergence)
```

## Running the Pipeline

The pipeline automatically runs with synteny enforcement:

```bash
bash label_regions.sh
```

The script will:
1. Run BLAST
2. Filter hits by chromosome matching
3. Check quality thresholds
4. Print detailed quality report to console
5. Save quality report to file
6. Generate annotation files (GFF3, BED, TSV)

## Output Files

After running, you'll have:

```
pretelomeric_labels/
├── pretelomeric_regions_7575.gff3            # Genome annotations
├── pretelomeric_regions_7575.bed             # BED format
├── pretelomeric_regions_7575.tsv             # Tabular data
├── pretelomeric_regions_7575_anchor_blast.txt    # Raw BLAST
├── pretelomeric_regions_7575_yprime_blast.txt    # Raw BLAST
└── pretelomeric_regions_7575_quality_report.txt  # Quality report ← NEW!
```

## Interpreting Results

### Scenario 1: Perfect Match
```
✅ All anchors found with ≥95% identity on correct chromosomes
   → Strains are very similar, high confidence in results
```

### Scenario 2: Some Low Quality
```
⚠️ 5 anchors with 90-95% identity
   → Some divergence, but reasonable for cross-strain comparison
   → Manually verify in IGV
```

### Scenario 3: Mismatches Detected
```
⚠️ CRITICAL: 2 anchors matched wrong chromosomes
   → Potential rearrangement or assembly issue
   → Investigate in IGV
   → Compare with known karyotypes
```

### Scenario 4: Many Missing
```
⚠️ 10+ anchors missing
   → Assembly may be incomplete
   → Try lowering stringency (--min-pident 70)
   → Check reference quality
```

## Troubleshooting

### Problem: Too many missing anchors

**Solution 1: Lower minimum identity**
Edit `label_regions.sh`:
```bash
MIN_PIDENT=70.0  # Lower from 75.0
```

**Solution 2: Check raw BLAST output**
```bash
less pretelomeric_regions_7575_anchor_blast.txt
# Look for hits that are close to threshold
```

**Solution 3: Check reference quality**
```bash
# Load reference in IGV
# Navigate to missing chromosome ends
# Check for assembly gaps or quality issues
```

### Problem: Many low quality matches (<95%)

**This may be normal** if:
- Strains are distantly related
- Known sequence polymorphisms exist
- Different genetic backgrounds

**This may be a problem** if:
- Expected closely related strains
- Assembly quality is poor
- Reference has many N's (gaps)

**Action:**
- Review each low-quality anchor in IGV
- Check alignment in BLAST output
- Consider if divergence is biologically meaningful

### Problem: Chromosome mismatches detected

**Possible explanations:**
1. **True biological rearrangement** - Validate with other methods
2. **Assembly error** - Check assembly quality at that locus
3. **Anchor duplication** - Ancient or recent duplication event
4. **Naming mismatch** - Check chromosome naming conventions

**Action:**
```bash
# Check the mismatched region in IGV
# Load both 6991 and 7575 references
# Compare synteny manually
```

## Advanced: Disabling Synteny Enforcement

If you want to see **all** BLAST hits regardless of chromosome (old behavior):

Edit `label_pretelomeric_regions.py` line 649-653:
```python
chr_end_regions, quality_report = assign_regions_to_chr_ends(
    anchor_df, yprime_df,
    enforce_synteny=False,  # ← Change to False
    min_high_quality_pident=95.0
)
```

This will:
- Keep all BLAST hits
- Still report mismatches in quality report
- Useful for detecting rearrangements

## Summary

**Key Benefits:**
1. ✅ **Accuracy**: Anchors only assigned to correct chromosomes
2. ✅ **Quality control**: Clear warnings for low-quality matches
3. ✅ **Transparency**: Detailed reporting of what was excluded and why
4. ✅ **Flexibility**: Can adjust thresholds based on your strains
5. ✅ **Documentation**: Quality report file for later review

**Expected Behavior:**
- High-quality anchors (≥95% identity) on correct chromosomes: **INCLUDED**
- Low-quality anchors (80-95%) on correct chromosomes: **INCLUDED with WARNING**
- Any anchor on wrong chromosome: **EXCLUDED with CRITICAL WARNING**
- Missing anchors: **REPORTED as WARNING**

This ensures your labeled reference has high-confidence annotations while giving you full visibility into any potential issues!
