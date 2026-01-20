# Y Prime Detection with Fragment Merging - Complete Guide

## Overview

Y prime detection is significantly more complex than anchor detection because:
1. **Variable number**: 0 to 10+ Y primes per chromosome end
2. **Tandem repeats**: Multiple adjacent Y primes (often identical)
3. **Fragmentation**: BLAST may find Y primes as 2-3 large fragments
4. **No chromosome constraint**: Any Y prime can occur on any chromosome end
5. **Size variation**: Short (5-6.4kb) vs Long (6.4-6.9kb) Y primes

## Key Challenges Solved

### Challenge 1: Distinguishing Real Fragments from Repetitive Noise

**Problem**: BLAST finds hundreds of tiny repetitive matches (96-500bp) within Y primes.

**Solution**: Filter for substantial hits only (>2.5kb, >80% identity, >50% coverage).

```
Before filtering: 8,519 BLAST hits
After filtering:    682 high-quality hits
Final Y primes:      37 complete regions
```

### Challenge 2: Merging Fragments Without Over-Merging

**Problem**: User wanted to merge 2-3 real fragments (e.g., two 3kb pieces → one 6kb Y prime), but naive merging grouped hundreds of tiny hits.

**Solution**: Strict merging rules:
- Maximum 3 fragments per merge
- Minimum 2.5kb per fragment
- Maximum 500bp gap between fragments
- Only merge if result is valid Y prime size

### Challenge 3: Tandem Repeat Detection

**Problem**: Some chromosome ends have 10+ nearly identical Y primes in tandem array.

**Solution**: Keep overlapping Y primes if they're:
- Spatially separated (>50% non-overlapping)
- High quality (>95% identity)
- Correct size (valid Y prime length)

## Detection Strategy

### Step 1: BLAST Search
```bash
blastn -query repeatmasker_6991_all_y_primes.fasta \
       -subject assembly_7575_reference.fasta \
       -evalue 1e-5 \
       -outfmt "6 qseqid sseqid pident length ..."
```

**Parameters**:
- E-value: 1e-5 (relaxed for cross-strain comparison)
- Min identity: 80% (Y primes more variable than anchors)
- Min coverage: 50% (allows fragment detection)

### Step 2: Filter High-Quality Hits

**Criteria**:
- Length ≥2.5kb (real fragments, not repetitive noise)
- Identity ≥80% (strain variation allowed)
- Coverage ≥50% OR length ≥4kb (allow large partial hits)

**Result**: 682 hits from 8,519 total (92% of noise filtered out)

### Step 3: Merge Adjacent Fragments

**Merging Rules**:
1. Group hits on same chromosome within 500bp of each other
2. Limit to 2-3 fragments maximum
3. Only merge if total length = valid Y prime size (5-6.9kb)
4. Keep complete Y primes as-is (no merging needed)

**Example** (chr5L):
```
Fragment 1: 3,826bp (identity: 96.2%)
Fragment 2: 4,688bp (identity: 96.8%)
Gap: 147bp (telomeric repeats)
-------------------------------------------
Merged: 6,482bp Long Y prime ✓
```

### Step 4: Filter Overlapping Y Primes

**Overlap Resolution**:
- If two Y primes overlap >50%, keep the best one
- Ranking: complete size > merged > high bitscore > high identity

**Use Case**: Avoid double-counting same Y prime from different query sequences

### Step 5: Assign to Chromosome Ends

**Assignment Logic**:
- Y prime must be telomere-proximal to anchor
- **L arm**: Y prime BEFORE anchor (smaller coordinates)
- **R arm**: Y prime AFTER anchor (larger coordinates)

**Example** (chr4R):
```
Anchor:     1,436,850 - 1,441,890
Y primes:   1,447,498 - 1,515,041 (10 tandem Y primes) ✓
```

## Size Classification

### Short Y Primes (5,000 - 6,400 bp)
- Examples: chr2L (5,976bp), chr8L (5,057bp), chr12L (5,192bp)
- Typically ~6kb
- Less complex internal structure

### Long Y Primes (6,400 - 6,900 bp)
- Examples: chr4R (6,606bp), chr14L (6,654bp), chr16L (6,738bp)
- Typically ~6.6kb
- More complex, often contain additional elements

### Fragments (<5,000 or >6,900 bp)
- Too small or large to be complete Y primes
- May indicate:
  - Partial Y prime at contig boundary
  - Degraded/truncated Y prime
  - Merged with non-Y prime sequence

## Real-World Results (Strain 7575)

### Summary Statistics
```
Total Y primes:        37
Short Y primes:        10
Long Y primes:         27
Fragments:              0
Merged from fragments: 35 (95%)
```

### Tandem Repeat Arrays

**chr4R**: 10 tandem Long Y primes
```
Position              Length   Identity  Type
1,447,498-1,454,103   6,606bp  97.46%    Long (merged 3 fragments)
1,454,161-1,460,766   6,606bp  97.46%    Long (merged 3 fragments)
1,460,824-1,467,429   6,606bp  97.46%    Long (merged 3 fragments)
... 7 more identical Y primes ...
```
**Spacing**: ~6.8kb per Y prime (6.6kb Y prime + 200bp gap)

**chr12R**: 5 tandem Y primes (mix of Long and other types)

**chr8R**: 4 tandem Y primes

**chr14L**: 3 tandem Y primes (Short type)

### Solo Y Primes

**chr2L**: Single Short Y prime (5,976bp, 100% identity)
- Only Y prime that didn't require fragment merging!
- Perfect match to reference

**chr5L, chr5R, chr6L, chr7R, chr13L, chr14R, chr15L, chr15R, chr16L, chr16R**: Single Y primes per chr end

### Chr Ends Without Y Primes (12 total)

```
chr1L, chr1R, chr3L, chr3R, chr4L, chr6R, chr7L, chr9L, chr9R, chr10R, chr11L, chr11R
```

**This is normal and expected** - not all chr ends have Y primes. Y primes are variable elements.

## Fragment Merging Examples

### Example 1: Two-Fragment Merge (chr5L)
```
Before Merging:
  Fragment 1: 1-3,826    (3,826bp, 96.2% identity)
  Fragment 2: 3,973-6,669 (2,697bp actual, 4,688bp hit)
  Gap: 147bp

After Merging:
  Y prime: 187-6,669 (6,482bp Long Y prime)
  Identity: 96.5% (weighted average)
  Class: Long
```

### Example 2: Three-Fragment Merge (chr4R, typical)
```
Before Merging:
  Fragment 1: position 1-4,990   (4,990bp from query A)
  Fragment 2: position 1-4,990   (4,990bp from query B, overlapping)
  Fragment 3: position 1-4,570   (4,570bp from query C, overlapping)

After Merging:
  Y prime: 1,447,498-1,454,103 (6,606bp)
  Identity: 97.46%
  Note: Fragments overlap significantly, indicating tandem structure
```

### Example 3: Three-Fragment Merge (chr8L)
```
Before Merging:
  Fragment 1: 1-3,673    (3,673bp, Short Y prime hit)
  Fragment 2: 1-3,673    (3,673bp, same Short Y prime, different query)
  Fragment 3: 1-5,057    (5,057bp, longer hit covering both)
  Gap: minimal overlap

After Merging:
  Y prime: 192-5,249 (5,057bp Short Y prime)
  Identity: 99.57%
  Class: Short
```

## Quality Metrics

### High-Quality Detection Indicators

✅ **Identity 95-100%**: Excellent match, little divergence
- chr2L: 100% (perfect)
- chr14L: 99.86% (near perfect)
- chr12L: 99.57%

✅ **Identity 90-95%**: Good match, some strain variation
- chr10L: 96.12%
- chr15L: 96.02%
- chr16R: 96.72%

⚠️ **Identity 85-90%**: Acceptable, more divergence
- (None in this dataset - all Y primes >95%)

### Fragment Count Interpretation

- **1 fragment**: Complete Y prime found in single BLAST hit (rare, only chr2L)
- **2 fragments**: Y prime split in two pieces, cleanly merged
- **3 fragments**: Y prime detected by multiple overlapping query sequences

**Note**: 3 fragments is normal when different Y prime reference sequences match overlapping regions.

## Configuration Parameters

### In `yprime_detection_v2.py`:

```python
# Size thresholds
YPRIME_SHORT_MIN = 5000   # Minimum short Y prime
YPRIME_SHORT_MAX = 6400   # Maximum short Y prime
YPRIME_LONG_MIN = 6400    # Minimum long Y prime
YPRIME_LONG_MAX = 6900    # Maximum long Y prime

# Quality filters
MIN_YPRIME_IDENTITY = 80.0  # Minimum identity (%)
MIN_YPRIME_COVERAGE = 50.0  # Minimum coverage (%)
MIN_HIT_LENGTH = 2500       # Minimum hit size (bp)

# Merging rules
MAX_FRAGMENT_GAP = 500      # Max gap between fragments (bp)
MAX_FRAGMENTS_TO_MERGE = 3  # Max fragments to merge
```

### Tuning Recommendations

**If too few Y primes found**:
- Lower `MIN_YPRIME_IDENTITY` to 75%
- Lower `MIN_HIT_LENGTH` to 2000bp
- Increase `MAX_FRAGMENT_GAP` to 1000bp

**If too many fragments (not merging)**:
- Increase `MAX_FRAGMENT_GAP` to 1000bp
- Lower `MIN_HIT_LENGTH` to 2000bp

**If over-merging (grouping unrelated hits)**:
- Increase `MIN_HIT_LENGTH` to 3000bp
- Decrease `MAX_FRAGMENT_GAP` to 300bp
- Decrease `MAX_FRAGMENTS_TO_MERGE` to 2

## Output Files

### 1. TSV File (test_yprimes_7575.tsv)

Example Y prime entries:
```
chr_end  chr              start    end      length  strand  type     source                                pident  bitscore  coverage
chr4R    chr4_extended    1447498  1454103  6606    +       y_prime  MERGED_3_fragments                    97.46   11892.0   148.77
chr2L    chr2_extended    6164     189      5976    -       y_prime  Y_Prime_chr2L1#Short/Solo/ID4_Green   100.00  10778.0   100.00
```

**Columns**:
- `source`: Either original query ID or "MERGED_N_fragments"
- `pident`: Weighted average identity if merged
- `coverage`: May be >100% if multiple overlapping fragments
- `strand`: +/- indicates orientation

### 2. GFF3 File (test_yprimes_7575.gff3)

Example:
```gff3
chr4_extended  TeloTracker  y_prime  1447498  1454103  .  +  .  ID=yprime_19;Name=MERGED_3_fragments;chr_end=chr4R;yprime_class=long;merged=True;fragment_count=3
```

**Attributes**:
- `yprime_class`: short, long, fragment, or unknown
- `merged`: True if created by merging fragments
- `fragment_count`: Number of fragments merged

### 3. BED File (test_yprimes_7575.bed)

Suitable for IGV visualization:
```bed
chr4_extended  1447498  1454103  MERGED_3_fragments  11892  +
```

## Verification Commands

### Count Y primes per chr end
```bash
awk '$7=="y_prime" {print $1}' test_yprimes_7575.tsv | sort | uniq -c
```

### Show Y prime sizes
```bash
awk '$7=="y_prime" {print $1, $5, $8}' test_yprimes_7575.tsv | sort
```

### Find tandem arrays (multiple Y primes per chr end)
```bash
awk '$7=="y_prime" {print $1}' test_yprimes_7575.tsv | sort | uniq -c | awk '$1 > 1'
```

### Show merged vs non-merged
```bash
awk '$7=="y_prime" {merged=($8 ~ /MERGED/); print merged, $1, $5}' test_yprimes_7575.tsv | sort
```

### Check for fragments (too small/large)
```bash
awk '$7=="y_prime" && ($5 < 5000 || $5 > 6900) {print $1, $5}' test_yprimes_7575.tsv
```

## Common Issues and Solutions

### Issue 1: Too many small fragments

**Symptom**: Hundreds of tiny (<1kb) Y prime hits

**Cause**: `MIN_HIT_LENGTH` too low, allowing repetitive matches

**Solution**: Increase `MIN_HIT_LENGTH` to 2500-3000bp

### Issue 2: No Y primes detected

**Symptom**: Zero Y primes found

**Causes**:
- Thresholds too strict
- Reference genome missing telomeric regions
- Significant strain divergence

**Solutions**:
1. Check BLAST output: `wc -l test_yprimes_7575_yprime_blast.txt`
2. Lower identity threshold to 75%
3. Lower hit length to 2000bp
4. Check reference assembly quality

### Issue 3: Y primes on wrong chr ends

**Symptom**: Y prime appears centromere-proximal (wrong side of anchor)

**Cause**: Assignment logic issue or misassembled reference

**Solution**: Manual verification in IGV - may be assembly error

### Issue 4: Over-merging (one huge "Y prime")

**Symptom**: 10-20kb "Y primes" merged from many fragments

**Cause**: `MAX_FRAGMENT_GAP` too large or `MIN_HIT_LENGTH` too small

**Solution**:
- Reduce `MAX_FRAGMENT_GAP` to 300bp
- Increase `MIN_HIT_LENGTH` to 3000bp

## Comparison to Anchor Detection

| Aspect | Anchors | Y Primes |
|--------|---------|----------|
| **Number per chr end** | Exactly 1 | 0 to 10+ |
| **Chromosome constraint** | Strict (synteny enforced) | None (can be anywhere) |
| **Identity threshold** | 95% (high) | 80% (lower, more variable) |
| **Coverage threshold** | 95% (complete) | 50% (allow fragments) |
| **Fragment merging** | None (keep best hit only) | Yes (2-3 fragments) |
| **Size variation** | Uniform (~5kb) | Bimodal (Short: 5-6.4kb, Long: 6.4-6.9kb) |
| **Tandem repeats** | Never | Common (4-10 copies) |
| **Missing** | Rare (assembly issue) | Common (normal variation) |

## Summary

Y prime detection successfully handles:
- ✅ Variable number per chr end (0-10+)
- ✅ Tandem repeat arrays (chr4R: 10 copies)
- ✅ Fragment merging (2-3 large pieces, not hundreds of tiny hits)
- ✅ Size classification (Short vs Long)
- ✅ No chromosome constraint (assigned by proximity to anchor)
- ✅ High quality results (96-100% identity for all detected Y primes)

**Key Innovation**: Filtering for substantial hits (>2.5kb) before merging eliminates 92% of repetitive noise while preserving real Y prime fragments.

## Next Steps

With anchors and Y primes detected, the pipeline can:
1. Infer X prime regions (between anchor and first Y prime)
2. Classify subtelomeric structure for each chr end
3. Compare subtelomeric architecture across strains
4. Identify structural variations in telomeric regions
