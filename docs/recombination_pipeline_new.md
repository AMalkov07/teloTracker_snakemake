# New Recombination Pipeline — How It Works

## Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         INPUTS (from steps 0-6)                        │
│                                                                        │
│  Per-chr-end FASTA files (anchored reads, telomere-confirmed only)     │
│  Day0 reference BED (labeled features per chr end)                     │
│  Y prime library, Spacer library, X element library                    │
└───────────────┬─────────────────────────────────────────────────────────┘
                │
                v
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 9.5: FILTER TELOMERE READS                                      │
│                                                                        │
│  Input:  per-chr-end FASTAs + post_telo_trimming.tsv                   │
│  Filter: only keep reads with repeat_length > 0                        │
│          (confirmed telomeric repeats detected)                        │
│  Result: ~28% of reads removed (partial reads without telomere)        │
│  Output: telomere_filtered_reads/{base}_{chr_end}_telomere_reads.fasta │
└───────────────┬─────────────────────────────────────────────────────────┘
                │
                v
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 10: MINIMAP2 ALIGNMENT (supplementary evidence only)             │
│                                                                        │
│  Purpose: NOT used for classification or gating.                       │
│           Only extracts supplementary alignment contigs.               │
│                                                                        │
│  minimap2 aligns each read to the full day0 reference.                 │
│  If minimap2 splits a read across two contigs (supplementary           │
│  alignment), the second contig name is recorded.                       │
│                                                                        │
│  Example: read from chr9R has supplementary to chr8_extended           │
│           → "chr8" recorded as potential recombination source          │
│                                                                        │
│  Output: *_alignment.tsv  (read_id, supplementary_contigs)             │
│          *.bam            (for manual IGV inspection)                   │
│                                                                        │
│  NOTE: ALL reads proceed to step 11 regardless of alignment results    │
└───────────────┬─────────────────────────────────────────────────────────┘
                │
                v
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 11: COMBINED FEATURE ANALYSIS (per chr end)                      │
│                                                                        │
│  Three analyses run on EVERY read — no gating, no skipping.            │
│  All three run independently, then results are reconciled.             │
│                                                                        │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐                          │
│  │ 11a:     │    │ 11b:     │    │ 11c:     │                          │
│  │ SPACER   │    │ X ELEM   │    │ Y PRIME  │                          │
│  │ CHUNKS   │    │ CHUNKS   │    │ BLAST    │                          │
│  └────┬─────┘    └────┬─────┘    └────┬─────┘                          │
│       │               │               │                                │
│       v               v               v                                │
│  ┌─────────────────────────────────────────────┐                       │
│  │ 11d: CROSS-FEATURE RECONCILIATION           │                       │
│  │      + CONFIDENCE SCORING                   │                       │
│  └─────────────────────────────────────────────┘                       │
│                                                                        │
│  Output: *_features.tsv (one row per read, all analysis combined)      │
└───────────────┬─────────────────────────────────────────────────────────┘
                │
                v
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 12: SUMMARY                                                      │
│                                                                        │
│  Reads all 32 *_features.tsv files                                     │
│  Produces one row per chr end with summary statistics                   │
│                                                                        │
│  Output: *_recombination_summary.tsv                                   │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Step 11a: Spacer Chunk Analysis (detailed)

```
For each read:

0. QUICK CHECK — does the reference spacer exist EXACTLY in this read?

   Uses a literal Python substring check:
     if reference_spacer_sequence in read_sequence:

   No aligner, no BLAST — just exact string matching.
   This checks for a 100% identity, 100% length match.

   If the reference spacer is found as an exact substring:
     → spacer_recombination = "no_change"
     → spacer_quick_check_skipped = True
     → SKIP chunk analysis for this read

   If not (expected for most/all reads due to sequencing errors):
     → Proceed to chunk analysis (steps 1-4 below)

   NOTE: Track what percentage of reads are resolved by this quick check.
   This metric is reported in the output as "spacer_quick_check_skipped"
   and in the summary TSV as "pct_spacer_quick_check_skipped".
   Expected: very few reads will pass (Nanopore ~5% error rate means
   even a perfect spacer will have errors). This check is intentionally
   strict — we want to measure how many reads have a truly perfect match.

1. CHOP the entire read into non-overlapping 250bp chunks

   Read: =====================================================>
   Chunks: [--1--][--2--][--3--][--4--][--5--][--6--]...[--N--]
           250bp  250bp  250bp  250bp  250bp  250bp     250bp

2. BLAST all chunks against the COMBINED spacer library
   (contains 250bp spacer chunks from ALL 32 chr ends)

3. For each chunk, record the BEST matching chr end + identity

4. WALK the chunks in order from start to end:

   NO RECOMBINATION — all chunks match the expected chr end:
   Chunk:    1     2     3     4     5     6     7     8
   Source:  chr9R chr9R chr9R chr9R chr9R chr9R chr9R chr9R
   Result:  spacer_recombination = "no_change"

   RECOMBINATION DETECTED — chunks switch from one source to another:
   Chunk:    1     2     3     4     5     6     7     8
   Source:  chr9R chr9R chr9R chr9R chr8L chr8L chr8L chr8L
                                    ^^^^ SWITCH HERE
   Result:  spacer_recombination = "switch_detected"
            spacer_source = "chr8L"
            spacer_switch_chunk = 4

   FULL SWITCH — all chunks match a DIFFERENT chr end:
   Chunk:    1     2     3     4     5     6     7     8
   Source:  chr8L chr8L chr8L chr8L chr8L chr8L chr8L chr8L
   Result:  spacer_recombination = "full_switch"
            spacer_source = "chr8L"

5. NOTE: Not all chunks will have BLAST hits (chunks in non-spacer
   regions like anchor, X element, Y primes won't match the spacer
   library). Only chunks that fall in a spacer region will produce hits.
   This is a natural filter — no need to know where the spacer is
   in advance.
```

---

## Step 11b: X Element Chunk Analysis (detailed)

```
Same approach as spacer, but BLASTs against the combined X element library.

X elements are shorter than spacers (typically 500-1500bp = 2-6 chunks)
so fewer chunks will have hits. The analysis logic is identical:

Chunk:    ...  12    13    14    15    ...
Source:        chr9R chr9R chr8L chr8L
                          ^^^^ SWITCH
```

---

## Step 11c: Y Prime Analysis (detailed)

```
1. BLAST the FULL READ (not chunked) against Y prime library
   - 80% identity minimum, merge adjacent hits, filter hits < 500bp
   - Deduplicate overlapping hits (keep best score)

2. Order Y prime hits anchor-to-telomere on the read
   (reverse order if telo_side = 'beginning')

3. Compare position-by-position against day0 reference:

   Reference: [ID2] [ID2] [ID2]     (chr4R has 3 Y primes, all ID2)
   Read:      [ID2] [ID1] [ID1]     (position 2 changed from ID2 to ID1)
              match  ^^^^ DIVERGENCE at position 2

   Categories:
   - "No Change"        — all Y primes match reference order
   - "1st Y' Change"    — position 1 (anchor-proximal) already differs
   - "Y' Recombination" — divergence at position 2+
   - "Y' Loss"          — read has fewer Y primes than reference
   - "Y' Gain"          — read has more Y primes than reference

4. Find COMPATIBLE REFERENCE ENDS:
   Which chr ends have a first Y prime matching the read's first Y prime?

   Read's first Y prime = ID1
   Chr ends with ID1: [chr5R, chr8L, chr8R, chr12L, chr16R]
   → These are all potential recombination sources from Y prime evidence alone

   This list is passed to cross-feature reconciliation (step 11d)
```

---

## Step 11d: Cross-Feature Reconciliation + Confidence Scoring

```
After all three analyses complete for a read, compare their conclusions:

CASE 1: ALL FEATURES AGREE
┌─────────────────────────────────────────────────────────────────┐
│  Spacer:    switch to chr8L at chunk 15                        │
│  X element: matches chr8L                                       │
│  Y prime:   compatible with chr8L (ID1)                         │
│  Minimap2:  supplementary alignment to chr8_extended            │
│                                                                 │
│  → recombination_source = "chr8L"                               │
│  → cross_feature_consistent = True                              │
│  → confidence = HIGH (base × 1.1 boost)                        │
└─────────────────────────────────────────────────────────────────┘

CASE 2: FEATURES PARTIALLY AGREE
┌─────────────────────────────────────────────────────────────────┐
│  Spacer:    switch to chr8L at chunk 15                        │
│  X element: no data (too few chunks matched)                    │
│  Y prime:   compatible with chr8L (ID1 matches 5 ends)          │
│                                                                 │
│  → recombination_source = "chr8L"                               │
│  → cross_feature_consistent = True (no contradiction)           │
│  → confidence = MODERATE                                        │
└─────────────────────────────────────────────────────────────────┘

CASE 3: FEATURES DISAGREE
┌─────────────────────────────────────────────────────────────────┐
│  Spacer:    switch to chr8L at chunk 15                        │
│  X element: matches chr12R                                      │
│  Y prime:   compatible with chr8L but NOT chr12R                │
│                                                                 │
│  → recombination_source = "chr8L" (most votes)                  │
│  → cross_feature_consistent = False                             │
│  → is_complex_event = True                                      │
│  → confidence = LOW (base × 0.3 penalty)                       │
│  → Note: may indicate compound recombination                    │
└─────────────────────────────────────────────────────────────────┘

CASE 4: NO RECOMBINATION
┌─────────────────────────────────────────────────────────────────┐
│  Spacer:    all chunks match reference chr end                  │
│  X element: all chunks match reference chr end                  │
│  Y prime:   "No Change"                                         │
│                                                                 │
│  → recombination_detected = False                               │
│  → confidence = 0.95                                            │
└─────────────────────────────────────────────────────────────────┘


CONFIDENCE SCORING — Three Factors:

  Factor 1: FEATURE DISTINCTIVENESS (how unique is this feature type?)
  ┌──────────────┬───────┬─────────────────────────────────────────┐
  │ Feature      │ Weight│ Reason                                  │
  ├──────────────┼───────┼─────────────────────────────────────────┤
  │ Spacer       │  0.9  │ Generally unique between chr ends       │
  │ Y prime      │  0.6  │ Distinct between IDs but ID assignment  │
  │              │       │ can be unreliable                       │
  │ X element    │  0.4  │ Often similar across chr ends           │
  └──────────────┴───────┴─────────────────────────────────────────┘

  Factor 2: SEPARATION (how clear is the best match?)
  - Best match 99%, second-best 50%  → separation = high  → confident
  - Best match 90%, second-best 85%  → separation = low   → ambiguous
  - Formula: (best_identity - second_best_identity) / best_identity

  Factor 3: COMPLEXITY (do features agree?)
  - All agree               → ×1.0 (or ×1.1 with Y prime support)
  - Partially agree         → ×0.7
  - Features disagree       → ×0.3 (large penalty)

  FINAL CONFIDENCE = Factor1 × Factor2 × Factor3
```

---

## Output: *_features.tsv (one row per read)

```
COLUMNS:
─────────────────────────────────────────────────────────────────

Read info:
  read_id, chr_end, read_length, telo_side

Spacer analysis:
  spacer_chunks_analyzed      — how many 250bp chunks had BLAST hits
  spacer_source               — best-matching chr end (or reference if no switch)
  spacer_switch_chunk         — chunk index where source changes (-1 if none)
  spacer_best_identity        — average identity of best-matching source
  spacer_second_best_identity — average identity of second-best source
  spacer_confidence           — per-feature confidence (0-1)
  spacer_recombination        — "no_change" | "switch_detected" | "full_switch" | "no_data"

X element analysis:
  (same columns as spacer, prefixed with x_element_)

Y prime analysis:
  y_prime_count_on_read              — how many Y primes BLAST found
  y_prime_observed_array             — e.g., "ID2,ID1,ID1" (anchor-to-telomere)
  y_prime_recombination_status       — "No Change" | "1st Y' Change" | etc.
  y_prime_divergence_idx             — position of first mismatch (-1 if none)
  y_prime_expected_at_divergence     — reference Y prime ID at that position
  y_prime_found_at_divergence        — actual Y prime ID found
  y_prime_downstream_consistent      — True if all post-divergence Y primes agree
  y_prime_compatible_ends            — chr ends consistent with observed Y prime pattern

Cross-feature reconciliation:
  recombination_detected      — True/False
  recombination_source        — best overall source chr end (or "ambiguous")
  overall_confidence           — 0-1, three-factor score
  cross_feature_consistent    — True if all features agree
  cross_feature_detail        — breakdown of what each feature found
  is_complex_event            — True if features disagree
  confidence_factors          — e.g., "base=0.85, complexity=1.10"

Flags:
  qc_flags                    — any quality warnings
```

---

## Comparison: Old Pipeline vs New Pipeline

```
┌─────────────────────────────┬──────────────────────────────────────────┐
│       OLD PIPELINE          │           NEW PIPELINE                   │
├─────────────────────────────┼──────────────────────────────────────────┤
│ Reads: ALL anchored reads   │ Reads: TELOMERE-CONFIRMED only (~72%)   │
│                             │                                          │
│ Gating: spacer/X analysis   │ Gating: NONE — all three analyses run   │
│ only if first Y prime       │ on every read unconditionally            │
│ is different                │                                          │
│                             │                                          │
│ Alignment: primary driver   │ Alignment: supplementary evidence only   │
│ (4-category classification) │ (no classification, no gating)           │
│                             │                                          │
│ Spacer: RepeatMasker on     │ Spacer: BLAST 250bp chunks of full read │
│ diverged tail only          │ against combined library (all chr ends)  │
│                             │                                          │
│ X element: RepeatMasker on  │ X element: BLAST 250bp chunks of full   │
│ diverged tail only          │ read against combined library            │
│                             │                                          │
│ Y prime: RepeatMasker       │ Y prime: BLAST full read against Y'     │
│ (heavy dependency)          │ library (no RepeatMasker dependency)     │
│                             │                                          │
│ Pairing libraries: per-pair │ Libraries: combined (all chr ends in 1)  │
│ subsets (ID1_and_ID2.fasta) │ No pairing files needed                  │
│                             │                                          │
│ Confidence: binary          │ Confidence: 0-1 three-factor score       │
│ (switch / no switch)        │ (distinctiveness × separation × complex) │
│                             │                                          │
│ Cross-feature: none         │ Cross-feature: consistency check with    │
│ (results placed side by     │ confidence boost/penalty, Y prime        │
│ side, no reconciliation)    │ compatible ends for validation           │
│                             │                                          │
│ Source ID: only from spacer │ Source ID: from spacer + X element +     │
│ (Y prime never identifies   │ Y prime compatible list + minimap2       │
│ donor chr end)              │ supplementary alignment                  │
│                             │                                          │
│ Complex events: not         │ Complex events: flagged with penalty,    │
│ considered                  │ detailed modeling deferred               │
│                             │                                          │
│ Dependencies: minimap2,     │ Dependencies: minimap2, BLAST, samtools  │
│ BLAST, RepeatMasker,        │ (RepeatMasker removed)                   │
│ samtools                    │                                          │
└─────────────────────────────┴──────────────────────────────────────────┘
```
