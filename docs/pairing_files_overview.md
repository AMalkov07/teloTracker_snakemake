# How Pairing Files Are Created and Used

## 0. Y Prime Identification and Color Groups

```
HOW Y PRIMES ARE IDENTIFIED AND GROUPED
========================================

The day0 reference genome contains subtelomeric regions at each chr end.
Y primes are repetitive elements found between the X element and the telomere.
A chr end can have 0 to 7+ Y primes.

  Subtelomeric structure (one chr end):
  [ANCHOR]──[SPACER]──[X element]──[Y' 1]──[Y' 2]──...──[Y' N]──[TELOMERE]
                                     ↑       ↑             ↑
                                     Each Y prime is assigned a color/ID

STEP 1: label_regions.sh runs RepeatMasker on the day0 reference
        against a known Y prime library. This identifies:
        - WHERE each Y prime is (start/end coordinates → BED file)
        - WHICH Y prime it is (sequence similarity → ID assignment)

STEP 2: Each Y prime gets classified into a "color group" (ID1–ID7)
        based on sequence similarity. The naming comes from the library:

        FASTA header format:
        >Y_Prime_chr4R1#Short/Solo/ID2_Green-Dark
                  ↑    ↑     ↑    ↑     ↑
                  |    |     |    |     └── Color name (visual label)
                  |    |     |    └── Y prime ID (ID1–ID7)
                  |    |     └── Array type (Solo/Tandem)
                  |    └── Size class (Short/Long)
                  └── Origin chr end + position number

        The ID (ID1–ID7) is what matters for recombination detection.
        Y primes with the SAME ID are structurally very similar —
        if a read recombines and picks up an ID1 Y prime, you can't
        tell which chr end it came from based on Y prime alone.

STEP 3: The extracted Y primes are written to:
        extracted_yprimes_{strain}.fasta
        (one sequence per Y prime instance across all chr ends)

RESULT: Each chr end gets a Y prime "signature":

  Strain 7302 Y prime assignments:
  ┌─────────┬────────────────────────┬───────────────────────────┐
  │ Chr End  │ Y Prime IDs            │ Notes                     │
  ├─────────┼────────────────────────┼───────────────────────────┤
  │ chr4R   │ [ID2,ID2,ID2,ID2,      │ 7 Y primes, all ID2      │
  │         │  ID2,ID2,ID2]          │                           │
  │ chr12R  │ [ID1,ID2,ID2]          │ 3 Y primes, mixed        │
  │ chr13L  │ [ID2,ID1]             │ 2 Y primes, mixed        │
  │ chr14L  │ [ID3,ID2,ID5]          │ 3 Y primes, all different │
  │ chr5R   │ [ID1]                  │ 1 Y prime                 │
  │ chr9L   │ [ID6]                  │ 1 Y prime                 │
  │ chr1L   │ (none)                 │ No Y primes at this end   │
  │ chr1R   │ (none)                 │ No Y primes at this end   │
  │ ...     │ ...                    │                           │
  └─────────┴────────────────────────┴───────────────────────────┘

  14 chr ends have NO Y primes in the BED → can't do Y prime
  comparison for those ends, rely on spacer/X element analysis.

  Many single-Y-prime chr ends share the same ID:
    ID1: chr5R, chr8L, chr8R, chr12L, chr16R  (5 chr ends)
    ID6: chr5L, chr9L, chr10L, chr14R          (4 chr ends)
    ID2: chr4R, chr15R                          (2 chr ends)

  → If a read's Y prime changes to ID1, it could be from any of
    5 different chr ends. This ambiguity is why spacer/X element
    analysis is needed to narrow down the source.

WHY THIS MATTERS FOR PAIRINGS:
  The Y prime IDs determine which chr ends get grouped together
  in the pairing files. If a read shows ID1 and ID2, only chr ends
  with ID1 or ID2 Y primes are relevant comparison targets.
```

## 1. From Day0 Reference to Pairing Libraries

```
DAY0 REFERENCE GENOME (assembly_{strain}_dorado_reference.fasta)
        |
        v
 label_regions.sh
        |
        +---> pretelomeric_regions_{strain}_simp.bed
        |        (BED file labeling: anchor, spacer, x_core, x_variable, y_prime per chr end)
        |
        +---> extracted_yprimes_{strain}.fasta
        |        (Y prime sequences with color/ID annotations)
        |
        v
 +--------------------------+     +----------------------------+
 | make_chopped_spacer_     |     | make_x_element_            |
 | sequences.py             |     | sequences.py               |
 |                          |     |                            |
 | For each chr end:        |     | For each chr end:          |
 |  1. Read BED spacer      |     |  1. Read BED x_core +     |
 |     coordinates          |     |     x_variable coordinates |
 |  2. Extract spacer       |     |  2. Extract X element      |
 |     from reference       |     |     from reference         |
 |  3. Chop into 250bp      |     |  3. Keep as whole sequence |
 |     non-overlapping      |     |                            |
 |     chunks               |     |                            |
 +--------------------------+     +----------------------------+
        |                                    |
        v                                    v
 {strain}_50kb_chopped_up_          {strain}_whole_x_regions_
 spacer_sequences.fasta             sequences.fasta
        |                                    |
        v                                    v
 +--------------------------+     +----------------------------+
 | make_databases_of_       |     | make_databases_of_         |
 | pairings_for_spacers.py  |     | pairings_for_x_elements.py |
 |                          |     |                            |
 | Each pairing file is a   |     | Each pairing file is a     |
 | SUBSET of the full 50kb  |     | SUBSET of the full X       |
 | spacer file above.       |     | regions file above.        |
 |                          |     |                            |
 | For each Y prime ID pair |     | For each Y prime ID pair   |
 | (e.g. ID1+ID2):          |     | (e.g. ID1+ID2):            |
 |  1. Find all chr ends    |     |  1. Find all chr ends      |
 |     that have ID1 or ID2 |     |     that have ID1 or ID2   |
 |     in ANY Y prime       |     |     in ANY Y prime         |
 |     position             |     |     position               |
 |  2. Extract just those   |     |  2. Extract just those     |
 |     chr ends' chunks     |     |     chr ends' sequences    |
 |     from the full file   |     |     from the full file     |
 |  3. Write to one FASTA   |     |  3. Write to one FASTA     |
 +--------------------------+     +----------------------------+
        |                                    |
        v                                    v
 pairings_for_spacers/              pairings_for_x_element_ends/
   {strain}_pairings/                 {strain}_pairings/
     7302_paired_ID1_and_ID2.fasta     7302_paired_ID1_and_ID2.fasta
     7302_paired_ID1_and_ID3.fasta     7302_paired_ID1_and_ID3.fasta
     7302_paired_chr4R_and_ID1.fasta   7302_paired_chr4R_and_ID1.fasta
     ...                               ...

 Each pairing file is a SUBSET of the full spacer/X element file,
 containing only the entries from chr ends that have the relevant
 Y prime IDs (in any position of their Y prime array).

 NOTE: Two types of pairing files are created:
   - ID_and_ID files (e.g., 7302_paired_ID1_and_ID2.fasta)
   - chr_and_ID files (e.g., 7302_paired_chr4R_and_ID1.fasta)
 However, only the chr_and_ID files are actually used by the
 current TeloTracker analysis pipeline. The ID_and_ID files are
 created as a byproduct but never loaded, because the read
 grouping step (make_pairings_from_y_primes_all_ends.py) only
 produces {chr_end}_and_{ID} read groups.
```

## 2. Why Pairings Exist — The Grouping Logic

```
EXAMPLE: Strain 7302, Y prime assignments (from extracted_yprimes)

  chr4R  → ID2       chr8L  → ID1       chr12L → ID1
  chr5R  → ID1       chr8R  → ID1       chr14L → ID3,ID2,ID5
  chr5L  → ID6       chr9L  → ID6       chr15R → ID2
  chr6L  → ID4       chr10L → ID6       chr16R → ID1
  chr7R  → ID5       chr12R → ID1,ID2,ID2  ...

                    ┌──────────────────────────────────────┐
                    │  Pairing: "ID1_and_ID2"              │
                    │                                      │
                    │  Which chr ends have ID1 or ID2?     │
                    │    ID1: chr5R, chr8L, chr8R, chr12L, │
                    │         chr12R, chr16R               │
                    │    ID2: chr4R, chr12R, chr14L, chr15R│
                    │                                      │
                    │  Combined spacer library contains:   │
                    │    chr4R  spacer chunks (250bp x N)  │
                    │    chr5R  spacer chunks (250bp x N)  │
                    │    chr8L  spacer chunks (250bp x N)  │
                    │    chr8R  spacer chunks (250bp x N)  │
                    │    chr12L spacer chunks (250bp x N)  │
                    │    chr12R spacer chunks (250bp x N)  │
                    │    chr14L spacer chunks (250bp x N)  │
                    │    chr15R spacer chunks (250bp x N)  │
                    │    chr16R spacer chunks (250bp x N)  │
                    └──────────────────────────────────────┘
```

## 3. How Pairings Are Used in Analysis

```
                    ANALYSIS PIPELINE
                    =================

Step 8 (RepeatMasker Y primes):
  For each chr end, RepeatMasker finds Y prime hits on reads.
  If a chr4R read shows Y prime ID1 (not its expected ID2):
    → This read is a recombination candidate
    → It goes into the "ID1_and_ID2" read group

         ┌─────────────────────────────────────────────────┐
         │  make_pairings_from_y_primes_all_ends.py        │
         │                                                 │
         │  *** ONLY LOOKS AT THE FIRST (ANCHOR-PROXIMAL)  │
         │  Y PRIME ON EACH READ ***                       │
         │                                                 │
         │  For each read:                                 │
         │    Get the read's 1st Y prime ID                │
         │    Compare to reference 1st Y prime             │
         │                                                 │
         │  1st Y prime MATCHES reference:                 │
         │    → Read is DROPPED entirely.                  │
         │    → NO spacer or X element analysis will run.  │
         │    → Even if positions 2,3,etc differ, the read │
         │      is excluded from all downstream analysis.  │
         │                                                 │
         │  1st Y prime DIFFERS from reference:            │
         │    → Read goes into group file:                 │
         │    paired_by_y_prime_reads/chr4R_and_ID1.fasta  │
         │    ({anchor_chr_end}_and_{unexpected_1st_ID})   │
         │                                                 │
         │  *** BLIND SPOT: Recombination that only affects│
         │  Y primes at positions 2+ is silently missed by │
         │  all spacer/X element follow-up analysis ***    │
         └─────────────────────────────────────────────────┘
                           |
                           v
Step 9-10 (RepeatMasker X elements + Spacers):
  ONLY runs on reads where 1st Y prime differed (from step above).
  For each read group (e.g., paired_by_y_prime_reads/chr4R_and_ID1.fasta):
    Run RepeatMasker against the MATCHING reference pairing library
    (7302_paired_chr4R_and_ID1.fasta from pairings_for_spacers/)

         ┌─────────────────────────────────────────────────┐
         │                                                 │
         │  Reads: paired_by_y_prime_reads/                │
         │    chr4R_and_ID1.fasta                          │
         │         │                                       │
         │         ├──→ RepeatMasker vs spacer pairing     │
         │         │    7302_paired_chr4R_and_ID1.fasta    │
         │         │    (spacer 250bp chunks from chr4R    │
         │         │     + all chr ends with ID1)          │
         │         │         │                             │
         │         │         v                             │
         │         │    Which chr end's chunks match?      │
         │         │    chr4R chunks → came from chr4R     │
         │         │    chr8L chunks → came from chr8L     │
         │         │                                       │
         │         └──→ RepeatMasker vs X element pairing  │
         │              7302_paired_chr4R_and_ID1.fasta    │
         │              (X elements from chr4R + all       │
         │               chr ends with ID1)               │
         │                   │                             │
         │                   v                             │
         │              Which chr end's X element matches? │
         └─────────────────────────────────────────────────┘
                           |
                           v
         ┌─────────────────────────────────────────────────┐
         │  get_recombination_switch_location.py           │
         │                                                 │
         │  For each read's spacer chunks:                 │
         │    Walk from anchor outward                     │
         │    Chunk 1: matches chr4R  ← expected (anchor)  │
         │    Chunk 2: matches chr4R  ← still same         │
         │    Chunk 3: matches chr8L  ← SWITCH!            │
         │    Chunk 4: matches chr8L  ← consistent         │
         │                                                 │
         │  → "Switch from chr4R to chr8L at chunk 3"      │
         │  → Switch distance = midpoint of chunk 3 coords │
         └─────────────────────────────────────────────────┘
```

## 4. Recombination Analysis — How It Works

### 4a. TeloTracker (GitHub) Approach — RepeatMasker-based

```
TELOTRACKER RECOMBINATION PIPELINE (Steps 8-10 of telomere_analysis.sh)
=======================================================================

Input: Per-chr-end reads from Step 3 (e.g., chr4R_anchor_reads.fasta)
       Day0 reference libraries (Y prime, spacer pairings, X element pairings)

┌─────────────────────────────────────────────────────────────────────────┐
│ STEP 8: Y PRIME ANALYSIS (RepeatMasker)                                │
│                                                                        │
│ For EACH of the 32 chr ends:                                           │
│                                                                        │
│   chr4R_anchor_reads.fasta                                             │
│         │                                                              │
│         v                                                              │
│   RepeatMasker -lib repeatmasker_{strain}_all_y_primes.fasta           │
│         │                                                              │
│         v                                                              │
│   .out file: which Y prime IDs are on each read, and where             │
│         │                                                              │
│         v                                                              │
│   ┌──────────────────────────────────────────────────────────┐         │
│   │ make_y_prime_repeatmasker_tsv.py                         │         │
│   │  - Parse .out file into clean TSV                        │         │
│   │  - For telo_side='beginning', reverse Y prime order      │         │
│   │    so index 0 = anchor-proximal                          │         │
│   └──────────────────────────────────────────────────────────┘         │
│         │                                                              │
│         v                                                              │
│   ┌──────────────────────────────────────────────────────────┐         │
│   │ get_stats_of_recombination.py                            │         │
│   │                                                          │         │
│   │  For each read, compare Y prime array to reference:      │         │
│   │                                                          │         │
│   │  Reference (chr4R): [ID2, ID2, ID2, ID2, ID2, ID2, ID2] │         │
│   │  Read found:        [ID2, ID2, ID1, ID1, ID1]           │         │
│   │                      pos1  pos2  pos3  pos4  pos5        │         │
│   │                             ↑                            │         │
│   │                      Position 3: expected ID2, found ID1 │         │
│   │                      → "Y' Recombination" at position 3  │         │
│   │                                                          │         │
│   │  Categories:                                             │         │
│   │    "No Change"        — all positions match reference     │         │
│   │    "1st Y' Change"    — position 1 already differs        │         │
│   │    "Y' Recombination" — divergence at position 2+         │         │
│   │    "Y' Loss"          — fewer Y primes than reference     │         │
│   │                                                          │         │
│   │  Also tracks: count changes (gain/loss), delta groups    │         │
│   └──────────────────────────────────────────────────────────┘         │
│         │                                                              │
│         v                                                              │
│   ┌──────────────────────────────────────────────────────────┐         │
│   │ make_pairings_from_y_primes_all_ends.py                  │         │
│   │                                                          │         │
│   │  *** ONLY CHECKS THE 1st (ANCHOR-PROXIMAL) Y PRIME ***  │         │
│   │                                                          │         │
│   │  For each read:                                          │         │
│   │    Compare read's 1st Y prime ID to reference            │         │
│   │                                                          │         │
│   │  If 1st Y prime MATCHES reference:                       │         │
│   │    → Read is DROPPED. Spacer and X element analysis      │         │
│   │      will NEVER run on this read, even if Y primes       │         │
│   │      at positions 2, 3, etc. are different.              │         │
│   │                                                          │         │
│   │  If 1st Y prime DIFFERS:                                 │         │
│   │    → chr4R read with 1st Y prime = ID1 (expected ID2)    │         │
│   │    → paired_by_y_prime_reads/chr4R_and_ID1.fasta         │         │
│   │                                                          │         │
│   │  Only these grouped reads get spacer/X element analysis  │         │
│   └──────────────────────────────────────────────────────────┘         │
└─────────────────────────────────────────────────────────────────────────┘
         │
         v
┌─────────────────────────────────────────────────────────────────────────┐
│ STEP 9: X ELEMENT ANALYSIS (RepeatMasker)                              │
│                                                                        │
│ For EACH read group (not per chr end — per Y prime pairing):           │
│                                                                        │
│   paired_by_y_prime_reads/chr4R_and_ID1.fasta                          │
│         │                                                              │
│         v                                                              │
│         v                                                              │
|   RepeatMasker -lib {strain}_paired_chr4R_and_ID1.fasta <reads.fasta>  |    
│         │                                                              │
│         v                                                              │
│   .out file: which chr end's X element matches each read               │
│         │                                                              │
│         v                                                              │
│   ┌──────────────────────────────────────────────────────────┐         │
│   │ make_x_element_ends_pairs_repeatmasker_tsv.py            │         │
│   │                                                          │         │
│   │  For each read:                                          │         │
│   │    If X element matches anchor chr end → "No Switch"     │         │
│   │    If X element matches different chr end → "Switch"     │         │
│   │                                                          │         │
│   │  Filters: SW_score ≥ 500, divergence ≤ 2%               │         │
│   └──────────────────────────────────────────────────────────┘         │
└─────────────────────────────────────────────────────────────────────────┘
         │
         v
┌─────────────────────────────────────────────────────────────────────────┐
│ STEP 10: SPACER ANALYSIS (RepeatMasker) + SWITCH LOCATION              │
│                                                                        │
│ For EACH read group:                                                   │
│                                                                        │
│   paired_by_y_prime_reads/chr4R_and_ID1.fasta                          │
│         │                                                              │
│         v                                                              │
|   RepeatMasker -lib {strain}_paired_chr4R_and_ID1.fasta <reads.fasta>  |    
│         │                                                              │
│         v                                                              │
│   .out file: which chr end's spacer chunks match each read position    │
│         │                                                              │
│         v                                                              │
│   ┌──────────────────────────────────────────────────────────┐         │
│   │ make_spacer_pairs_repeatmasker_tsv.py                    │         │
│   │  - Parse RepeatMasker results into TSV                   │         │
│   └──────────────────────────────────────────────────────────┘         │
│         │                                                              │
│         v                                                              │
│   ┌──────────────────────────────────────────────────────────┐         │
│   │ get_recombination_switch_location.py                     │         │
│   │                                                          │         │
│   │  For each read, walk through spacer chunk matches:       │         │
│   │                                                          │         │
│   │  Read anchored to chr4R:                                 │         │
│   │    Chunk  1: matches chr4R_spacer_section_1  ← same      │         │
│   │    Chunk  2: matches chr4R_spacer_section_2  ← same      │         │
│   │    Chunk  3: matches chr4R_spacer_section_3  ← same      │         │
│   │    Chunk  4: matches chr8L_spacer_section_4  ← SWITCH!   │         │
│   │    Chunk  5: matches chr8L_spacer_section_5  ← consistent│         │
│   │    Chunk  6: matches chr8L_spacer_section_6  ← consistent│         │
│   │                                                          │         │
│   │  Outputs:                                                │         │
│   │    Switch pair: "chr4R → chr8L"                          │         │
│   │    Conservative distance: bp from anchor to first switch │         │
│   │    Aggressive distance: bp to furthest switch            │         │
│   │    X element switch: yes/no                              │         │
│   │    ITS in donor: yes/no                                  │         │
│   └──────────────────────────────────────────────────────────┘         │
└─────────────────────────────────────────────────────────────────────────┘
```

### 4b. Snakemake Pipeline Approach — Alignment + BLAST-based

```
SNAKEMAKE RECOMBINATION PIPELINE (Steps 10-13)
===============================================

Input: Per-chr-end reads from Step 3 (e.g., chr4R_anchor_reads.fasta)
       Day0 reference genome + BED + Y prime/spacer/X element libraries

┌─────────────────────────────────────────────────────────────────────────┐
│ STEP 10: ALIGNMENT & CLASSIFICATION (minimap2)                         │
│ Script: run_alignment.py                                               │
│                                                                        │
│ For EACH of the 32 chr ends:                                           │
│                                                                        │
│   chr4R_anchor_reads.fasta                                             │
│         │                                                              │
│         v                                                              │
│   minimap2 -x map-ont → align full reads to day0 reference genome      │
│         │                                                              │
│         v                                                              │
│   Parse BAM: examine primary + supplementary alignments                │
│         │                                                              │
│         v                                                              │
│   CLASSIFY each read into one of 4 categories:                         │
│                                                                        │
│   ┌─────────────────────────────────────────────────────────────┐      │
│   │                                                             │      │
│   │  no_recombination:  Entire read aligns to expected chr end  │      │
│   │                     Identity ≥ 95%, Coverage ≥ 95%          │      │
│   │                     No supplementary alignments             │      │
│   │                                                             │      │
│   │  two_clean_halves:  Read has 2+ alignment segments that     │      │
│   │                     chain cleanly across the read           │      │
│   │                     Segment 1 (anchor-proximal): chr9       │      │
│   │                     Segment 2 (telomere-proximal): chr14    │      │
│   │                     → Breakpoint between segments           │      │
│   │                     → Source contig directly identified      │      │
│   │                                                             │      │
│   │  one_clean_half:    Good alignment from anchor side, but    │      │
│   │                     telomere side unaligned or low quality   │      │
│   │                     → Breakpoint at end of good alignment   │      │
│   │                     → Diverged tail extracted for BLAST      │      │
│   │                     → Source unknown (determined in Step 11) │      │
│   │                                                             │      │
│   │  no_clean_halves:   No good alignment to expected chr end   │      │
│   │                     → Suspicious, flagged as QC issue       │      │
│   │                                                             │      │
│   └─────────────────────────────────────────────────────────────┘      │
│                                                                        │
│   For each breakpoint, record:                                         │
│     - Position on read and reference                                   │
│     - Which BED feature it falls in (spacer? X element? Y prime?)     │
│     - Whether it's mid-element or at a boundary                        │
│     - The source contig (if two_clean_halves)                         │
│                                                                        │
│   Output: *_breakpoints.tsv, *.bam, *_diverged_tails.fasta            │
└─────────────────────────────────────────────────────────────────────────┘
         │
         v
┌─────────────────────────────────────────────────────────────────────────┐
│ STEP 11: COMBINED FEATURE ANALYSIS (BLAST)                             │
│ Script: analyze_features.py                                            │
│                                                                        │
│ For EACH of the 32 chr ends:                                           │
│                                                                        │
│   BATCH PREPROCESSING (run once per chr end):                          │
│   ┌─────────────────────────────────────────────────────────────┐      │
│   │ 1. BLAST full reads vs Y prime library                      │      │
│   │    (replaces RepeatMasker — finds Y prime positions + IDs)  │      │
│   │    - Merge adjacent hits from same Y prime (≤500bp gap)     │      │
│   │    - Filter hits < 500bp                                    │      │
│   │    - Deduplicate overlapping hits (keep best score)         │      │
│   │                                                             │      │
│   │ 2. BLAST diverged tails vs spacer pairing library           │      │
│   │    (identifies spacer source for reads with spacer breaks)  │      │
│   │                                                             │      │
│   │ 3. BLAST diverged tails vs X element pairing library        │      │
│   │    (identifies X element source for reads with X breaks)    │      │
│   └─────────────────────────────────────────────────────────────┘      │
│         │                                                              │
│         v                                                              │
│   PER-READ STRUCTURAL WALK (anchor → telomere):                        │
│   ┌─────────────────────────────────────────────────────────────┐      │
│   │                                                             │      │
│   │ For each read, use alignment coordinates + BED features:    │      │
│   │                                                             │      │
│   │ 1. FEATURE REACH — which features does the read cover?      │      │
│   │    Read aligned to ref positions 420,000 → 437,000          │      │
│   │    BED says:                                                │      │
│   │      420,920-425,960 = anchor        ✓ covered              │      │
│   │      425,961-435,915 = spacer        ✓ covered              │      │
│   │      435,916-436,380 = x_core        ✓ covered              │      │
│   │      436,380-436,661 = x_variable    ✓ covered              │      │
│   │    → feature_reach = "anchor,spacer,x_core,x_variable"     │      │
│   │                                                             │      │
│   │ 2. SPACER CHECK — if alignment covers spacer region:        │      │
│   │    spacer_found = True                                      │      │
│   │    If breakpoint is IN spacer → BLAST diverged tail         │      │
│   │    against spacer library → identify source chr end         │      │
│   │    Also: 250bp chunk analysis around breakpoint             │      │
│   │                                                             │      │
│   │ 3. X ELEMENT CHECK — if alignment covers X element:         │      │
│   │    x_element_found = True                                   │      │
│   │    If breakpoint is IN X element → BLAST diverged tail      │      │
│   │    against X element library → identify source chr end      │      │
│   │                                                             │      │
│   │ 4. Y PRIME ANALYSIS — using BLAST hits from preprocessing:  │      │
│   │    Position-by-position comparison vs reference Y prime      │      │
│   │    order (from BED, not hardcoded):                         │      │
│   │                                                             │      │
│   │    Reference: [ID2, ID2, ID2, ID2, ID2, ID2, ID2]          │      │
│   │    Read BLAST: [ID2, ID2, ID1, ID1, ID1]                   │      │
│   │                          ↑                                  │      │
│   │    Divergence at position 3 → recombination                │      │
│   │                                                             │      │
│   │    Color-group ambiguity: if found ID is shared by          │      │
│   │    multiple chr ends → confidence split equally             │      │
│   │                                                             │      │
│   │ 5. POST-Y-PRIME TELOMERE CHECK:                             │      │
│   │    After last Y prime, is the remaining sequence            │      │
│   │    telomeric? If not → flagged (structural issue)           │      │
│   │                                                             │      │
│   │ 6. UNIFIED HYPOTHESIS with confidence scoring:              │      │
│   │    Structural evidence (from step 10): weight 0.50          │      │
│   │    Y prime evidence: weight 0.35                            │      │
│   │    Spacer/X element evidence: weight 0.15                   │      │
│   │    Cross-element consistency boost: ×1.3                    │      │
│   │    Output: up to 5 ranked hypotheses per read              │      │
│   │                                                             │      │
│   └─────────────────────────────────────────────────────────────┘      │
│                                                                        │
│   Output: *_features.tsv (one row per read, all evidence combined)     │
└─────────────────────────────────────────────────────────────────────────┘
         │
         v
┌─────────────────────────────────────────────────────────────────────────┐
│ STEP 12: CROSS-CHR-END SUMMARY                                         │
│ Script: aggregate_recombination.py --summarize                         │
│                                                                        │
│   Reads all 32 *_features.tsv files                                    │
│   Produces one summary row per chr end:                                │
│     - Total reads                                                      │
│     - % no_recombination / two_clean_halves / one_clean_half           │
│     - Most common recombination source                                 │
│     - Mean confidence of top hypothesis                                │
│                                                                        │
│   Output: *_recombination_summary.tsv                                  │
└─────────────────────────────────────────────────────────────────────────┘
```

### 4c. Key Differences Between the Two Approaches

```
┌───────────────────────┬──────────────────────┬──────────────────────────┐
│ Aspect                │ TeloTracker (GitHub)  │ Snakemake Pipeline       │
├───────────────────────┼──────────────────────┼──────────────────────────┤
│                       │                      │                          │
│ Breakpoint detection  │ IMPLICIT             │ EXPLICIT                 │
│                       │ Inferred from where  │ minimap2 supplementary   │
│                       │ RepeatMasker matches │ alignments show exactly  │
│                       │ change source        │ where the read splits    │
│                       │                      │                          │
│ Primary tool          │ RepeatMasker         │ minimap2 + BLAST         │
│                       │ (for all 3 features) │ (minimap2 for breakpoint │
│                       │                      │  BLAST for features)     │
│                       │                      │                          │
│ Analysis unit         │ Read GROUPS          │ Individual READS         │
│                       │ (grouped by Y prime  │ (each read analyzed      │
│                       │  pairing first)      │  independently per       │
│                       │                      │  chr end)                │
│                       │                      │                          │
│ Y prime reference     │ HARDCODED dicts      │ DYNAMIC from BED file    │
│                       │ per strain in Python │ (works for any strain    │
│                       │                      │  without code changes)   │
│                       │                      │                          │
│ Spacer analysis       │ RepeatMasker matches │ BLAST of diverged tail + │
│                       │ against 250bp chunks │ 250bp chunk analysis     │
│                       │ in pairing library   │ around breakpoint        │
│                       │                      │                          │
│ Processing order      │ Y prime FIRST →      │ All features analyzed    │
│                       │ groups reads →       │ TOGETHER per read in     │
│                       │ THEN spacer/X elem   │ one structural walk      │
│                       │                      │                          │
│ Output format         │ Categorical          │ Quantitative             │
│                       │ "No Change"          │ 0-1 confidence scores    │
│                       │ "Switch to chr8L"    │ Multiple ranked          │
│                       │                      │ hypotheses               │
│                       │                      │                          │
│ Ambiguity handling    │ Reports the match    │ Distributes confidence   │
│                       │ as-is (no scoring)   │ across same-color-group  │
│                       │                      │ sources equally          │
│                       │                      │                          │
│ Dependencies          │ RepeatMasker         │ minimap2, BLAST,         │
│                       │                      │ samtools, pysam          │
│                       │                      │                          │
│ Parallelism           │ Per Y prime pairing  │ Per chr end (32 jobs)    │
│                       │ group                │                          │
│                       │                      │                          │
└───────────────────────┴──────────────────────┴──────────────────────────┘
```

## 5. Differences: GitHub vs Local Pairing Files

```
┌───────────────────────────────────┬───────────────────────────────────┐
│  LOCAL (label_regions.sh output)  │  GITHUB (TeloTracker shipped)     │
├───────────────────────────────────┼───────────────────────────────────┤
│                                   │                                   │
│  SPACER PAIRINGS                  │  SPACER PAIRINGS                  │
│  119 files                        │  245 files                        │
│  21 ID_and_ID + 98 chr_and_ID     │  21 ID_and_ID + 224 chr_and_ID   │
│                                   │                                   │
│  Format:                          │  Format:                          │
│  >chr4R_50kb_spacer_section_1_    │  >chr4R_50kb_spacer_section_1_    │
│   from_repeat_to_plus_50kb#1      │   from_repeat_to_plus_50kb#1      │
│  GCCCCTATAACGTCGATCTT...          │  TTTAGGATATTGCGGTTAGC...          │
│  (250bp chunks) ✓                 │  (250bp chunks) ✓                 │
│                                   │                                   │
│  Same header format ✓             │  Same header format ✓             │
│  Same chr ends per file ✓         │  Same chr ends per file ✓         │
│  DIFFERENT sequences ✗            │  DIFFERENT sequences ✗            │
│  (different reference genome)     │  (different reference genome)     │
│                                   │                                   │
├───────────────────────────────────┼───────────────────────────────────┤
│                                   │                                   │
│  X ELEMENT PAIRINGS               │  X ELEMENT PAIRINGS               │
│  119 files                        │  245 files                        │
│                                   │                                   │
│  Format:                          │  Format:                          │
│  >chr12L_x_ends                   │  >chr4R_x_ends_region             │
│  CAATATGTTTATTTCGTAAA...          │  GGTAGAACAATAGTATGGTG...          │
│  (VARIABLE length: 375-754bp)     │  (FIXED 250bp chunks)             │
│                                   │                                   │
│  Different header suffix ✗        │  Different header suffix ✗        │
│  Whole X regions (not chunked)    │  Chunked into 250bp pieces        │
│                                   │                                   │
├───────────────────────────────────┼───────────────────────────────────┤
│                                   │                                   │
│  MISSING chr_and_ID pairings:     │  EXTRA chr_and_ID pairings:       │
│  (only 98 of possible 224)        │  chr10L_and_IDx                   │
│  Skips chr ends without Y primes  │  chr12L_and_IDx                   │
│  in the BED file                  │  chr13L_and_IDx  ... etc          │
│                                   │  (pairs ALL chr ends, even those  │
│                                   │   without Y primes in BED)        │
│                                   │                                   │
└───────────────────────────────────┴───────────────────────────────────┘

KEY DIFFERENCES SUMMARY:
  1. File count:    GitHub has ~2x more files (includes all chr_end pairings)
  2. Sequences:     Different (built from different reference genome assemblies)
  3. X elements:    GitHub chunks to 250bp; local keeps whole regions
  4. Header naming: X elements differ ("x_ends_region" vs "x_ends")
  5. Coverage:      GitHub pairs every chr end; local skips those without Y primes
```
