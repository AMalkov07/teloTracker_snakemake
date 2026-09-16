# Per-sample chr-end mismatch counts

For every sample: how many times each chromosome end had a read whose Y' did not match
its expected group (recipient), and how many times that end's Y' element was the one
matched instead (donor / "mistaken-for"). Counts include ALL evidence categories --
strong, weak, FAILS, no junction, and reference defect -- so this is raw mismatch volume,
not confirmed recombination. See the per-sample annotated_summary.tsv for the evidence
breakdown of any individual row.

## 6991_day0  (4 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr10L | 1 | 0 |
| chr13L | 1 | 1 |
| chr14L | 1 | 0 |
| chr2L | 0 | 1 |
| chr4R | 0 | 1 |
| chr6L | 1 | 0 |
| chr7R | 0 | 1 |

## 6991_day0_TeloTag  (6 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr13L | 1 | 2 |
| chr14L | 0 | 2 |
| chr14R | 3 | 0 |
| chr16R | 1 | 0 |
| chr5R | 0 | 1 |
| chr7R | 1 | 0 |
| chr8R | 0 | 1 |

## 6991_day0_TeloTag_with_selection  (11 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr12R | 0 | 4 |
| chr13L | 1 | 0 |
| chr14L | 1 | 2 |
| chr14R | 0 | 2 |
| chr15R | 0 | 1 |
| chr16L | 1 | 1 |
| chr16R | 2 | 0 |
| chr2L | 1 | 1 |
| chr5L | 2 | 0 |
| chr7R | 1 | 0 |
| chr8L | 1 | 0 |
| chr8R | 1 | 0 |

## 6991_day0_reference  (27 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr10L | 0 | 2 |
| chr12L | 3 | 0 |
| chr12R | 0 | 1 |
| chr13L | 4 | 3 |
| chr14L | 0 | 5 |
| chr14R | 5 | 0 |
| chr15R | 1 | 1 |
| chr16L | 2 | 0 |
| chr16R | 1 | 1 |
| chr2L | 3 | 5 |
| chr5L | 2 | 2 |
| chr5R | 0 | 1 |
| chr6L | 4 | 0 |
| chr8L | 1 | 0 |
| chr8R | 1 | 2 |
| chr9L | 0 | 4 |

## 6991_day0_reference_promethion  (12 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr13L | 4 | 3 |
| chr14L | 0 | 1 |
| chr14R | 2 | 0 |
| chr16R | 2 | 0 |
| chr2L | 1 | 3 |
| chr5R | 1 | 0 |
| chr6L | 0 | 2 |
| chr8L | 1 | 0 |
| chr8R | 1 | 1 |
| chr9L | 0 | 2 |

## 6991_day0_with_selection  (439 total mismatched copies) **388 of chr14L's 391 "as recipient" and chr7R's 391 "as donor" are the known chr14L-1 assembly defect (5,720 bp vs 6,654), not real mismatches.**

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr10L | 5 | 5 |
| chr12L | 2 | 0 |
| chr13L | 18 | 2 |
| chr14L | 388 | 5 |
| chr14R | 8 | 0 |
| chr15R | 0 | 1 |
| chr16L | 2 | 6 |
| chr16R | 2 | 0 |
| chr2L | 3 | 12 |
| chr5R | 3 | 2 |
| chr6L | 3 | 4 |
| chr7R | 1 | 391 |
| chr8L | 2 | 1 |
| chr8R | 2 | 5 |
| chr9L | 0 | 5 |

## 6991_day0_with_selection_repeat  (18 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr10L | 3 | 1 |
| chr13L | 9 | 1 |
| chr14L | 0 | 3 |
| chr14R | 3 | 0 |
| chr16L | 0 | 2 |
| chr16R | 1 | 0 |
| chr2L | 1 | 2 |
| chr5L | 1 | 0 |
| chr6L | 0 | 6 |
| chr8L | 0 | 1 |
| chr8R | 0 | 1 |
| chr9L | 0 | 1 |

## 6991_day0_with_selection_repeat2  (13 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr10L | 2 | 0 |
| chr12L | 2 | 0 |
| chr13L | 3 | 2 |
| chr14L | 0 | 2 |
| chr14R | 1 | 0 |
| chr15R | 1 | 0 |
| chr16R | 2 | 1 |
| chr2L | 2 | 0 |
| chr6L | 0 | 3 |
| chr7R | 0 | 2 |
| chr8R | 0 | 2 |
| chr9L | 0 | 1 |

## 7172_day0_with_selection  (4 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr13L | 0 | 1 |
| chr14L | 0 | 3 |
| chr16R | 2 | 0 |
| chr5R | 1 | 0 |
| chr6L | 1 | 0 |

## 7302_day0_with_selection  (14 total mismatched copies)

| end | as recipient (mismatched) | as donor (mistaken-for) |
|---|---|---|
| chr10L | 1 | 1 |
| chr14L | 0 | 6 |
| chr14R | 2 | 3 |
| chr16L | 0 | 1 |
| chr2L | 1 | 1 |
| chr5R | 2 | 0 |
| chr6L | 7 | 0 |
| chr8R | 1 | 2 |
# Combined (all 10 samples)

548 total mismatched copies pooled across all ten samples (160 excluding
the 388 chr14L-1 reference-defect rows from 6991_day0_with_selection).

| end | recipient (all) | recipient (excl. defect) | donor (all) | donor (excl. defect) |
|---|---|---|---|---|
| chr10L | 12 | 12 | 9 | 9 |
| chr12L | 7 | 7 | 0 | 0 |
| chr12R | 0 | 0 | 5 | 5 |
| chr13L | 41 | 41 | 15 | 15 |
| chr14L | 390 | 2 | 29 | 29 |
| chr14R | 24 | 24 | 5 | 5 |
| chr15R | 2 | 2 | 3 | 3 |
| chr16L | 5 | 5 | 10 | 7 |
| chr16R | 13 | 13 | 2 | 2 |
| chr2L | 12 | 12 | 25 | 25 |
| chr4R | 0 | 0 | 1 | 1 |
| chr5L | 5 | 5 | 2 | 2 |
| chr5R | 7 | 7 | 4 | 4 |
| chr6L | 16 | 16 | 15 | 15 |
| chr7R | 3 | 3 | 394 | 9 |
| chr8L | 5 | 5 | 2 | 2 |
| chr8R | 6 | 6 | 14 | 14 |
| chr9L | 0 | 0 | 13 | 13 |

| **TOTAL** | **548** | **160** | **548** | **160** |# 6991 combined, excluding 6991_day0_with_selection

Pooled across the seven 6991 day-0 samples other than 6991_day0_with_selection, whose
chr14L-1 is mis-assembled (5,720 bp vs 6,654 bp) and produces 388 artefactual rows.
Samples included: 6991_day0, 6991_day0_TeloTag, 6991_day0_TeloTag_with_selection, 6991_day0_reference, 6991_day0_reference_promethion, 6991_day0_with_selection_repeat, 6991_day0_with_selection_repeat2

91 total mismatched copies.

| end | recipient (mismatched) | donor (mistaken-for) |
|---|---|---|
| chr10L | 6 | 3 |
| chr12L | 5 | 0 |
| chr12R | 0 | 5 |
| chr13L | 23 | 12 |
| chr14L | 2 | 15 |
| chr14R | 14 | 2 |
| chr15R | 2 | 2 |
| chr16L | 3 | 3 |
| chr16R | 9 | 2 |
| chr2L | 8 | 12 |
| chr4R | 0 | 1 |
| chr5L | 5 | 2 |
| chr5R | 1 | 2 |
| chr6L | 5 | 11 |
| chr7R | 2 | 3 |
| chr8L | 3 | 1 |
| chr8R | 3 | 7 |
| chr9L | 0 | 8 |
| **TOTAL** | **91** | **91** |