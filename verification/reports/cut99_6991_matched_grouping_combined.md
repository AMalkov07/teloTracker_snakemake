# 6991 combined: the five samples sharing an identical cut99 grouping

Restricted to samples whose cut99 partition (34 elements -> 12 groups) is EXACTLY
identical, verified by comparing group membership directly (not just group count).

## Which samples qualify

| sample | elements | groups | included |
|---|---|---|---|
| 6991_day0 | 34 | 12 | yes |
| 6991_day0_TeloTag_with_selection | 34 | 12 | yes |
| 6991_day0_reference_promethion | 34 | 12 | yes |
| 6991_day0_with_selection_repeat | 34 | 12 | yes |
| 6991_day0_with_selection_repeat2 | 34 | 12 | yes |
| 6991_day0_reference | 33 | 12 | **no** -- one fewer element (missing chr12R-7); same 12 groups otherwise but not an exact match |
| 6991_day0_TeloTag | 27 | 13 | **no** -- missing the whole chr4R array (structurally different) |
| 6991_day0_with_selection | 34 | 13 | **no** -- the chr14L-1 defect (5,720 bp) additionally splits off as its own group |

**5 of 8 samples qualify** (not 6 -- 6991_day0_reference is close but not exact).

## Group key (identical across all 5 samples)

| group | size | members |
|---|---|---|
| G8 | 18 | chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-1, chr14L-2, chr15R-1, chr16L-1, chr4R-1, chr4R-2, chr4R-3, chr4R-4, chr4R-5, chr4R-6, chr4R-7, chr7R-1 |
| G2 | 4 | chr13L-1, chr14L-3, chr14L-4, chr14L-5 |
| G7 | 2 | chr10L-1, chr9L-1 |
| G1 | 2 | chr2L-1, chr6L-1 |
| G4 | 1 | chr12L-1 |
| G11 | 1 | chr12R-1 |
| G10 | 1 | chr14R-1 |
| G6 | 1 | chr16R-1 |
| G12 | 1 | chr5L-1 |
| G9 | 1 | chr5R-1 |
| G5 | 1 | chr8L-1 |
| G3 | 1 | chr8R-1 |

## Pooled counts across the 5 samples (58 mismatched copies)

### As recipient

| end | count |
|---|---|
| chr13L | 18 |
| chr16R | 7 |
| chr10L | 6 |
| chr14R | 6 |
| chr2L | 5 |
| chr5L | 3 |
| chr8R | 2 |
| chr8L | 2 |
| chr12L | 2 |
| chr14L | 2 |
| chr6L | 1 |
| chr7R | 1 |
| chr15R | 1 |
| chr16L | 1 |
| chr5R | 1 |

### As donor, by group (unambiguous -- group labels are identical across all 5)

| group | times used as donor | group size |
|---|---|---|
| G1 | 18 | 2 |
| G2 | 13 | 4 |
| G8 | 10 | 18 |
| G7 | 5 | 2 |
| G11 | 4 | 1 |
| G3 | 4 | 1 |
| G10 | 2 | 1 |
| G6 | 1 | 1 |
| G5 | 1 | 1 |