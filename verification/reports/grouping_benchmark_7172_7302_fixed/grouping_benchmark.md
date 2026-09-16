# Y' grouping schemes scored on 6991 day-0 reads

Day-0 populations should carry no recombination, so every mismatch below is the
grouping getting it wrong. A FALSE CALL is a read whose Y' array string does not
match its own anchor's reference array -- what the pipeline would report as
## POOLED
11583 non-recombinant reads across 2 day-0 populations

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 36 | 11583 | 15219 | 2743 | 18.02 | 1000 | 8.63 | 546 | 4.71 |
| condensed | 17 | 11583 | 15219 | 29 | 0.19 | 29 | 0.25 | 26 | 0.22 |
| cut99 | 12 | 11583 | 15219 | 18 | 0.12 | 18 | 0.16 | 18 | 0.16 |
| cut97 | 10 | 11583 | 15219 | 18 | 0.12 | 18 | 0.16 | 18 | 0.16 |
| curated_variant | 12 | 11583 | 15219 | 23 | 0.15 | 23 | 0.20 | 23 | 0.20 |
| silhouette | 8 | 11583 | 15219 | 18 | 0.12 | 18 | 0.16 | 18 | 0.16 |
| curated_family | 7 | 11583 | 15219 | 18 | 0.12 | 18 | 0.16 | 18 | 0.16 |

### where element sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr12R | chr4R | 431 |
| chr4R | chr14L | 172 |
| chr12R | chr14L | 161 |
| chr14L | chr4R | 157 |
| chr16L | chr4R | 143 |
| chr4R | chr16L | 120 |
| chr12R | chr16L | 46 |
| chr14L | chr16L | 46 |

### where condensed sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr6L | chr14L | 5 |
| chr13L | chr12R|chr14L|chr4R | 3 |
| chr15R | chr12R|chr14L|chr16L|chr4R | 2 |
| chr16R | chr14L | 2 |
| chr14R | chr14L|chr16L|chr7R | 2 |
| chr13L | chr14L | 1 |
| chr16L | chr15R | 1 |
| chr5R | chr13L | 1 |

### where cut99 sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr6L | chr13L|chr14L | 5 |
| chr16R | chr13L|chr14L | 2 |
| chr14R | chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr7R | 2 |
| chr5R | chr13L|chr14L | 1 |
| chr10L | chr14R | 1 |
| chr2L | chr8R | 1 |
| chr5R | chr10L|chr9L | 1 |
| chr5R | chr14R | 1 |

### where cut97 sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr6L | chr13L|chr14L|chr8R | 6 |
| chr16R | chr13L|chr14L|chr8R | 2 |
| chr14R | chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 2 |
| chr5R | chr13L|chr14L|chr8R | 1 |
| chr10L | chr14R | 1 |
| chr2L | chr13L|chr14L|chr8R | 1 |
| chr5R | chr10L|chr9L | 1 |
| chr5R | chr14R | 1 |

### where curated_variant sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr6L | chr14L | 5 |
| chr15R | chr12R|chr13L|chr14L|chr4R | 4 |
| chr16R | chr14L | 2 |
| chr14R | chr14L | 2 |
| chr13L | chr14L | 1 |
| chr10L | chr14R | 1 |
| chr15R | chr16L|chr7R | 1 |
| chr2L | chr12L|chr12R|chr16R|chr5R|chr8L|chr8R | 1 |

### where silhouette sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr6L | chr12L|chr13L|chr14L|chr8L|chr8R | 6 |
| chr16R | chr12L|chr13L|chr14L|chr8L|chr8R | 2 |
| chr14R | chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 2 |
| chr5R | chr12L|chr13L|chr14L|chr8L|chr8R | 1 |
| chr10L | chr14R | 1 |
| chr2L | chr12L|chr13L|chr14L|chr8L|chr8R | 1 |
| chr5R | chr10L|chr9L | 1 |
| chr5R | chr14R | 1 |

### where curated_family sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr6L | chr14L | 5 |
| chr16R | chr14L | 2 |
| chr14R | chr14L|chr16L|chr7R | 2 |
| chr5R | chr10L|chr14R|chr5L|chr9L | 2 |
| chr13L | chr14L | 1 |
| chr15R | chr14L|chr16L|chr7R | 1 |
| chr2L | chr12L|chr12R|chr16R|chr5R|chr8L|chr8R | 1 |
| chr6L | chr14L|chr16L|chr7R | 1 |

recombination on a read where nothing happened.

## 7172_day0_with_selection
4471 non-recombinant reads scored; 4846 set aside (end has no reference Y' array: 4766, no Y' copy detected: 22, copy count 2 vs reference 1: 9)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 35 | 4471 | 5630 | 861 | 15.29 | 365 | 8.16 | 284 | 6.35 |
| condensed | 18 | 4471 | 5630 | 9 | 0.16 | 9 | 0.20 | 8 | 0.18 |
| cut99 | 12 | 4471 | 5630 | 4 | 0.07 | 4 | 0.09 | 4 | 0.09 |
| cut97 | 10 | 4471 | 5630 | 4 | 0.07 | 4 | 0.09 | 4 | 0.09 |
| curated_variant | 9 | 4471 | 5630 | 4 | 0.07 | 4 | 0.09 | 4 | 0.09 |
| silhouette | 8 | 4471 | 5630 | 4 | 0.07 | 4 | 0.09 | 4 | 0.09 |
| curated_family | 6 | 4471 | 5630 | 4 | 0.07 | 4 | 0.09 | 4 | 0.09 |

## 7302_day0_with_selection
7112 non-recombinant reads scored; 7090 set aside (end has no reference Y' array: 7004, no Y' copy detected: 30, copy count 2 vs reference 1: 8)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 36 | 7112 | 9589 | 1882 | 19.63 | 635 | 8.93 | 262 | 3.68 |
| condensed | 17 | 7112 | 9589 | 20 | 0.21 | 20 | 0.28 | 18 | 0.25 |
| cut99 | 12 | 7112 | 9589 | 14 | 0.15 | 14 | 0.20 | 14 | 0.20 |
| cut97 | 10 | 7112 | 9589 | 14 | 0.15 | 14 | 0.20 | 14 | 0.20 |
| curated_variant | 12 | 7112 | 9589 | 19 | 0.20 | 19 | 0.27 | 19 | 0.27 |
| silhouette | 8 | 7112 | 9589 | 14 | 0.15 | 14 | 0.20 | 14 | 0.20 |
| curated_family | 7 | 7112 | 9589 | 14 | 0.15 | 14 | 0.20 | 14 | 0.20 |


---

# Before / after the chr16L boundary fix — 7172 and 7302 day-0

Same locus, same 77 bp anchor-proximal over-extension found in all eight 6991 day-0
references, confirmed independently on both strains' day-0 populations
(`7172_day0_with_selection`, `7302_day0_with_selection` — the only day-0 population
available on Argon for each). 11,583 non-recombinant reads pooled across the two:

| scheme | groups before -> after | false call % before | false call % after |
|---|---|---|---|
| element | 36 -> 36 | 14.19 | 8.63 |
| condensed | 18 -> 17 | 7.93 | 0.25 |
| cut99 | 13 -> 12 | 7.84 | 0.16 |
| cut97 | 10 -> 10 | 0.16 | 0.16 |
| curated_variant | 12 -> 12 | 2.31 | **0.20** |
| silhouette | 8 -> 8 | 0.16 | 0.16 |
| curated_family | 7 -> 7 | 0.22 | 0.16 |

Same pattern as 6991: `condensed` and `cut99` drop by a factor of ~30-50 and land in the
same band as every other scheme; `cut97`/`silhouette` are unchanged since they already
grouped chr7R/chr14L/chr16L together.

**`curated_variant` fully recovers here, unlike 6991.** Both strains' curated libraries carry
the same `chr7R1`/`chr16L1 = ID5_Blue-Dark`, `chr14L1 = ID5_Blue-Light` split as 6991's, so the
label disagreement the boundary fix cannot touch is present in exactly the same form. But in
6991 that split alone left 442 residual chr14L -> chr16L/chr7R miscalls (the dominant term in
`curated_variant`'s 0.81 %); here the equivalent confusion barely registers (0-1 copies) --
7172/7302's chr14L arrays evidently don't put many anchored reads' first copy through this
particular Blue-Light/Blue-Dark ambiguity, so the fix takes `curated_variant` down to the
same 0.16-0.20 % floor as everything else on these two strains.
