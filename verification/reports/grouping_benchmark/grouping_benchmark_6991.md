# Y' grouping schemes scored on 6991 day-0 reads

Day-0 populations should carry no recombination, so every mismatch below is the
grouping getting it wrong. A FALSE CALL is a read whose Y' array string does not
match its own anchor's reference array -- what the pipeline would report as
## POOLED
52182 non-recombinant reads across 7 day-0 populations

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 52182 | 59132 | 9540 | 16.13 | 4510 | 8.64 | 4269 | 8.18 |
| condensed | 18 | 52182 | 59132 | 3659 | 6.19 | 3658 | 7.01 | 3639 | 6.97 |
| cut99 | 13 | 52182 | 59132 | 3611 | 6.11 | 3611 | 6.92 | 3610 | 6.92 |
| cut97 | 10 | 52182 | 59132 | 87 | 0.15 | 87 | 0.17 | 86 | 0.16 |
| curated_variant | 10 | 52182 | 59132 | 669 | 1.13 | 669 | 1.28 | 668 | 1.28 |
| silhouette | 8 | 52182 | 59132 | 82 | 0.14 | 82 | 0.16 | 81 | 0.16 |
| curated_family | 6 | 52182 | 59132 | 140 | 0.24 | 140 | 0.27 | 139 | 0.27 |

### where element sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr7R | chr16L | 2968 |
| chr12R | chr4R | 1002 |
| chr4R | chr14L | 969 |
| chr14L | chr16L | 493 |
| chr12R | chr14L | 485 |
| chr14L | chr4R | 376 |
| chr15R | chr16L | 59 |
| chr14L | chr12R | 38 |

### where condensed sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr7R | chr16L | 2968 |
| chr14L | chr16L | 493 |
| chr15R | chr16L | 59 |
| chr13L | chr6L | 10 |
| chr13L | chr14L | 9 |
| chr13L | chr2L | 7 |
| chr15R | chr12R|chr14L|chr4R | 7 |
| chr2L | chr8R | 6 |

### where cut99 sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr7R | chr16L | 2968 |
| chr14L | chr16L | 493 |
| chr15R | chr16L | 59 |
| chr13L | chr2L|chr6L | 17 |
| chr16R | chr13L|chr14L | 9 |
| chr2L | chr8R | 6 |
| chr6L | chr13L|chr14L | 5 |
| chr14R | chr10L|chr9L | 5 |

### where cut97 sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr13L | chr2L|chr6L | 17 |
| chr16R | chr13L|chr14L|chr8R | 9 |
| chr2L | chr13L|chr14L|chr8R | 7 |
| chr10L | chr12R|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 6 |
| chr6L | chr13L|chr14L|chr8R | 5 |
| chr14R | chr10L|chr9L | 5 |
| chr5L | chr12R | 3 |
| chr12L | chr13L|chr14L|chr8R | 3 |

### where curated_variant sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr14L | chr16L|chr7R | 507 |
| chr15R | chr16L|chr7R | 59 |
| chr13L | chr2L|chr6L | 17 |
| chr13L | chr14L | 9 |
| chr15R | chr12R|chr14L|chr4R | 7 |
| chr2L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | 6 |
| chr10L | chr16L|chr7R | 5 |
| chr16R | chr14L | 5 |

### where silhouette sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr13L | chr2L|chr6L | 17 |
| chr16R | chr12L|chr13L|chr14L|chr8L|chr8R | 9 |
| chr2L | chr12L|chr13L|chr14L|chr8L|chr8R | 7 |
| chr10L | chr12R|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 6 |
| chr6L | chr12L|chr13L|chr14L|chr8L|chr8R | 5 |
| chr14R | chr10L|chr9L | 5 |
| chr5L | chr12R | 3 |
| chr14R | chr12R|chr14L|chr15R|chr16L|chr5R|chr7R | 2 |

### where curated_family sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr15R | chr14L|chr16L|chr7R | 59 |
| chr13L | chr2L|chr6L | 17 |
| chr13L | chr14L | 9 |
| chr10L | chr14L|chr16L|chr7R | 6 |
| chr2L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | 6 |
| chr16R | chr14L | 5 |
| chr5L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | 4 |
| chr14R | chr14L|chr16L|chr7R | 3 |

recombination on a read where nothing happened.

## 6991_day0
3489 non-recombinant reads scored; 3756 set aside (end has no reference Y' array: 3692, no Y' copy detected: 23, copy count 2 vs reference 1: 6)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 3489 | 4361 | 1048 | 24.03 | 384 | 11.01 | 366 | 10.49 |
| condensed | 18 | 3489 | 4361 | 264 | 6.05 | 264 | 7.57 | 263 | 7.54 |
| cut99 | 13 | 3489 | 4361 | 261 | 5.98 | 261 | 7.48 | 260 | 7.45 |
| cut97 | 10 | 3489 | 4361 | 4 | 0.09 | 4 | 0.11 | 3 | 0.09 |
| curated_variant | 10 | 3489 | 4361 | 61 | 1.40 | 61 | 1.75 | 60 | 1.72 |
| silhouette | 8 | 3489 | 4361 | 4 | 0.09 | 4 | 0.11 | 3 | 0.09 |
| curated_family | 6 | 3489 | 4361 | 11 | 0.25 | 11 | 0.32 | 10 | 0.29 |

## 6991_day0_TeloTag
7428 non-recombinant reads scored; 8999 set aside (end has no reference Y' array: 8918, no Y' copy detected: 40, copy count 2 vs reference 1: 13)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 27 | 7428 | 7790 | 722 | 9.27 | 480 | 6.46 | 476 | 6.41 |
| condensed | 19 | 7428 | 7790 | 460 | 5.91 | 459 | 6.18 | 441 | 5.94 |
| cut99 | 14 | 7428 | 7790 | 434 | 5.57 | 434 | 5.84 | 434 | 5.84 |
| cut97 | 10 | 7428 | 7790 | 4 | 0.05 | 4 | 0.05 | 4 | 0.05 |
| curated_variant | 10 | 7428 | 7790 | 72 | 0.92 | 72 | 0.97 | 72 | 0.97 |
| silhouette | 8 | 7428 | 7790 | 4 | 0.05 | 4 | 0.05 | 4 | 0.05 |
| curated_family | 6 | 7428 | 7790 | 23 | 0.30 | 23 | 0.31 | 23 | 0.31 |

## 6991_day0_TeloTag_with_selection
8765 non-recombinant reads scored; 9316 set aside (end has no reference Y' array: 9194, no Y' copy detected: 71, copy count 2 vs reference 1: 12)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 8765 | 9355 | 1032 | 11.03 | 624 | 7.12 | 600 | 6.85 |
| condensed | 18 | 8765 | 9355 | 550 | 5.88 | 550 | 6.27 | 550 | 6.27 |
| cut99 | 13 | 8765 | 9355 | 547 | 5.85 | 547 | 6.24 | 547 | 6.24 |
| cut97 | 10 | 8765 | 9355 | 10 | 0.11 | 10 | 0.11 | 10 | 0.11 |
| curated_variant | 10 | 8765 | 9355 | 63 | 0.67 | 63 | 0.72 | 63 | 0.72 |
| silhouette | 8 | 8765 | 9355 | 10 | 0.11 | 10 | 0.11 | 10 | 0.11 |
| curated_family | 6 | 8765 | 9355 | 22 | 0.24 | 22 | 0.25 | 22 | 0.25 |

## 6991_day0_reference
7511 non-recombinant reads scored; 7881 set aside (end has no reference Y' array: 7684, no Y' copy detected: 68, copy count 2 vs reference 1: 40)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 33 | 7511 | 8573 | 1434 | 16.73 | 673 | 8.96 | 635 | 8.45 |
| condensed | 18 | 7511 | 8573 | 532 | 6.21 | 532 | 7.08 | 532 | 7.08 |
| cut99 | 13 | 7511 | 8573 | 523 | 6.10 | 523 | 6.96 | 523 | 6.96 |
| cut97 | 10 | 7511 | 8573 | 27 | 0.31 | 27 | 0.36 | 27 | 0.36 |
| curated_variant | 10 | 7511 | 8573 | 124 | 1.45 | 124 | 1.65 | 124 | 1.65 |
| silhouette | 8 | 7511 | 8573 | 26 | 0.30 | 26 | 0.35 | 26 | 0.35 |
| curated_family | 6 | 7511 | 8573 | 38 | 0.44 | 38 | 0.51 | 38 | 0.51 |

## 6991_day0_reference_promethion
7181 non-recombinant reads scored; 8540 set aside (end has no reference Y' array: 8393, no Y' copy detected: 73, copy count 2 vs reference 1: 20)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 7181 | 8509 | 1634 | 19.20 | 646 | 9.00 | 609 | 8.48 |
| condensed | 18 | 7181 | 8509 | 462 | 5.43 | 462 | 6.43 | 462 | 6.43 |
| cut99 | 13 | 7181 | 8509 | 456 | 5.36 | 456 | 6.35 | 456 | 6.35 |
| cut97 | 10 | 7181 | 8509 | 11 | 0.13 | 11 | 0.15 | 11 | 0.15 |
| curated_variant | 10 | 7181 | 8509 | 89 | 1.05 | 89 | 1.24 | 89 | 1.24 |
| silhouette | 8 | 7181 | 8509 | 10 | 0.12 | 10 | 0.14 | 10 | 0.14 |
| curated_family | 6 | 7181 | 8509 | 13 | 0.15 | 13 | 0.18 | 13 | 0.18 |

## 6991_day0_with_selection_repeat
12417 non-recombinant reads scored; 11932 set aside (end has no reference Y' array: 11731, no Y' copy detected: 68, copy count 2 vs reference 1: 39)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 12417 | 14461 | 2688 | 18.59 | 1223 | 9.85 | 1132 | 9.12 |
| condensed | 18 | 12417 | 14461 | 991 | 6.85 | 991 | 7.98 | 991 | 7.98 |
| cut99 | 13 | 12417 | 14461 | 990 | 6.85 | 990 | 7.97 | 990 | 7.97 |
| cut97 | 10 | 12417 | 14461 | 18 | 0.12 | 18 | 0.14 | 18 | 0.14 |
| curated_variant | 10 | 12417 | 14461 | 190 | 1.31 | 190 | 1.53 | 190 | 1.53 |
| silhouette | 8 | 12417 | 14461 | 17 | 0.12 | 17 | 0.14 | 17 | 0.14 |
| curated_family | 6 | 12417 | 14461 | 18 | 0.12 | 18 | 0.14 | 18 | 0.14 |

## 6991_day0_with_selection_repeat2
5391 non-recombinant reads scored; 5473 set aside (end has no reference Y' array: 5402, no Y' copy detected: 31, copy count 2 vs reference 1: 16)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 5391 | 6083 | 982 | 16.14 | 480 | 8.90 | 451 | 8.37 |
| condensed | 18 | 5391 | 6083 | 400 | 6.58 | 400 | 7.42 | 400 | 7.42 |
| cut99 | 13 | 5391 | 6083 | 400 | 6.58 | 400 | 7.42 | 400 | 7.42 |
| cut97 | 10 | 5391 | 6083 | 13 | 0.21 | 13 | 0.24 | 13 | 0.24 |
| curated_variant | 10 | 5391 | 6083 | 70 | 1.15 | 70 | 1.30 | 70 | 1.30 |
| silhouette | 8 | 5391 | 6083 | 11 | 0.18 | 11 | 0.20 | 11 | 0.20 |
| curated_family | 6 | 5391 | 6083 | 15 | 0.25 | 15 | 0.28 | 15 | 0.28 |


---

# The dominant error: chr7R / chr14L / chr16L

Almost every false call in the finer schemes comes from three long Y' elements that are the
same sequence to within 1-2 bp:

| pair | identity over the shared 6654 bp | note |
|---|---|---|
| chr7R_1 vs chr14L_1 | 99.970 % (~2 bp) | aligned end to end |
| chr7R_1 vs chr16L_1 | 99.985 % (~1 bp) | chr16L_1 carries **77 extra bp at its 5' end** |
| chr14L_1 vs chr16L_1 | 99.955 % (~3 bp) | same 77 bp offset |

The 77 bp is an annotation boundary difference, not biology: the BED calls the chr16L element
77 bp further toward the anchor than it calls chr7R and chr14L. Two things follow.

**1. Reads go to the longest entry.** Those 77 bp are on the reads as well, so the longer
library entry wins the RepeatMasker SW score: 806/806 chr7R reads and 192/204 chr14L
anchor-proximal copies match `E-chr16L-1`, never their own element. This is systematic, not a
coin flip between near-identical entries.

**2. The similarity formula turns the offset into a split.** Similarity is coverage-penalised,
so 77 unaligned bp drop chr16L_1 to 98.8 % against the other two -- below a 99 % cut. `cut99`
and `condensed` therefore put chr16L in its own group, and that phantom group then captures
essentially every chr7R and chr14L read.

How each scheme carves the three, and what it costs:

| scheme | chr7R_1 | chr16L_1 | chr14L_1 | consequence |
|---|---|---|---|---|
| condensed / cut99 | A | **B** | A | chr7R 100 % wrong, chr14L copy 1 94 % wrong |
| cut97 / silhouette | A | A | A | no error at this locus |
| curated_variant | ID5_Blue-Dark | ID5_Blue-Dark | **ID5_Blue-Light** | chr7R fine, chr14L wrong (507 copies) |

The curated row is the chr14L `Blue-Dark` / `Blue-Light` artefact already seen in the
Supplementary Data 5 comparison, where it produced most of the spurious wild-type switch
calls. Its mechanism is now known: the curated scheme asserts a Blue-Light/Blue-Dark
distinction that rests on 1-2 bp across 6.6 kb -- below ONT read accuracy, so it is not
reliably readable from a read no matter how the grouping is done.

## False-call rate with that locus set aside

| scheme | groups | all reads | excluding chr7R/chr14L/chr16L |
|---|---|---|---|
| element | 34 | 8.64 % | 1.97 % |
| condensed | 18 | 7.01 % | 0.41 % |
| cut99 | 13 | 6.92 % | **0.31 %** |
| cut97 | 10 | 0.17 % | 0.18 % |
| curated_variant | 10 | 1.28 % | 0.33 % |
| silhouette | 8 | 0.16 % | 0.17 % |
| curated_family | 6 | 0.27 % | 0.29 % |

Away from that one locus every scheme sits between 0.17 % and 0.41 %, and `cut99` -- 13 groups
against silhouette's 8, and 12 of 18 chromosome ends uniquely identified against silhouette's
6 -- is the best of the finer ones. The extra resolution is close to free; it is the boundary
inconsistency, not the finer grouping, that makes it expensive.
