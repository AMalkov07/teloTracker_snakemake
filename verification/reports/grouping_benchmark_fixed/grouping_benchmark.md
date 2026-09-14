# Y' grouping schemes scored on 6991 day-0 reads

Day-0 populations should carry no recombination, so every mismatch below is the
grouping getting it wrong. A FALSE CALL is a read whose Y' array string does not
match its own anchor's reference array -- what the pipeline would report as
## POOLED
79593 non-recombinant reads across 8 day-0 populations

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 79593 | 91359 | 11383 | 12.46 | 3187 | 4.00 | 2612 | 3.28 |
| condensed | 17 | 79593 | 91359 | 602 | 0.66 | 599 | 0.75 | 580 | 0.73 |
| cut99 | 12 | 79593 | 91359 | 530 | 0.58 | 530 | 0.67 | 141 | 0.18 |
| cut97 | 10 | 79593 | 91359 | 527 | 0.58 | 527 | 0.66 | 138 | 0.17 |
| curated_variant | 10 | 79593 | 91359 | 650 | 0.71 | 648 | 0.81 | 647 | 0.81 |
| silhouette | 8 | 79593 | 91359 | 132 | 0.14 | 132 | 0.17 | 131 | 0.16 |
| curated_family | 6 | 79593 | 91359 | 134 | 0.15 | 134 | 0.17 | 133 | 0.17 |

### where element sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr12R | chr4R | 2081 |
| chr4R | chr14L | 1324 |
| chr12R | chr14L | 760 |
| chr14L | chr4R | 671 |
| chr14L | chr7R | 439 |
| chr7R | chr16L | 354 |
| chr16L | chr7R | 317 |
| chr14L | chr12R | 38 |

### where condensed sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr14L | chr16L|chr7R | 388 |
| chr13L | chr14L | 16 |
| chr13L | chr2L | 15 |
| chr13L | chr6L | 14 |
| chr4R | chr15R | 10 |
| chr14R | chr9L | 9 |
| chr15R | chr12R|chr14L|chr4R | 9 |
| chr2L | chr8R | 9 |

### where cut99 sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr13L | chr2L|chr6L | 29 |
| chr16R | chr13L|chr14L | 11 |
| chr10L | chr12R|chr14L|chr15R|chr16L|chr4R|chr7R | 9 |
| chr14R | chr10L|chr9L | 9 |
| chr2L | chr8R | 9 |
| chr6L | chr13L|chr14L | 8 |
| chr14R | chr12R|chr14L|chr15R|chr16L|chr4R|chr7R | 6 |
| chr13L | chr12R|chr14L|chr15R|chr16L|chr4R|chr7R | 5 |

### where cut97 sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr13L | chr2L|chr6L | 29 |
| chr16R | chr13L|chr14L|chr8R | 11 |
| chr10L | chr12R|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 10 |
| chr2L | chr13L|chr14L|chr8R | 10 |
| chr14R | chr10L|chr9L | 9 |
| chr6L | chr13L|chr14L|chr8R | 8 |
| chr13L | chr12R|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 6 |
| chr14R | chr12R|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 6 |

### where curated_variant sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr14L | chr16L|chr7R | 442 |
| chr13L | chr2L|chr6L | 29 |
| chr7R | chr14L | 21 |
| chr13L | chr14L | 16 |
| chr16L | chr14L | 15 |
| chr4R | chr15R | 10 |
| chr14R | chr10L|chr9L | 9 |
| chr15R | chr12R|chr14L|chr4R | 9 |

### where silhouette sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr13L | chr2L|chr6L | 29 |
| chr16R | chr12L|chr13L|chr14L|chr8L|chr8R | 11 |
| chr10L | chr12R|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 10 |
| chr2L | chr12L|chr13L|chr14L|chr8L|chr8R | 10 |
| chr14R | chr10L|chr9L | 9 |
| chr6L | chr12L|chr13L|chr14L|chr8L|chr8R | 8 |
| chr13L | chr12R|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 6 |
| chr14R | chr12R|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | 6 |

### where curated_family sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr13L | chr2L|chr6L | 29 |
| chr13L | chr14L | 16 |
| chr10L | chr14L|chr16L|chr7R | 9 |
| chr2L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | 9 |
| chr14R | chr14L|chr16L|chr7R | 7 |
| chr16R | chr14L | 6 |
| chr6L | chr14L | 6 |
| chr5L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | 4 |

recombination on a read where nothing happened.

## 6991_day0
3489 non-recombinant reads scored; 3756 set aside (end has no reference Y' array: 3692, no Y' copy detected: 23, copy count 2 vs reference 1: 6)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 3489 | 4361 | 833 | 19.10 | 186 | 5.33 | 164 | 4.70 |
| condensed | 17 | 3489 | 4361 | 7 | 0.16 | 7 | 0.20 | 6 | 0.17 |
| cut99 | 12 | 3489 | 4361 | 4 | 0.09 | 4 | 0.11 | 3 | 0.09 |
| cut97 | 10 | 3489 | 4361 | 4 | 0.09 | 4 | 0.11 | 3 | 0.09 |
| curated_variant | 10 | 3489 | 4361 | 42 | 0.96 | 42 | 1.20 | 41 | 1.18 |
| silhouette | 8 | 3489 | 4361 | 4 | 0.09 | 4 | 0.11 | 3 | 0.09 |
| curated_family | 6 | 3489 | 4361 | 5 | 0.11 | 5 | 0.14 | 4 | 0.11 |

## 6991_day0_TeloTag
7428 non-recombinant reads scored; 8999 set aside (end has no reference Y' array: 8918, no Y' copy detected: 40, copy count 2 vs reference 1: 13)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 27 | 7428 | 7790 | 595 | 7.64 | 382 | 5.14 | 365 | 4.91 |
| condensed | 18 | 7428 | 7790 | 32 | 0.41 | 31 | 0.42 | 13 | 0.18 |
| cut99 | 13 | 7428 | 7790 | 6 | 0.08 | 6 | 0.08 | 6 | 0.08 |
| cut97 | 11 | 7428 | 7790 | 5 | 0.06 | 5 | 0.07 | 5 | 0.07 |
| curated_variant | 10 | 7428 | 7790 | 40 | 0.51 | 40 | 0.54 | 40 | 0.54 |
| silhouette | 8 | 7428 | 7790 | 4 | 0.05 | 4 | 0.05 | 4 | 0.05 |
| curated_family | 6 | 7428 | 7790 | 7 | 0.09 | 7 | 0.09 | 7 | 0.09 |

## 6991_day0_TeloTag_with_selection
8765 non-recombinant reads scored; 9316 set aside (end has no reference Y' array: 9194, no Y' copy detected: 71, copy count 2 vs reference 1: 12)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 8765 | 9355 | 519 | 5.55 | 150 | 1.71 | 108 | 1.23 |
| condensed | 17 | 8765 | 9355 | 15 | 0.16 | 15 | 0.17 | 15 | 0.17 |
| cut99 | 12 | 8765 | 9355 | 11 | 0.12 | 11 | 0.13 | 11 | 0.13 |
| cut97 | 10 | 8765 | 9355 | 11 | 0.12 | 11 | 0.13 | 11 | 0.13 |
| curated_variant | 10 | 8765 | 9355 | 17 | 0.18 | 17 | 0.19 | 17 | 0.19 |
| silhouette | 8 | 8765 | 9355 | 11 | 0.12 | 11 | 0.13 | 11 | 0.13 |
| curated_family | 6 | 8765 | 9355 | 13 | 0.14 | 13 | 0.15 | 13 | 0.15 |

## 6991_day0_reference
7511 non-recombinant reads scored; 7881 set aside (end has no reference Y' array: 7684, no Y' copy detected: 68, copy count 2 vs reference 1: 40)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 33 | 7511 | 8573 | 1217 | 14.20 | 530 | 7.06 | 472 | 6.28 |
| condensed | 17 | 7511 | 8573 | 37 | 0.43 | 37 | 0.49 | 37 | 0.49 |
| cut99 | 12 | 7511 | 8573 | 27 | 0.31 | 27 | 0.36 | 27 | 0.36 |
| cut97 | 10 | 7511 | 8573 | 27 | 0.31 | 27 | 0.36 | 27 | 0.36 |
| curated_variant | 10 | 7511 | 8573 | 40 | 0.47 | 40 | 0.53 | 40 | 0.53 |
| silhouette | 8 | 7511 | 8573 | 26 | 0.30 | 26 | 0.35 | 26 | 0.35 |
| curated_family | 6 | 7511 | 8573 | 24 | 0.28 | 24 | 0.32 | 24 | 0.32 |

## 6991_day0_reference_promethion
7180 non-recombinant reads scored; 8541 set aside (end has no reference Y' array: 8393, no Y' copy detected: 73, copy count 2 vs reference 1: 21)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 7180 | 8508 | 1206 | 14.17 | 288 | 4.01 | 231 | 3.22 |
| condensed | 17 | 7180 | 8508 | 18 | 0.21 | 18 | 0.25 | 18 | 0.25 |
| cut99 | 12 | 7180 | 8508 | 12 | 0.14 | 12 | 0.17 | 12 | 0.17 |
| cut97 | 10 | 7180 | 8508 | 11 | 0.13 | 11 | 0.15 | 11 | 0.15 |
| curated_variant | 10 | 7180 | 8508 | 18 | 0.21 | 18 | 0.25 | 18 | 0.25 |
| silhouette | 8 | 7180 | 8508 | 10 | 0.12 | 10 | 0.14 | 10 | 0.14 |
| curated_family | 6 | 7180 | 8508 | 10 | 0.12 | 10 | 0.14 | 10 | 0.14 |

## 6991_day0_with_selection
27409 non-recombinant reads scored; 26064 set aside (end has no reference Y' array: 25609, no Y' copy detected: 152, copy count 2 vs reference 1: 77)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 27409 | 32225 | 4680 | 14.52 | 1061 | 3.87 | 887 | 3.24 |
| condensed | 18 | 27409 | 32225 | 461 | 1.43 | 459 | 1.67 | 459 | 1.67 |
| cut99 | 13 | 27409 | 32225 | 439 | 1.36 | 439 | 1.60 | 51 | 0.19 |
| cut97 | 11 | 27409 | 32225 | 438 | 1.36 | 438 | 1.60 | 50 | 0.18 |
| curated_variant | 10 | 27409 | 32225 | 454 | 1.41 | 452 | 1.65 | 452 | 1.65 |
| silhouette | 8 | 27409 | 32225 | 49 | 0.15 | 49 | 0.18 | 49 | 0.18 |
| curated_family | 6 | 27409 | 32225 | 50 | 0.16 | 50 | 0.18 | 50 | 0.18 |

## 6991_day0_with_selection_repeat
12420 non-recombinant reads scored; 11929 set aside (end has no reference Y' array: 11731, no Y' copy detected: 68, copy count 2 vs reference 1: 36)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 12420 | 14464 | 1731 | 11.97 | 435 | 3.50 | 277 | 2.23 |
| condensed | 17 | 12420 | 14464 | 19 | 0.13 | 19 | 0.15 | 19 | 0.15 |
| cut99 | 12 | 12420 | 14464 | 18 | 0.12 | 18 | 0.14 | 18 | 0.14 |
| cut97 | 10 | 12420 | 14464 | 18 | 0.12 | 18 | 0.14 | 18 | 0.14 |
| curated_variant | 10 | 12420 | 14464 | 26 | 0.18 | 26 | 0.21 | 26 | 0.21 |
| silhouette | 8 | 12420 | 14464 | 17 | 0.12 | 17 | 0.14 | 17 | 0.14 |
| curated_family | 6 | 12420 | 14464 | 15 | 0.10 | 15 | 0.12 | 15 | 0.12 |

## 6991_day0_with_selection_repeat2
5391 non-recombinant reads scored; 5473 set aside (end has no reference Y' array: 5402, no Y' copy detected: 31, copy count 2 vs reference 1: 16)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 34 | 5391 | 6083 | 602 | 9.90 | 155 | 2.88 | 108 | 2.00 |
| condensed | 17 | 5391 | 6083 | 13 | 0.21 | 13 | 0.24 | 13 | 0.24 |
| cut99 | 12 | 5391 | 6083 | 13 | 0.21 | 13 | 0.24 | 13 | 0.24 |
| cut97 | 10 | 5391 | 6083 | 13 | 0.21 | 13 | 0.24 | 13 | 0.24 |
| curated_variant | 10 | 5391 | 6083 | 13 | 0.21 | 13 | 0.24 | 13 | 0.24 |
| silhouette | 8 | 5391 | 6083 | 11 | 0.18 | 11 | 0.20 | 11 | 0.20 |
| curated_family | 6 | 5391 | 6083 | 10 | 0.16 | 10 | 0.19 | 10 | 0.19 |


---

# Before / after the chr16L boundary fix (same 8 samples, 79,590-93 reads)

| scheme | groups before -> after | false call % before | false call % after |
|---|---|---|---|
| element | 34 -> 34 | 8.95 | 4.00 |
| condensed | 18 -> 17 | 7.21 | 0.75 |
| cut99 | 13 -> 12 | 7.12 | 0.67 |
| cut97 | 10 -> 10 | 0.66 | 0.66 |
| curated_variant | 10 -> 10 | 1.44 | 0.81 |
| silhouette | 8 -> 8 | 0.16 | 0.17 |
| curated_family | 6 -> 6 | 0.27 | 0.17 |

`condensed` and `cut99` -- the two schemes whose false calls were almost entirely the
chr7R/chr14L -> chr16L confusion -- drop by a factor of ~9-10, landing at 0.75 % / 0.67 %, in
the same band as every other scheme. `cut97` and `silhouette` are unchanged, as expected: they
already grouped the three together before the fix. `element` (unlabelled per-element matching,
no grouping at all) drops from 8.95 % to 4.00 % -- the fix corrects the matching itself, not
just how matches are grouped afterwards, but element level is inherently the floor (see the
POOLED confusion tables above for what still confuses at that resolution).

`curated_variant` improves (1.44 % -> 0.81 %) but does not reach the other schemes' floor,
because its remaining error is a different, unrelated locus: it labels chr7R1/chr16L1
`ID5_Blue-Dark` and chr14L1 `ID5_Blue-Light`, a distinction rests on 1-2 bp across 6.6 kb (see
`verification/curated_refs/README.md`) that the boundary fix does not touch and that no
grouping scheme can read reliably from an ONT read.
