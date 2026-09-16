# Y' grouping schemes scored on 6991 day-0 reads

Day-0 populations should carry no recombination, so every mismatch below is the
grouping getting it wrong. A FALSE CALL is a read whose Y' array string does not
match its own anchor's reference array -- what the pipeline would report as
## POOLED
11583 non-recombinant reads across 2 day-0 populations

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 36 | 11583 | 15219 | 3624 | 23.81 | 1644 | 14.19 | 1232 | 10.64 |
| condensed | 18 | 11583 | 15219 | 918 | 6.03 | 918 | 7.93 | 915 | 7.90 |
| cut99 | 13 | 11583 | 15219 | 908 | 5.97 | 908 | 7.84 | 908 | 7.84 |
| cut97 | 10 | 11583 | 15219 | 18 | 0.12 | 18 | 0.16 | 18 | 0.16 |
| curated_variant | 12 | 11583 | 15219 | 267 | 1.75 | 267 | 2.31 | 267 | 2.31 |
| silhouette | 8 | 11583 | 15219 | 18 | 0.12 | 18 | 0.16 | 18 | 0.16 |
| curated_family | 7 | 11583 | 15219 | 26 | 0.17 | 26 | 0.22 | 26 | 0.22 |

### where element sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr7R | chr16L | 644 |
| chr12R | chr4R | 431 |
| chr14L | chr16L | 283 |
| chr4R | chr14L | 172 |
| chr12R | chr14L | 161 |
| chr14L | chr4R | 157 |
| chr16L | chr4R | 143 |
| chr4R | chr16L | 120 |

### where condensed sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr7R | chr16L | 644 |
| chr14L | chr16L | 237 |
| chr15R | chr16L | 9 |
| chr6L | chr14L | 5 |
| chr13L | chr12R|chr14L|chr4R | 3 |
| chr15R | chr12R|chr14L|chr16L|chr4R | 2 |
| chr16R | chr14L | 2 |
| chr14R | chr14L|chr7R | 2 |

### where cut99 sends the foreign-looking copies

| anchored end | apparent donor group lives at | copies |
|---|---|---|
| chr7R | chr16L | 644 |
| chr14L | chr16L | 237 |
| chr15R | chr16L | 9 |
| chr6L | chr13L|chr14L | 5 |
| chr16R | chr13L|chr14L | 2 |
| chr14R | chr12R|chr13L|chr14L|chr15R|chr4R|chr7R | 2 |
| chr5R | chr13L|chr14L | 1 |
| chr10L | chr14R | 1 |

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
| chr14L | chr16L|chr7R | 237 |
| chr15R | chr16L|chr7R | 9 |
| chr6L | chr14L | 5 |
| chr15R | chr12R|chr13L|chr14L|chr4R | 3 |
| chr16R | chr14L | 2 |
| chr14R | chr14L | 2 |
| chr13L | chr14L | 1 |
| chr10L | chr14R | 1 |

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
| chr15R | chr14L|chr16L|chr7R | 9 |
| chr6L | chr14L | 5 |
| chr16R | chr14L | 2 |
| chr14R | chr14L|chr16L|chr7R | 2 |
| chr5R | chr10L|chr14R|chr5L|chr9L | 2 |
| chr13L | chr14L | 1 |
| chr2L | chr12L|chr12R|chr16R|chr5R|chr8L|chr8R | 1 |
| chr6L | chr14L|chr16L|chr7R | 1 |

recombination on a read where nothing happened.

## 7172_day0_with_selection
4471 non-recombinant reads scored; 4846 set aside (end has no reference Y' array: 4766, no Y' copy detected: 22, copy count 2 vs reference 1: 9)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 35 | 4471 | 5630 | 1218 | 21.63 | 595 | 13.31 | 514 | 11.50 |
| condensed | 19 | 4471 | 5630 | 367 | 6.52 | 367 | 8.21 | 366 | 8.19 |
| cut99 | 13 | 4471 | 5630 | 362 | 6.43 | 362 | 8.10 | 362 | 8.10 |
| cut97 | 10 | 4471 | 5630 | 4 | 0.07 | 4 | 0.09 | 4 | 0.09 |
| curated_variant | 9 | 4471 | 5630 | 133 | 2.36 | 133 | 2.97 | 133 | 2.97 |
| silhouette | 8 | 4471 | 5630 | 4 | 0.07 | 4 | 0.09 | 4 | 0.09 |
| curated_family | 6 | 4471 | 5630 | 6 | 0.11 | 6 | 0.13 | 6 | 0.13 |

## 7302_day0_with_selection
7112 non-recombinant reads scored; 7090 set aside (end has no reference Y' array: 7004, no Y' copy detected: 30, copy count 2 vs reference 1: 8)

| scheme | groups | reads | copies | copy errors | copy err % | FALSE CALLS | false call % | of those, foreign-donor | foreign % |
|---|---|---|---|---|---|---|---|---|---|
| element | 36 | 7112 | 9589 | 2406 | 25.09 | 1049 | 14.75 | 718 | 10.10 |
| condensed | 18 | 7112 | 9589 | 551 | 5.75 | 551 | 7.75 | 549 | 7.72 |
| cut99 | 13 | 7112 | 9589 | 546 | 5.69 | 546 | 7.68 | 546 | 7.68 |
| cut97 | 10 | 7112 | 9589 | 14 | 0.15 | 14 | 0.20 | 14 | 0.20 |
| curated_variant | 12 | 7112 | 9589 | 134 | 1.40 | 134 | 1.88 | 134 | 1.88 |
| silhouette | 8 | 7112 | 9589 | 14 | 0.15 | 14 | 0.20 | 14 | 0.20 |
| curated_family | 7 | 7112 | 9589 | 20 | 0.21 | 20 | 0.28 | 20 | 0.28 |

