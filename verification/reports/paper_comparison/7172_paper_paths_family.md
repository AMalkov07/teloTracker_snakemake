# Summary

| PD | reads | paper_Y | caught | missed | extra | agreement_pct |
|---|---|---|---|---|---|---|
| 28.0 | 60.0 | 9.0 | 6.0 | 3.0 | 6.0 | 85.0 |
| 33.0 | 125.0 | 24.0 | 23.0 | 1.0 | 14.0 | 88.0 |


# Our path parser applied to the paper's own Y' calls — strain 7172

reference: `7172_day0_with_selection`, curated 7172 library; orientation inferred from the data.

## PD 28 — 60 reads
* paper flags 9 switching reads; we catch 6, miss 3, and flag 6 it calls N (85.0 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr6L | ID3,ID3,ID5 | chr14L | chr14L[3]:ID3 > chr14L[3]:ID3 > chr14L|chr16L|chr7R:ID5 |
| chr1R | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L[3]:ID3 |
| chr15L | ID5,ID2,ID3,ID4 | chr14L > chr2L|chr6L | chr14L|chr16L|chr7R:ID5 > chr14L[2-3]:ID2,ID3 > chr2L|chr6L:ID4 |
| chr11R | ID4,ID3,ID3,ID3,ID3,ID3 | chr2L|chr6L > chr14L | chr2L|chr6L:ID4 > chr14L[3]:ID3,ID3,ID3,ID3,ID3(circ x5.0 strong) |
| chr11L | ID6,ID3,ID3,ID3,ID3 | chr10L|chr14R|chr5L|chr9L > chr14L | chr10L|chr14R|chr5L|chr9L:ID6 > chr14L[3]:ID3 > chr14L[3]:ID3,ID3,ID3(circ x3.0 moderate) |
| chr16R | ID1,ID2,ID2,ID6 | chr12R|chr4R > chr10L|chr14R|chr5L|chr9L | chr12R|chr4R:ID2,ID2 > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr9L | ID6,ID5,ID3,ID3,ID5 | chr14L | chr14L|chr16L:ID5 > chr14L[3]:ID3 > chr14L[3]:ID3 > chr14L|chr16L|chr7R:ID5 |
| chr1L | ID5,ID5,ID5 | chr14L|chr16L | chr14L|chr16L:ID5,ID5,ID5(circ x3.0) |
| chr13R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID6 | chr12R > chr10L|chr14R|chr5L|chr9L | chr12R[1-7]:ID1,ID2,ID2,ID2,ID2,ID2,ID2 > chr12R|chr4R:ID2,ID2 > chr10L|chr14R|chr5L|chr9L:ID6 |

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr16L | ID5,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr16L > chr12R > chr16L |
| chr13R | ID5,ID2,ID2,ID2 | chr16L|chr7R > chr4R |
| chr15L | ID5,ID2,ID3,ID3,ID3 | chr16L > chr14L |
| chr15R | ID2,ID3,ID2 | chr14L > chr15R |
| chr6L | ID4,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |
| chr14L | ID5,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |


## PD 33 — 125 reads
* paper flags 24 switching reads; we catch 23, miss 1, and flag 14 it calls N (88.0 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr14L | ID5,ID2,ID4,ID4,ID2,ID4,ID2,ID4 | chr2L|chr6L > chr14L > chr2L|chr6L > chr14L > chr2L|chr6L | chr2L|chr6L:ID4 > chr2L|chr6L:ID4 > self[2]:ID2 > chr2L|chr6L:ID4 > self[2]:ID2 > chr2L|chr6L:ID4 |
| chr11L | ID5,ID3 | chr14L | chr14L|chr16L:ID5 > chr14L[3]:ID3 |
| chr15R | ID2,ID2,ID4 | chr12R|chr14L|chr4R > chr2L|chr6L | chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 |
| chr12R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID4 | chr12R > chr2L|chr6L | self[6-7]:ID2,ID2,ID2,ID2(circ x2.0 moderate | alt chr12R[2] + chr12R[6-8]) > chr2L|chr6L:ID4 |
| chr10R | ID5,ID2,ID4 | chr14L > chr2L|chr6L | chr14L[1-2]:ID5,ID2 > chr2L|chr6L:ID4 |
| chr3R | ID5,ID2,ID4,ID2 | chr14L > chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr14L[1-2]:ID5,ID2 > chr2L|chr6L:ID4 > chr12R|chr14L|chr15R|chr16L|chr4R:ID2 |
| chr14R | ID6,ID4,ID2,ID4,ID2,ID4,ID2,ID4,ID2,ID4 | chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L | chr2L|chr6L:ID4 > chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 > chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 > chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 > chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 |
| chr8L | ID1,ID2,ID4,ID2 | chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 > chr12R|chr14L|chr15R|chr16L|chr4R:ID2 |
| chr7R | ID5,ID4,ID3,ID4 | chr2L|chr6L > chr14L > chr2L|chr6L | chr2L|chr6L:ID4 > chr14L[3]:ID3 > chr2L|chr6L:ID4 |
| chr11R | ID5,ID2,ID4 | chr14L > chr2L|chr6L | chr14L[1-2]:ID5,ID2 > chr2L|chr6L:ID4 |
| chr4L | ID5,ID2,ID4 | chr14L > chr2L|chr6L | chr14L[1-2]:ID5,ID2 > chr2L|chr6L:ID4 |
| chr13L | ID4,ID2,ID2 | chr2L|chr6L > chr16L|chr4R | chr2L|chr6L:ID4 > chr16L|chr4R:ID2,ID2 |
| chr1L | ID5,ID6 | chr14L|chr16L|chr7R > chr10L|chr14R|chr5L|chr9L | chr14L|chr16L|chr7R:ID5 > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr11L | ID6,ID2,ID2,ID2 | chr10L|chr14R|chr5L|chr9L > chr12R | chr10L|chr14R|chr5L|chr9L:ID6 > chr12R[5-7]:ID2,ID2,ID2 |
| chr9R | ID6,ID2,ID3 | chr10L|chr14R|chr5L|chr9L > chr14L | chr10L|chr14R|chr5L|chr9L:ID6 > chr14L[2-3]:ID2,ID3 |
| chr12L | ID1,ID2,ID4,ID2,ID4 | chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L | chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 > chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 |
| chr7R | ID5,ID2,ID4,ID2 | chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 > chr12R|chr14L|chr15R|chr16L|chr4R:ID2 |
| chr5L | ID6,ID2,ID4,ID6 | chr12R|chr14L|chr4R > chr2L|chr6L > chr5L | chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 > self[1]:ID6 |
| chr15L | ID5,ID3,ID4 | chr14L > chr2L|chr6L | chr14L|chr16L:ID5 > chr14L[3]:ID3 > chr2L|chr6L:ID4 |
| chr2L | ID4,ID2,ID2,ID2,ID2,ID2,ID6 | chr12R|chr4R > chr10L|chr14R|chr5L|chr9L | chr12R|chr4R:ID2,ID2,ID2,ID2,ID2(circ x1.7) > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr7R | ID5,ID2,ID4 | chr12R|chr14L|chr4R > chr2L|chr6L | chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 |
| chr5R | ID1,ID2,ID4 | chr12R|chr14L|chr4R > chr2L|chr6L | chr12R|chr14L|chr4R:ID2 > chr2L|chr6L:ID4 |
| chr1L | ID4,ID2 | chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr2L|chr6L:ID4 > chr12R|chr14L|chr15R|chr16L|chr4R:ID2 |
| chr11R | ID5,ID2,ID4 | chr14L > chr2L|chr6L | chr14L[1-2]:ID5,ID2 > chr2L|chr6L:ID4 |

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr2L | ID4,ID2,ID2,ID3 | chr4R > chr14L |
| chr2L | ID3,ID2,ID2 | chr14L > chr12R|chr4R |
| chr2R | ID1,ID2,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr5R | ID1,ID3,ID3,ID2,ID2 | chr14L > chr16L|chr4R |
| chr10R | ID1,ID2,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr4R | ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr12R > chr4R |
| chr4L | ID5,ID2,ID2,ID2,ID2 | chr14L|chr16L|chr7R > chr12R |
| chr6L | ID1,ID2,ID2,ID2 | chr12L|chr13L|chr16R|chr5R|chr8L|chr8R > chr4R |
| chr2L | ID3,ID2,ID2 | chr14L > chr12R |
| chr7L | ID5,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |
| chr6L | ID1,ID2,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr16R | ID1,ID1,ID2,ID2,ID2,ID2,ID2,ID2 | chr16R > chr12R|chr4R |
| chr2R | ID3,ID2,ID2,ID2 | chr14L > chr4R |
| chr15L | ID5,ID2,ID3,ID3,ID3 | chr16L > chr14L |

