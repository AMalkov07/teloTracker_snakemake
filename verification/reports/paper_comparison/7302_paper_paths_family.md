# Summary

| PD | reads | paper_Y | caught | missed | extra | agreement_pct |
|---|---|---|---|---|---|---|
| 28.0 | 103.0 | 9.0 | 4.0 | 5.0 | 3.0 | 92.2 |
| 33.0 | 86.0 | 7.0 | 5.0 | 2.0 | 3.0 | 94.2 |


# Our path parser applied to the paper's own Y' calls — strain 7302

reference: `7302_day0_with_selection`, curated 7302 library; orientation inferred from the data.

## PD 28 — 103 reads
* paper flags 9 switching reads; we catch 4, miss 5, and flag 3 it calls N (92.2 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr13L | ID7,ID2,ID8,ID2,ID8,ID2,ID8 | chr13L > chr12R|chr14L|chr4R > chr13L > chr12R|chr4R > chr13L | self[4]:ID8 > chr12R|chr14L|chr4R:ID2 > self[4]:ID8 > chr12R|chr4R:ID2 > self[4]:ID8 |
| chr13R | ID4,ID2 | chr2L|chr6L > chr12R|chr13L|chr14L|chr15R|chr4R | chr2L|chr6L:ID4 > chr12R|chr13L|chr14L|chr15R|chr4R:ID2 |
| chr14R | ID5,ID3,ID3,ID3 | chr14L | chr14L?[1]:ID5 > chr14L[3-5]:ID3,ID3,ID3 |
| chr13R | ID5,ID7 | chr14L|chr16L|chr7R > chr13L | chr14L|chr16L|chr7R:ID5 > chr13L[1]:ID7 |
| chr7R | ID5,ID5,ID3,ID5,ID5,ID3 | chr14L | chr14L?[1]:ID5 > chr14L[3]:ID3 > chr14L?[1]:ID5 > chr14L?[1]:ID5 > chr14L[3]:ID3 |
| chr13L | ID7,ID8,ID7,ID2,ID8,ID7,ID2,ID8,ID7,ID2 | chr13L | self[4]:ID8 > self[1-2]:ID7,ID2 > self[4]:ID8 > self[1-2]:ID7,ID2 > self[4]:ID8 > self[1-2]:ID7,ID2 |
| chr15L | ID6,ID3,ID3 | chr10L|chr14R|chr5L|chr9L > chr14L | chr10L|chr14R|chr5L|chr9L:ID6 > chr14L[3-4]:ID3,ID3 |
| chr16R | ID1,ID3,ID3,ID2 | chr14L | chr14L[3-4]:ID3,ID3 > chr12R|chr13L|chr14L|chr15R|chr4R:ID2 |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID2 | chr14L | self[3]:ID3 > self[2]:ID2 |

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr13L | ID7,ID8,ID7,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr13L > chr4R |
| chr12L | ID3,ID2,ID2,ID2 | chr14L > chr12R|chr4R |
| chr16R | ID1,ID2,ID3 | chr13L > chr14L |


## PD 33 — 86 reads
* paper flags 7 switching reads; we catch 5, miss 2, and flag 3 it calls N (94.2 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr3R | ID6,ID2,ID3,ID3 | chr10L|chr14R|chr5L|chr9L > chr14L | chr10L|chr14R|chr5L|chr9L:ID6 > chr14L[2-4]:ID2,ID3,ID3 |
| chr13L | ID7,ID8,ID7,ID8,ID8,ID7,ID2,ID2,ID2,ID2,ID2 | chr13L | self[3-4]:ID8,ID7,ID8(circ x1.5 moderate | alt chr13L[4] + chr13L[3-4]) > self[4]:ID8 > self[3]:ID7 > self[2]:ID2,ID2,ID2,ID2,ID2(circ x5.0 strong) |
| chr16R | ID3,ID7,ID8 | chr14L > chr13L | chr14L[3]:ID3 > chr13L[3-4]:ID7,ID8 |
| chr15R | ID2,ID2,ID8,ID8 | chr13L | chr13L?[2]:ID2 > chr13L[4]:ID8 > chr13L[4]:ID8 |
| chr16R | ID3,ID7,ID8 | chr14L > chr13L | chr14L[3]:ID3 > chr13L[1]:ID7 > chr13L[4]:ID8 |
| chr1R | ID5,ID7,ID3 | chr14L|chr16L|chr7R > chr13L > chr14L | chr14L|chr16L|chr7R:ID5 > chr13L[1]:ID7 > chr14L[3]:ID3 |
| chr3R | ID5,ID8,ID7,ID8 | chr14L > chr13L | chr14L?[1]:ID5 > chr13L[3-4]:ID8,ID7,ID8(circ x1.5 moderate | alt chr13L[4] + chr13L[3-4]) |

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr2L | ID3,ID2,ID2 | chr14L > chr12R|chr4R |
| chr12R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID3,ID3 | chr13L > chr14L |
| chr2R | ID3,ID2,ID2,ID2,ID2 | chr14L > chr12R |

