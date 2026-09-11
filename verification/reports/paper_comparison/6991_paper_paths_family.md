# Summary

| PD | reads | paper_Y | caught | missed | extra | agreement_pct |
|---|---|---|---|---|---|---|
| 28.0 | 240.0 | 42.0 | 24.0 | 18.0 | 11.0 | 87.9 |
| 33.0 | 244.0 | 59.0 | 36.0 | 23.0 | 11.0 | 86.1 |


# Our path parser applied to the paper's own Y' calls — strain 6991

reference: `6991_day0_with_selection_repeat`, curated 6991 library; orientation inferred from the data.

## PD 28 — 240 reads
* paper flags 42 switching reads; we catch 24, miss 18, and flag 11 it calls N (87.9 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr7R | ID5,ID2,ID2 | chr4R | chr4R[1-2]:ID2,ID2 |
| chr13R | ID5,ID2,ID2 | chr14L | chr14L[1-2]:ID5,ID2 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr2R | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L[3]:ID3 |
| chr15R | ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2,ID2 |
| chr11R | ID2,ID3 | chr14L | chr14L[2-3]:ID2,ID3 |
| chr12L | ID1,ID2,ID2,ID6 | chr12R|chr14L|chr15R|chr4R > chr10L|chr14R|chr5L|chr9L | chr12R|chr14L|chr15R|chr4R:ID2 > chr12R|chr14L|chr15R|chr4R:ID2 > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr12R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID3 | chr12R > chr14L | self[2-3]:ID2,ID2 > chr14L[3]:ID3 |
| chr2L | ID4,ID3,ID4 | chr14L > chr2L | chr14L[3]:ID3 > self[1]:ID4 |
| chr7L | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1 > chr14L[3]:ID3 |
| chr15R | ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2,ID2 |
| chr11L | ID6,ID3 | chr10L|chr14R|chr5L|chr9L > chr14L | chr10L|chr14R|chr5L|chr9L:ID6 > chr14L[3]:ID3 |
| chr5R | ID1,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2,ID2 |
| chr10R | ID3,ID4 | chr14L > chr2L|chr6L | chr14L[3]:ID3 > chr2L|chr6L:ID4 |
| chr11L | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L[3]:ID3 |
| chr15R | ID2,ID2,ID2,ID6 | chr4R > chr10L|chr14R|chr5L|chr9L | chr4R[1-2]:ID2,ID2 > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr14L | ID5,ID2,ID3,ID3,ID5,ID3,ID3,ID3 | chr14L | self[1]:ID5 > self[3-5]:ID3,ID3,ID3 |
| chr3R | ID6,ID2,ID2 | chr10L|chr14R|chr5L|chr9L > chr12R|chr4R | chr10L|chr14R|chr5L|chr9L:ID6 > chr12R|chr4R:ID2,ID2 |
| chr3R | ID5,ID2 | chr14L | chr14L[1-2]:ID5,ID2 |
| chr8R | ID1,ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2,ID2,ID2 |
| chr12R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr12R | self[2]:ID2 > self[2-4]:ID2,ID2,ID2,ID2(circ x1.3 moderate | alt chr12R[2] + chr12R[2-4]) |
| chr10L | ID5,ID2,ID2,ID2,ID2 | chr14L > chr4R | chr14L[1-2]:ID5,ID2 > chr4R[1-3]:ID2,ID2,ID2 |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID2,ID3,ID3,ID3,ID3,ID3,ID2 | chr14L | self[2-3]:ID2,ID3 > self[3]:ID3,ID3,ID3,ID3(circ x4.0 moderate | alt chr14L[3-5] + chr14L[3]) > self[2]:ID2 |
| chr8L | ID1,ID2,ID6,ID2 | chr12R|chr14L|chr4R > chr10L|chr14R|chr5L|chr9L > chr12R|chr14L|chr15R|chr4R | chr12R|chr14L|chr4R:ID2 > chr10L|chr14R|chr5L|chr9L:ID6 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr13R | ID5,ID3,ID3,ID3,ID3,ID2,ID3,ID3,ID3,ID3 | chr14L | chr14L?[1]:ID5 > chr14L[3-5]:ID3,ID3,ID3 > chr14L[2-3]:ID3,ID2,ID3(circ x1.5 moderate | alt chr14L[3] + chr14L[2-3]) > chr14L[3-5]:ID3,ID3,ID3 |
| chr9R | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1 > chr14L[3]:ID3 |
| chr11L | ID5,ID3 | chr14L | chr14L?[1]:ID5 > chr14L[3]:ID3 |
| chr2R | ID3,ID4 | chr14L > chr2L|chr6L | chr14L[3]:ID3 > chr2L|chr6L:ID4 |
| chr16L | ID5,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr14L|chr15R|chr4R:ID2 > chr12R|chr14L|chr15R|chr4R:ID2 > chr12R|chr4R:ID2,ID2,ID2(circ x1.5) > chr12R|chr14L|chr15R|chr4R:ID2 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr7R | ID5,ID5,ID2 | chr7R > chr12R|chr14L|chr15R|chr4R | self[1]:ID5 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr2R | ID1,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L[3]:ID3,ID3,ID3,ID3,ID3,ID3,ID3(circ x7.0 strong) |
| chr2L | ID3,ID1,ID1 | chr14L > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | chr14L[3]:ID3 > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID2,ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2,ID2,ID2,ID2(circ x2.0) |
| chr11L | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L[3]:ID3 |
| chr15L | ID4,ID2,ID3 | chr2L|chr6L > chr14L | chr2L|chr6L:ID4 > chr14L[2-3]:ID2,ID3 |
| chr7R | ID5,ID5,ID2,ID2 | chr14L | chr14L[1-2]:ID5,ID2 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr11R | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1 > chr14L[3]:ID3 |
| chr15L | ID5,ID3 | chr14L | chr14L?[1]:ID5 > chr14L[3]:ID3 |
| chr13L | ID1,ID3,ID1 | chr14L > chr13L | chr14L[3]:ID3 > self[1]:ID1 |
| chr3R | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1 > chr14L[3]:ID3 |
| chr13R | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1 > chr14L[3]:ID3 |

_(2 more rows in the TSV)_

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr4L | ID5,ID2,ID2,ID2,ID2 | chr14L|chr16L|chr7R > chr12R |
| chr8R | ID3,ID2,ID2,ID2 | chr14L > chr4R |
| chr13R | ID3,ID2,ID2 | chr14L > chr12R|chr4R |
| chr11R | ID1,ID2,ID2,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr4R |
| chr3R | ID1,ID3,ID3,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L |
| chr13L | ID1,ID1,ID3,ID2 | chr12R > chr14L |
| chr5R | ID1,ID2,ID3,ID3,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |
| chr13L | ID3,ID2,ID2,ID2 | chr14L > chr12R|chr4R |
| chr5R | ID1,ID1,ID1 | chr12R > chr5R |
| chr5R | ID1,ID2,ID2,ID1 | chr12R|chr4R > chr5R |
| chr11R | ID3,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr4R |


## PD 33 — 244 reads
* paper flags 59 switching reads; we catch 36, miss 23, and flag 11 it calls N (86.1 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr10R | ID1,ID2,ID2 | chr12R | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr12R|chr4R:ID2,ID2 |
| chr1L | ID2,ID5 | chr14L | chr12R|chr14L|chr4R:ID2 > chr14L|chr16L|chr7R:ID5 |
| chr8R | ID1,ID3,ID1,ID1,ID1,ID1 | chr14L > chr12R > chr8R | chr14L[3]:ID3 > chr12R?[1]:ID1 > self[1]:ID1,ID1,ID1(circ x3.0 strong) |
| chr5R | ID1,ID3,ID1,ID3,ID1,ID3,ID1,ID3,ID1,ID3,ID1,ID3,ID1 | chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R | chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 |
| chr7L | ID1,ID3,ID3 | chr12R > chr14L | chr12R?[1]:ID1 > chr14L[3-4]:ID3,ID3 |
| chr13L | ID1,ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1 > chr14L[3]:ID3 |
| chr13R | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L[3]:ID3 |
| chr3R | ID5,ID3 | chr14L | chr14L?[1]:ID5 > chr14L[3]:ID3 |
| chr13R | ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2,ID2,ID2 |
| chr4L | ID2,ID6 | chr12R|chr14L|chr15R|chr4R > chr10L|chr14R|chr5L|chr9L | chr12R|chr14L|chr15R|chr4R:ID2 > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr11L | ID6,ID2,ID2,ID2,ID2,ID2,ID2 | chr10L|chr14R|chr5L|chr9L > chr12R | chr10L|chr14R|chr5L|chr9L:ID6 > chr12R[2-7]:ID2,ID2,ID2,ID2,ID2,ID2 |
| chr1R | ID1,ID4 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr2L|chr6L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr2L|chr6L:ID4 |
| chr2R | ID1,ID5 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L|chr16L|chr7R | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L|chr16L|chr7R:ID5 |
| chr4L | ID3,ID2 | chr14L | chr14L[3]:ID3 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr4L | ID5,ID2 | chr14L | chr14L[1-2]:ID5,ID2 |
| chr1L | ID6,ID3 | chr10L|chr14R|chr5L|chr9L > chr14L | chr10L|chr14R|chr5L|chr9L:ID6 > chr14L[3]:ID3 |
| chr14R | ID6,ID2,ID6 | chr4R > chr14R | chr4R?[1]:ID2 > self[1]:ID6 |
| chr11R | ID6,ID2,ID3,ID3 | chr10L|chr14R|chr5L|chr9L > chr14L | chr10L|chr14R|chr5L|chr9L:ID6 > chr14L[2-4]:ID2,ID3,ID3 |
| chr10R | ID1,ID1,ID1,ID1,ID1,ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1,ID1,ID1,ID1,ID1,ID1(circ x6.0 moderate) > chr14L[3]:ID3 |
| chr13R | ID5,ID2 | chr14L | chr14L[1-2]:ID5,ID2 |
| chr10R | ID6,ID2 | chr10L|chr14R|chr5L|chr9L > chr12R|chr14L|chr15R|chr4R | chr10L|chr14R|chr5L|chr9L:ID6 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr7L | ID1,ID2,ID3,ID3 | chr12R > chr14L | chr12R[1-2]:ID1,ID2 > chr14L[3-4]:ID3,ID3 |
| chr5R | ID1,ID1,ID6 | chr5R > chr10L|chr14R|chr5L|chr9L | self[1]:ID1 > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr15L | ID6,ID1 | chr10L|chr14R|chr5L|chr9L > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | chr10L|chr14R|chr5L|chr9L:ID6 > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 |
| chr7L | ID1,ID6 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr10L|chr14R|chr5L|chr9L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr7L | ID5,ID3,ID6,ID6,ID3,ID6,ID6 | chr14L > chr10L|chr14R|chr5L|chr9L > chr14L > chr10L|chr14R|chr5L|chr9L | chr14L?[1]:ID5 > chr14L[3]:ID3 > chr10L|chr14R|chr5L|chr9L:ID6 > chr10L|chr14R|chr5L|chr9L:ID6 > chr14L[3]:ID3 > chr10L|chr14R|chr5L|chr9L:ID6 > chr10L|chr14R|chr5L|chr9L:ID6 |
| chr9R | ID4,ID2 | chr2L|chr6L > chr12R|chr14L|chr15R|chr4R | chr2L|chr6L:ID4 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr7L | ID2,ID2,ID2 | chr4R | chr4R[1-3]:ID2,ID2,ID2 |
| chr3L | ID5,ID2,ID2,ID2 | chr14L > chr4R | chr14L[1-2]:ID5,ID2 > chr4R[1-2]:ID2,ID2 |
| chr4L | ID1,ID3,ID3,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L[3-5]:ID3,ID3,ID3 |
| chr5R | ID1,ID3,ID1,ID3,ID1,ID3,ID5,ID1,ID3,ID1,ID3,ID1,ID3,ID1 | chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R | chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > chr14L|chr16L|chr7R:ID5 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 > chr14L[3]:ID3 > self[1]:ID1 |
| chr5R | ID1,ID2,ID6,ID2 | chr12R|chr14L|chr15R|chr4R > chr10L|chr14R|chr5L|chr9L > chr12R|chr14L|chr15R|chr4R | chr12R|chr14L|chr15R|chr4R:ID2 > chr10L|chr14R|chr5L|chr9L:ID6 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr13R | ID5,ID2,ID3,ID3,ID3 | chr14L | chr14L[1-2]:ID5,ID2 > chr14L[3-5]:ID3,ID3,ID3 |
| chr9L | ID6,ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2,ID2,ID2 |
| chr3L | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr14L[3]:ID3 |
| chr3R | ID1,ID2 | chr12R | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1 > chr12R|chr14L|chr15R|chr4R:ID2 |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID2,ID2,ID3,ID3,ID3,ID3,ID3 | chr12R|chr4R > chr14L | chr12R|chr4R:ID2,ID2 > self[3]:ID3,ID3,ID3,ID3,ID3(circ x5.0 moderate | alt chr14L[3-5] + chr14L[3-4]) |
| chr8R | ID1,ID3,ID1 | chr14L > chr8R | chr14L[3]:ID3 > self[1]:ID1 |
| chr10R | ID5,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr14L | chr14L?[1]:ID5 > chr14L[3]:ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3(circ x9.0 strong) |
| chr15R | ID2,ID2,ID2,ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2,ID2,ID2,ID2,ID2(circ x1.2) |

_(19 more rows in the TSV)_

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr12L | ID1,ID2,ID2,ID2,ID3 | chr12R|chr4R > chr14L |
| chr12L | ID3,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |
| chr7R | ID5,ID5,ID5 | chr14L > chr7R |
| chr16R | ID1,ID2,ID2,ID1 | chr4R > chr16R |
| chr7R | ID5,ID5,ID5 | chr14L > chr7R |
| chr11L | ID1,ID2,ID2 | chr12L|chr13L|chr16R|chr5R|chr8L|chr8R > chr4R |
| chr7R | ID5,ID5,ID5 | chr14L > chr7R |
| chr5R | ID1,ID2,ID1 | chr12R|chr4R > chr5R |
| chr13R | ID1,ID2,ID2,ID2 | chr12L|chr13L|chr16R|chr5R|chr8L|chr8R > chr12R |
| chr14R | ID6,ID6,ID2 | chr14R > chr12R|chr14L|chr15R|chr4R |
| chr5R | ID1,ID1,ID1,ID1 | chr5R > chr12R > chr5R |

