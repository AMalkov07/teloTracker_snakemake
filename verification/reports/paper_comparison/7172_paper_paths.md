# Summary

| PD | reads | paper_Y | caught | missed | extra | agreement_pct |
|---|---|---|---|---|---|---|
| 28.0 | 60.0 | 9.0 | 9.0 | 0.0 | 7.0 | 88.3 |
| 33.0 | 125.0 | 24.0 | 24.0 | 0.0 | 15.0 | 88.0 |


# Our path parser applied to the paper's own Y' calls — strain 7172

reference: `7172_day0_with_selection`, curated 7172 library; orientation inferred from the data.

**Orientation heuristic: 179/179 reads (100.0 %) match the true telo_side.**

## PD 28 — 60 reads
* paper flags 9 switching reads; we catch 9, miss 0, and flag 7 it calls N (88.3 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr6L | ID3,ID3,ID5 | chr14L > chr16L|chr7R | chr14L[3]:ID3_Orange > chr14L[3]:ID3_Orange > chr16L|chr7R:ID5_Blue-Dark |
| chr1R | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14L[3]:ID3_Orange |
| chr15L | ID5,ID2,ID3,ID4 | chr16L|chr7R > chr14L > chr2L|chr6L | chr16L|chr7R:ID5_Blue-Dark > chr14L[2-3]:ID2_Red,ID3_Orange > chr2L|chr6L:ID4_Green-Light |
| chr11R | ID4,ID3,ID3,ID3,ID3,ID3 | chr2L|chr6L > chr14L | chr2L|chr6L:ID4_Green-Light > chr14L[3]:ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange(circ x5.0 strong) |
| chr11L | ID6,ID3,ID3,ID3,ID3 | chr14R > chr14L | chr14R[1]:ID6_Purple-Neutral > chr14L[3]:ID3_Orange > chr14L[3]:ID3_Orange,ID3_Orange,ID3_Orange(circ x3.0 moderate) |
| chr16R | ID1,ID2,ID2,ID6 | chr12R|chr4R > chr10L|chr9L | chr12R|chr4R:ID2_Red,ID2_Red > chr10L|chr9L:ID6_Purple-Dark |
| chr9L | ID6,ID5,ID3,ID3,ID5 | chr16L > chr14L | chr16L?[1]:ID5_Blue-Dark > chr14L[3]:ID3_Orange > chr14L[3]:ID3_Orange > chr14L[1]:ID5_Blue-Light |
| chr1L | ID5,ID5,ID5 | chr16L > chr14L | chr16L?[1]:ID5_Blue-Dark > chr14L[1]:ID5_Blue-Light > chr14L[1]:ID5_Blue-Light |
| chr13R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID6 | chr12R > chr10L|chr9L | chr12R[1-7]:ID1_Gray,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red > chr12R|chr4R:ID2_Red,ID2_Red > chr10L|chr9L:ID6_Purple-Dark |

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr16L | ID5,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr16L > chr12R > chr16L |
| chr13R | ID5,ID2,ID2,ID2 | chr14L > chr12R|chr4R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID2 | chr16L > chr14L |
| chr15L | ID5,ID2,ID3 | chr16L|chr7R > chr14L |
| chr15R | ID2,ID3,ID2 | chr14L > chr15R |
| chr6L | ID4,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |
| chr14L | ID5,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |


## PD 33 — 125 reads
* paper flags 24 switching reads; we catch 24, miss 0, and flag 15 it calls N (88.0 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr14L | ID5,ID2,ID4,ID4,ID2,ID4,ID2,ID4 | chr16L > chr2L|chr6L > chr14L > chr2L|chr6L > chr14L > chr2L|chr6L | chr16L[1-2]:ID5_Blue-Dark,ID2_Red > chr2L|chr6L:ID4_Green-Light > chr2L|chr6L:ID4_Green-Light > self[2]:ID2_Red > chr2L|chr6L:ID4_Green-Light > self[2]:ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr11L | ID5,ID3 | chr16L > chr14L | chr16L?[1]:ID5_Blue-Dark > chr14L[3]:ID3_Orange |
| chr15R | ID2,ID2,ID4 | chr12R|chr14L|chr4R > chr2L|chr6L | chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr12R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID4 | chr12R > chr2L|chr6L | self[6-7]:ID2_Red,ID2_Red,ID2_Red,ID2_Red(circ x2.0 moderate | alt chr12R[2] + chr12R[6-8]) > chr2L|chr6L:ID4_Green-Light |
| chr10R | ID5,ID2,ID4 | chr16L > chr2L|chr6L | chr16L[1-2]:ID5_Blue-Dark,ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr3R | ID5,ID2,ID4,ID2 | chr16L > chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr16L[1-2]:ID5_Blue-Dark,ID2_Red > chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr15R|chr16L|chr4R:ID2_Red |
| chr14R | ID6,ID4,ID2,ID4,ID2,ID4,ID2,ID4,ID2,ID4 | chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L | chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr8L | ID1,ID2,ID4,ID2 | chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr15R|chr16L|chr4R:ID2_Red |
| chr7R | ID5,ID4,ID3,ID4 | chr2L|chr6L > chr14L > chr2L|chr6L | chr2L|chr6L:ID4_Green-Light > chr14L[3]:ID3_Orange > chr2L|chr6L:ID4_Green-Light |
| chr11R | ID5,ID2,ID4 | chr16L > chr2L|chr6L | chr16L[1-2]:ID5_Blue-Dark,ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr4L | ID5,ID2,ID4 | chr14L > chr2L|chr6L | chr14L[1-2]:ID5_Blue-Light,ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr13L | ID4,ID2,ID2 | chr2L|chr6L > chr16L|chr4R | chr2L|chr6L:ID4_Green-Light > chr16L|chr4R:ID2_Red,ID2_Red |
| chr1L | ID5,ID6 | chr16L|chr7R > chr10L|chr9L | chr16L|chr7R:ID5_Blue-Dark > chr10L|chr9L:ID6_Purple-Dark |
| chr11L | ID6,ID2,ID2,ID2 | chr10L|chr9L > chr12R | chr10L|chr9L:ID6_Purple-Dark > chr12R[5-7]:ID2_Red,ID2_Red,ID2_Red |
| chr9R | ID6,ID2,ID3 | chr14R > chr14L | chr14R[1]:ID6_Purple-Neutral > chr14L[2-3]:ID2_Red,ID3_Orange |
| chr12L | ID1,ID2,ID4,ID2,ID4 | chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L | chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr7R | ID5,ID2,ID4,ID2 | chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr15R|chr16L|chr4R:ID2_Red |
| chr5L | ID6,ID2,ID4,ID6 | chr12R|chr14L|chr4R > chr2L|chr6L > chr5L | chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light > self[1]:ID6_Purple-Light |
| chr15L | ID5,ID3,ID4 | chr16L > chr14L > chr2L|chr6L | chr16L?[1]:ID5_Blue-Dark > chr14L[3]:ID3_Orange > chr2L|chr6L:ID4_Green-Light |
| chr2L | ID4,ID2,ID2,ID2,ID2,ID2,ID6 | chr12R|chr4R > chr10L|chr9L | chr12R|chr4R:ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red(circ x1.7) > chr10L|chr9L:ID6_Purple-Dark |
| chr7R | ID5,ID2,ID4 | chr12R|chr14L|chr4R > chr2L|chr6L | chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr5R | ID1,ID2,ID4 | chr12R|chr14L|chr4R > chr2L|chr6L | chr12R|chr14L|chr4R:ID2_Red > chr2L|chr6L:ID4_Green-Light |
| chr1L | ID4,ID2 | chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr2L|chr6L:ID4_Green-Light > chr12R|chr14L|chr15R|chr16L|chr4R:ID2_Red |
| chr11R | ID5,ID2,ID4 | chr16L > chr2L|chr6L | chr16L[1-2]:ID5_Blue-Dark,ID2_Red > chr2L|chr6L:ID4_Green-Light |

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr2L | ID4,ID2,ID2,ID3 | chr4R > chr14L |
| chr2L | ID3,ID2,ID2 | chr14L > chr12R|chr4R |
| chr2R | ID1,ID2,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr5R | ID1,ID3,ID3,ID2,ID2 | chr14L > chr16L|chr4R |
| chr10R | ID1,ID2,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr11L | ID5,ID2,ID3 | chr16L|chr7R > chr14L |
| chr4R | ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr12R > chr4R |
| chr4L | ID5,ID2,ID2,ID2,ID2 | chr16L|chr7R > chr12R |
| chr13R | ID5,ID2,ID2 | chr7R > chr12R|chr4R |
| chr6L | ID1,ID2,ID2,ID2 | chr12L|chr13L|chr16R|chr5R|chr8L|chr8R > chr4R |
| chr2L | ID3,ID2,ID2 | chr14L > chr12R |
| chr7L | ID5,ID2,ID2,ID2,ID2,ID2 | chr7R > chr12R |
| chr6L | ID1,ID2,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr16R | ID1,ID1,ID2,ID2,ID2,ID2,ID2,ID2 | chr16R > chr12R|chr4R |
| chr2R | ID3,ID2,ID2,ID2 | chr14L > chr4R |

