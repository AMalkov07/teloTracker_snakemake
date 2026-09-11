# Summary

| PD | reads | paper_Y | caught | missed | extra | agreement_pct |
|---|---|---|---|---|---|---|
| 28.0 | 240.0 | 42.0 | 41.0 | 1.0 | 25.0 | 89.2 |
| 33.0 | 244.0 | 59.0 | 59.0 | 0.0 | 17.0 | 93.0 |


# Our path parser applied to the paper's own Y' calls — strain 6991

reference: `6991_day0_with_selection_repeat`, curated 6991 library; orientation inferred from the data.

## PD 28 — 240 reads
* paper flags 42 switching reads; we catch 41, miss 1, and flag 25 it calls N (89.2 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr7R | ID5,ID2,ID2 | chr4R > chr15R | chr4R?[1]:ID2_Red-Light > chr15R[1]:ID2_Red-Dark |
| chr13R | ID5,ID2,ID2 | chr16L|chr7R > chr15R | chr16L|chr7R:ID5_Blue-Dark > chr15R[1]:ID2_Red-Dark > chr15R[1]:ID2_Red-Dark |
| chr2R | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14L[3]:ID3_Orange |
| chr15R | ID2,ID2,ID2 | chr12R|chr4R | chr12R|chr4R:ID2_Red-Light,ID2_Red-Light |
| chr11R | ID2,ID3 | chr15R > chr14L | chr15R[1]:ID2_Red-Dark > chr14L[3]:ID3_Orange |
| chr12L | ID1,ID2,ID2,ID6 | chr12R|chr14L|chr4R > chr15R > chr14R | chr12R|chr14L|chr4R:ID2_Red-Light > chr15R[1]:ID2_Red-Dark > chr14R[1]:ID6_Purple-Neutral |
| chr12R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID3 | chr15R > chr14L | chr15R[1]:ID2_Red-Dark > chr14L[2-3]:ID2_Red-Light,ID3_Orange |
| chr2L | ID4,ID3,ID4 | chr14L > chr2L | chr14L[3]:ID3_Orange > self[1]:ID4_Green-Light |
| chr7L | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1_Gray > chr14L[3]:ID3_Orange |
| chr15R | ID2,ID2,ID2 | chr12R|chr14L|chr4R > chr15R | chr12R|chr14L|chr4R:ID2_Red-Light > self[1]:ID2_Red-Dark |
| chr11L | ID6,ID3 | chr14R > chr14L | chr14R[1]:ID6_Purple-Neutral > chr14L[3]:ID3_Orange |
| chr5R | ID1,ID2,ID2 | chr15R > chr12R|chr14L|chr4R | chr15R[1]:ID2_Red-Dark > chr12R|chr14L|chr4R:ID2_Red-Light |
| chr10R | ID3,ID4 | chr14L > chr2L|chr6L | chr14L[3]:ID3_Orange > chr2L|chr6L:ID4_Green-Light |
| chr11L | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14L[3]:ID3_Orange |
| chr15R | ID2,ID2,ID2,ID6 | chr4R > chr15R > chr10L|chr9L | chr4R?[1]:ID2_Red-Light > self[1]:ID2_Red-Dark > chr10L|chr9L:ID6_Purple-Dark |
| chr14L | ID5,ID2,ID3,ID3,ID5,ID3,ID3,ID3 | chr16L|chr7R > chr14L | chr16L|chr7R:ID5_Blue-Dark > self[3-5]:ID3_Orange,ID3_Orange,ID3_Orange |
| chr3R | ID6,ID2,ID2 | chr14R > chr12R|chr4R | chr14R[1]:ID6_Purple-Neutral > chr12R|chr4R:ID2_Red-Light,ID2_Red-Light |
| chr3R | ID5,ID2 | chr16L|chr7R > chr12R|chr14L|chr4R | chr16L|chr7R:ID5_Blue-Dark > chr12R|chr14L|chr4R:ID2_Red-Light |
| chr8R | ID1,ID2,ID2,ID2 | chr12R|chr4R > chr15R | chr12R|chr4R:ID2_Red-Light,ID2_Red-Light > chr15R[1]:ID2_Red-Dark |
| chr12R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr15R > chr12R | chr15R[1]:ID2_Red-Dark > self[2]:ID2_Red-Light > self[2]:ID2_Red-Light > self[2-4]:ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light(circ x1.3 moderate | alt chr12R[2] + chr12R[2-4]) |
| chr10L | ID5,ID2,ID2,ID2,ID2 | chr16L|chr7R > chr4R | chr16L|chr7R:ID5_Blue-Dark > chr4R[1-4]:ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID2,ID3,ID3,ID3,ID3,ID3,ID2 | chr14L > chr15R | self[2-3]:ID2_Red-Light,ID3_Orange > self[3]:ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange(circ x4.0 moderate | alt chr14L[3-5] + chr14L[3]) > chr15R[1]:ID2_Red-Dark |
| chr8L | ID1,ID2,ID6,ID2 | chr12R|chr14L|chr4R > chr14R > chr12R|chr14L|chr4R | chr12R|chr14L|chr4R:ID2_Red-Light > chr14R[1]:ID6_Purple-Neutral > chr12R|chr14L|chr4R:ID2_Red-Light |
| chr13R | ID5,ID3,ID3,ID3,ID3,ID2,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L | chr16L|chr7R:ID5_Blue-Dark > chr14L[3-5]:ID3_Orange,ID3_Orange,ID3_Orange > chr14L[2-3]:ID3_Orange,ID2_Red-Light,ID3_Orange(circ x1.5 moderate | alt chr14L[3] + chr14L[2-3]) > chr14L[3-5]:ID3_Orange,ID3_Orange,ID3_Orange |
| chr9R | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1_Gray > chr14L[3]:ID3_Orange |
| chr11L | ID5,ID3 | chr16L|chr7R > chr14L | chr16L|chr7R:ID5_Blue-Dark > chr14L[3]:ID3_Orange |
| chr2R | ID3,ID4 | chr14L > chr2L|chr6L | chr14L[3]:ID3_Orange > chr2L|chr6L:ID4_Green-Light |
| chr16L | ID5,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr15R > chr12R|chr14L|chr4R > chr15R | chr15R[1]:ID2_Red-Dark > chr15R[1]:ID2_Red-Dark > chr12R|chr14L|chr4R:ID2_Red-Light > chr15R[1]:ID2_Red-Dark > chr15R[1]:ID2_Red-Dark > chr15R[1]:ID2_Red-Dark > chr15R[1]:ID2_Red-Dark |
| chr7R | ID5,ID5,ID2 | chr7R > chr12R|chr14L|chr4R | self[1]:ID5_Blue-Dark > chr12R|chr14L|chr4R:ID2_Red-Light |
| chr2R | ID1,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14L[3]:ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange(circ x7.0 strong) |
| chr2L | ID3,ID1,ID1 | chr14L > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | chr14L[3]:ID3_Orange > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID2,ID2,ID2,ID2 | chr16L|chr7R > chr14L > chr15R > chr14L | chr16L|chr7R:ID5_Blue-Dark > self[2-3]:ID2_Red-Light,ID3_Orange > self[3-4]:ID3_Orange,ID3_Orange > self[2]:ID2_Red-Light > chr15R[1]:ID2_Red-Dark > chr15R[1]:ID2_Red-Dark > self[2]:ID2_Red-Light |
| chr11L | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14L[3]:ID3_Orange |
| chr15L | ID4,ID2,ID3 | chr2L|chr6L > chr14L | chr2L|chr6L:ID4_Green-Light > chr14L[2-3]:ID2_Red-Light,ID3_Orange |
| chr7R | ID5,ID5,ID2,ID2 | chr7R > chr12R|chr4R | self[1]:ID5_Blue-Dark > chr12R|chr4R:ID2_Red-Light,ID2_Red-Light |
| chr11R | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1_Gray > chr14L[3]:ID3_Orange |
| chr15L | ID5,ID3 | chr16L|chr7R > chr14L | chr16L|chr7R:ID5_Blue-Dark > chr14L[3]:ID3_Orange |
| chr13L | ID1,ID3,ID1 | chr14L > chr13L | chr14L[3]:ID3_Orange > self[1]:ID1_Gray |
| chr3R | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1_Gray > chr14L[3]:ID3_Orange |
| chr13R | ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1_Gray > chr14L[3]:ID3_Orange |

_(2 more rows in the TSV)_

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr14L | ID5,ID2,ID3,ID2,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr7L | ID5,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID2,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr4L | ID5,ID2,ID2,ID2,ID2 | chr16L|chr7R > chr12R |
| chr8R | ID3,ID2,ID2,ID2 | chr14L > chr4R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr13R | ID3,ID2,ID2 | chr14L > chr12R|chr4R |
| chr11R | ID1,ID2,ID2,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr4R |
| chr3R | ID1,ID3,ID3,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr13L | ID1,ID1,ID3,ID2 | chr12R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr5R | ID1,ID2,ID3,ID3,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID2,ID3,ID3 | chr16L|chr7R > chr14L |
| chr13L | ID3,ID2,ID2,ID2 | chr14L > chr12R|chr4R |
| chr5R | ID1,ID1,ID1 | chr12R > chr5R |
| chr5R | ID1,ID2,ID2,ID1 | chr12R|chr4R > chr5R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr11R | ID3,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr4R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |


## PD 33 — 244 reads
* paper flags 59 switching reads; we catch 59, miss 0, and flag 17 it calls N (93.0 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr10R | ID1,ID2,ID2 | chr12R > chr15R | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr12R|chr14L|chr4R:ID2_Red-Light > chr15R[1]:ID2_Red-Dark |
| chr1L | ID2,ID5 | chr12R|chr14L|chr4R > chr16L|chr7R | chr12R|chr14L|chr4R:ID2_Red-Light > chr16L|chr7R:ID5_Blue-Dark |
| chr8R | ID1,ID3,ID1,ID1,ID1,ID1 | chr14L > chr12R > chr8R | chr14L[3]:ID3_Orange > chr12R?[1]:ID1_Gray > self[1]:ID1_Gray,ID1_Gray,ID1_Gray(circ x3.0 strong) |
| chr5R | ID1,ID3,ID1,ID3,ID1,ID3,ID1,ID3,ID1,ID3,ID1,ID3,ID1 | chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R | chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray |
| chr7L | ID1,ID3,ID3 | chr12R > chr14L | chr12R?[1]:ID1_Gray > chr14L[3-4]:ID3_Orange,ID3_Orange |
| chr13L | ID1,ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1_Gray > chr14L[3]:ID3_Orange |
| chr13R | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14L[3]:ID3_Orange |
| chr3R | ID5,ID3 | chr16L|chr7R > chr14L | chr16L|chr7R:ID5_Blue-Dark > chr14L[3]:ID3_Orange |
| chr13R | ID2,ID2,ID2 | chr15R > chr12R|chr4R | chr15R[1]:ID2_Red-Dark > chr12R|chr4R:ID2_Red-Light,ID2_Red-Light |
| chr4L | ID2,ID6 | chr12R|chr14L|chr4R > chr10L|chr9L | chr12R|chr14L|chr4R:ID2_Red-Light > chr10L|chr9L:ID6_Purple-Dark |
| chr11L | ID6,ID2,ID2,ID2,ID2,ID2,ID2 | chr10L|chr9L > chr15R > chr12R | chr10L|chr9L:ID6_Purple-Dark > chr15R[1]:ID2_Red-Dark > chr12R[2-6]:ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light |
| chr1R | ID1,ID4 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr2L|chr6L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr2L|chr6L:ID4_Green-Light |
| chr2R | ID1,ID5 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr7R | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr16L|chr7R:ID5_Blue-Dark |
| chr4L | ID3,ID2 | chr14L > chr15R | chr14L[3]:ID3_Orange > chr15R[1]:ID2_Red-Dark |
| chr4L | ID5,ID2 | chr16L|chr7R > chr12R|chr14L|chr4R | chr16L|chr7R:ID5_Blue-Dark > chr12R|chr14L|chr4R:ID2_Red-Light |
| chr1L | ID6,ID3 | chr14R > chr14L | chr14R[1]:ID6_Purple-Neutral > chr14L[3]:ID3_Orange |
| chr14R | ID6,ID2,ID6 | chr4R > chr14R | chr4R?[1]:ID2_Red-Light > self[1]:ID6_Purple-Neutral |
| chr11R | ID6,ID2,ID3,ID3 | chr10L|chr9L > chr14L | chr10L|chr9L:ID6_Purple-Dark > chr14L[2-4]:ID2_Red-Light,ID3_Orange,ID3_Orange |
| chr10R | ID1,ID1,ID1,ID1,ID1,ID1,ID3 | chr12R > chr14L | chr12R?[1]:ID1_Gray,ID1_Gray,ID1_Gray,ID1_Gray,ID1_Gray,ID1_Gray(circ x6.0 moderate) > chr14L[3]:ID3_Orange |
| chr13R | ID5,ID2 | chr16L|chr7R > chr15R | chr16L|chr7R:ID5_Blue-Dark > chr15R[1]:ID2_Red-Dark |
| chr10R | ID6,ID2 | chr10L|chr9L > chr12R|chr14L|chr4R | chr10L|chr9L:ID6_Purple-Dark > chr12R|chr14L|chr4R:ID2_Red-Light |
| chr7L | ID1,ID2,ID3,ID3 | chr12R > chr14L | chr12R[1-2]:ID1_Gray,ID2_Red-Light > chr14L[3-4]:ID3_Orange,ID3_Orange |
| chr5R | ID1,ID1,ID6 | chr5R > chr5L | self[1]:ID1_Gray > chr5L[1]:ID6_Purple-Light |
| chr15L | ID6,ID1 | chr14R > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R | chr14R[1]:ID6_Purple-Neutral > chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray |
| chr7L | ID1,ID6 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14R | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14R[1]:ID6_Purple-Neutral |
| chr7L | ID5,ID3,ID6,ID6,ID3,ID6,ID6 | chr16L|chr7R > chr14L > chr10L|chr9L > chr14L > chr10L|chr9L | chr16L|chr7R:ID5_Blue-Dark > chr14L[3]:ID3_Orange > chr10L|chr9L:ID6_Purple-Dark > chr10L|chr9L:ID6_Purple-Dark > chr14L[3]:ID3_Orange > chr10L|chr9L:ID6_Purple-Dark > chr10L|chr9L:ID6_Purple-Dark |
| chr9R | ID4,ID2 | chr2L|chr6L > chr15R | chr2L|chr6L:ID4_Green-Light > chr15R[1]:ID2_Red-Dark |
| chr7L | ID2,ID2,ID2 | chr4R > chr15R | chr4R[1-2]:ID2_Red-Light,ID2_Red-Light > chr15R[1]:ID2_Red-Dark |
| chr3L | ID5,ID2,ID2,ID2 | chr16L|chr7R > chr4R | chr16L|chr7R:ID5_Blue-Dark > chr4R[1-3]:ID2_Red-Light,ID2_Red-Light,ID2_Red-Light |
| chr4L | ID1,ID3,ID3,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14L[3-5]:ID3_Orange,ID3_Orange,ID3_Orange |
| chr5R | ID1,ID3,ID1,ID3,ID1,ID3,ID5,ID1,ID3,ID1,ID3,ID1,ID3,ID1 | chr14L > chr5R > chr14L > chr5R > chr14L > chr16L|chr7R > chr5R > chr14L > chr5R > chr14L > chr5R > chr14L > chr5R | chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > chr16L|chr7R:ID5_Blue-Dark > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray > chr14L[3]:ID3_Orange > self[1]:ID1_Gray |
| chr5R | ID1,ID2,ID6,ID2 | chr15R > chr10L|chr9L > chr15R | chr15R[1]:ID2_Red-Dark > chr10L|chr9L:ID6_Purple-Dark > chr15R[1]:ID2_Red-Dark |
| chr13R | ID5,ID2,ID3,ID3,ID3 | chr16L|chr7R > chr14L | chr16L|chr7R:ID5_Blue-Dark > chr12R|chr14L|chr4R:ID2_Red-Light > chr14L[3-5]:ID3_Orange,ID3_Orange,ID3_Orange |
| chr9L | ID6,ID2,ID2,ID2 | chr15R > chr12R|chr4R | chr15R[1]:ID2_Red-Dark > chr12R|chr4R:ID2_Red-Light,ID2_Red-Light |
| chr3L | ID1,ID3 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr14L[3]:ID3_Orange |
| chr3R | ID1,ID2 | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr15R | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R:ID1_Gray > chr15R[1]:ID2_Red-Dark |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID2,ID2,ID3,ID3,ID3,ID3,ID3 | chr15R > chr14L | chr15R[1]:ID2_Red-Dark > self[2-5]:ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange > self[3-4]:ID3_Orange,ID3_Orange |
| chr8R | ID1,ID3,ID1 | chr14L > chr8R | chr14L[3]:ID3_Orange > self[1]:ID1_Gray |
| chr10R | ID5,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L | chr16L|chr7R:ID5_Blue-Dark > chr14L[3]:ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange(circ x9.0 strong) |
| chr15R | ID2,ID2,ID2,ID2,ID2,ID2 | chr4R > chr15R > chr12R|chr14L|chr4R > chr15R | chr4R?[1]:ID2_Red-Light > self[1]:ID2_Red-Dark > chr12R|chr14L|chr4R:ID2_Red-Light > self[1]:ID2_Red-Dark > self[1]:ID2_Red-Dark |

_(19 more rows in the TSV)_

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr14L | ID5,ID2,ID3,ID3,ID2,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr12L | ID1,ID2,ID2,ID2,ID3 | chr12R|chr4R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID2,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr12L | ID3,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr14L > chr12R |
| chr16R | ID1,ID2,ID2,ID1 | chr4R > chr16R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID2,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr11L | ID1,ID2,ID2 | chr12L|chr13L|chr16R|chr5R|chr8L|chr8R > chr4R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr5R | ID1,ID2,ID1 | chr12R|chr4R > chr5R |
| chr13R | ID1,ID2,ID2,ID2 | chr12L|chr13L|chr16R|chr5R|chr8L|chr8R > chr12R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr14R | ID6,ID6,ID2 | chr14R > chr12R|chr14L|chr4R |
| chr5R | ID1,ID1,ID1,ID1 | chr5R > chr12R > chr5R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |

