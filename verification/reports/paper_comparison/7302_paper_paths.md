# Summary

| PD | reads | paper_Y | caught | missed | extra | agreement_pct |
|---|---|---|---|---|---|---|
| 28.0 | 103.0 | 9.0 | 8.0 | 1.0 | 4.0 | 95.1 |
| 33.0 | 86.0 | 7.0 | 6.0 | 1.0 | 5.0 | 93.0 |


# Our path parser applied to the paper's own Y' calls — strain 7302

reference: `7302_day0_with_selection`, curated 7302 library; orientation inferred from the data.

**Orientation heuristic: 187/187 reads (100.0 %) match the true telo_side.**

## PD 28 — 103 reads
* paper flags 9 switching reads; we catch 8, miss 1, and flag 4 it calls N (95.1 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr13L | ID7,ID2,ID8,ID2,ID8,ID2,ID8 | chr15R > chr13L > chr15R > chr13L > chr12R|chr4R > chr13L | chr15R[1]:ID2_Red-Dark > self[4]:ID8_Brown > chr15R[1]:ID2_Red-Dark > self[4]:ID8_Brown > chr12R|chr4R:ID2_Red-Light > self[4]:ID8_Brown |
| chr13R | ID4,ID2 | chr2L > chr12R|chr13L|chr14L|chr4R | chr2L[1]:ID4_Green-Light > chr12R|chr13L|chr14L|chr4R:ID2_Red-Light |
| chr14R | ID5,ID3,ID3,ID3 | chr16L|chr7R > chr14L | chr16L|chr7R:ID5_Blue-Dark > chr14L[3-5]:ID3_Orange,ID3_Orange,ID3_Orange |
| chr13R | ID5,ID7 | chr16L|chr7R > chr13L | chr16L|chr7R:ID5_Blue-Dark > chr13L[1]:ID7_Yellow |
| chr7R | ID5,ID5,ID3,ID5,ID5,ID3 | chr7R > chr14L > chr7R > chr14L | self[1]:ID5_Blue-Dark > chr14L[3]:ID3_Orange > self[1]:ID5_Blue-Dark > self[1]:ID5_Blue-Dark > chr14L[3]:ID3_Orange |
| chr13L | ID7,ID8,ID7,ID2,ID8,ID7,ID2,ID8,ID7,ID2 | chr13L | self[4]:ID8_Brown > self[1-2]:ID7_Yellow,ID2_Red-Light > self[4]:ID8_Brown > self[1-2]:ID7_Yellow,ID2_Red-Light > self[4]:ID8_Brown > self[1-2]:ID7_Yellow,ID2_Red-Light |
| chr15L | ID6,ID3,ID3 | chr10L|chr9L > chr14L | chr10L|chr9L:ID6_Purple-Dark > chr14L[3-4]:ID3_Orange,ID3_Orange |
| chr16R | ID1,ID3,ID3,ID2 | chr14L > chr15R | chr14L[3-4]:ID3_Orange,ID3_Orange > chr15R[1]:ID2_Red-Dark |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID2 | chr14L > chr15R | self[3]:ID3_Orange > chr15R[1]:ID2_Red-Dark |

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr13L | ID7,ID8,ID7,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | chr13L > chr4R |
| chr12L | ID3,ID2,ID2,ID2 | chr14L > chr12R|chr4R |
| chr16R | ID1,ID2,ID3 | chr13L > chr14L |
| chr13R | ID5,ID2,ID3,ID3,ID3,ID3 | ? > chr14L |


## PD 33 — 86 reads
* paper flags 7 switching reads; we catch 6, miss 1, and flag 5 it calls N (93.0 % agreement)

### the paper's switching reads as our parser reads them
| chr_end | ids | blocks | path |
|---|---|---|---|
| chr3R | ID6,ID2,ID3,ID3 | chr14R > chr14L | chr14R[1]:ID6_Purple-Neutral > chr14L[2-4]:ID2_Red-Light,ID3_Orange,ID3_Orange |
| chr13L | ID7,ID8,ID7,ID8,ID8,ID7,ID2,ID2,ID2,ID2,ID2 | chr13L | self[3-4]:ID8_Brown,ID7_Yellow,ID8_Brown(circ x1.5 moderate | alt chr13L[4] + chr13L[3-4]) > self[4]:ID8_Brown > self[3]:ID7_Yellow > self[2]:ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light(circ x5.0 strong) |
| chr16R | ID3,ID7,ID8 | chr14L > chr13L | chr14L[3]:ID3_Orange > chr13L[3-4]:ID7_Yellow,ID8_Brown |
| chr15R | ID2,ID2,ID8,ID8 | chr15R > chr13L | self[1]:ID2_Red-Dark > chr13L[4]:ID8_Brown > chr13L[4]:ID8_Brown |
| chr16R | ID3,ID7,ID8 | chr14L > chr13L | chr14L[3]:ID3_Orange > chr13L[1]:ID7_Yellow > chr13L[4]:ID8_Brown |
| chr1R | ID5,ID7,ID3 | chr16L|chr7R > chr13L > chr14L | chr16L|chr7R:ID5_Blue-Dark > chr13L[1]:ID7_Yellow > chr14L[3]:ID3_Orange |
| chr3R | ID5,ID8,ID7,ID8 | chr16L|chr7R > chr13L | chr16L|chr7R:ID5_Blue-Dark > chr13L[3-4]:ID8_Brown,ID7_Yellow,ID8_Brown(circ x1.5 moderate | alt chr13L[4] + chr13L[3-4]) |

### reads we call a switch and the paper does not
| chr_end | ids | blocks |
|---|---|---|
| chr2L | ID3,ID2,ID2 | chr14L > chr12R|chr4R |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | chr16L|chr7R > chr14L |
| chr12R | ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID3,ID3 | chr13L > chr14L |
| chr2R | ID3,ID2,ID2,ID2,ID2 | chr14L > chr12R |
| chr15R | ID2,ID2,ID2 | chr15R > chr12R|chr13L|chr14L|chr4R |

