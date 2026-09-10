# Y' path report -- verification/snapshot_curated_path

## 7302_day4_with_selection
271 gain-like reads, 124 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 66 |
| unresolved | 22 |
| unique donor by IDs+ITS | 19 |
| self (tandem amplification) | 17 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 35 |
| chr14L | 30 |
| chr4R | 14 |
| chr12R | 6 |

### Circles from the path (29 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr13L[1-2] | ID7,ID2 | 3 | 3 | 2.5 | 2 | 0 | chr2R,chr7R,chr10R |
| chr14L[3-4] | ID3,ID3 | 2 | 2 | 3.0 | 1 | 0 | chr8R,chr12L |
| chr13L[1-2] | ID2,ID7 | 2 | 2 | 2.5 | 1 | 0 | chr5R,chr8L |
| chr12R|chr4R | ID2,ID2 | 2 | 2 | 1.5 | 0 | 0 | chr7R,chr14L |
| chr12R|chr13L | ID2 | 2 | 2 | 3.0 | 0 | 0 | chr5L,chr15R |
| self[3] | ID3 | 2 | 1 | 7.0 | 1 | 0 | chr14L |
| self[2] | ID2 | 2 | 1 | 4.0 | 0 | 0 | chr13L |
| chr12R[2-4] | ID2,ID2,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr7L |
| chr12R[2-3] | ID2,ID2 | 1 | 1 | 2.0 | 0 | 0 | chr15R |
| chr12R|chr13L|chr14L|chr4R | ID2 | 1 | 1 | 6.0 | 0 | 0 | chr12L |
| chr12R[2-6] | ID2,ID2,ID2,ID2,ID2 | 1 | 1 | 1.4 | 1 | 0 | chr15R |
| chr12R|chr4R | ID2,ID2,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr9L |
| chr13L[2-3] | ID7,ID2 | 1 | 1 | 2.5 | 1 | 0 | chr4L |
| chr14L[3-5] | ID3,ID3,ID3 | 1 | 1 | 1.3 | 0 | 0 | chr14R |
| chr14L[1-4] | ID3,ID3,ID5,ID2 | 1 | 1 | 1.2 | 0 | 0 | chr8L |
| chr14L[3] | ID3 | 1 | 1 | 4.0 | 0 | 0 | chr14R |
| chr4R[1-4] | ID2,ID2,ID2,ID2 | 1 | 1 | 1.2 | 0 | 0 | chr7R |
| chr4R[3-5] | ID2,ID2,ID2 | 1 | 1 | 1.7 | 0 | 0 | chr4L |
| chr4R[3-4] | ID2,ID2 | 1 | 1 | 3.5 | 0 | 0 | chr12L |
| self[1-3] | ID7,ID7,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr13L |
| self[1-2] | ID5,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr14L |

### Composite paths (>= 2 segments): 32 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10R | chr14L|chr16L|chr7R:ID5 > chr14L[2-5]:ID2,ID3,ID3,ID3 > chr14L[3]:ID3 |
| chr11L | chr13L[1]:ID7 > chr12R|chr13L|chr14L|chr15R|chr4R:ID2 |
| chr12L | chr14L[3]:ID3 > chr12R|chr13L|chr14L|chr4R:ID2,ID2,ID2,ID2,ID2,ID2(circ x6.0) |
| chr12L | chr13L?[2]:ID2 > chr12R|chr13L|chr14L|chr15R|chr4R:ID2 |
| chr12R | chr13L?[2]:ID2 > chr14L[3-4]:ID3,ID3 |
| chr13L | self[2-4]:ID2,ID7,ID2 > self[2]:ID2,ID2,ID2,ID2(circ x4.0 moderate) |
| chr13L | chr14L?[2]:ID2 > self[2]:ID2,ID2,ID2(circ x3.0 moderate) |
| chr14L | self[1]:ID5 > self[1-2]:ID5,ID2,ID5(circ x1.5 moderate | alt chr14L[1-2] + chr14L[1]) > self[1-2]:ID5,ID2 |
| chr14L | self[3]:ID3 > self[2]:ID2 |
| chr14R | self[1]:ID6 > self[1]:ID6 |
| chr15R | chr14L[1-2]:ID5,ID2 > self[1]:ID2 |
| chr16R | chr14L[3]:ID3 > chr12R|chr4R:ID2,ID2 |


## 7302_day5_with_selection
392 gain-like reads, 203 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 117 |
| unresolved | 37 |
| unique donor by IDs+ITS | 26 |
| self (tandem amplification) | 23 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 56 |
| chr14L | 51 |
| chr4R | 27 |
| chr12R | 9 |

### Circles from the path (47 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr13L[1-2] | ID7,ID2 | 4 | 4 | 2.5 | 1 | 1 | chr1R,chr2R,chr5R,chr15R |
| chr14L[3] | ID3 | 4 | 4 | 6.0 | 1 | 1 | chr2R,chr8L,chr11L,chr16R |
| chr13L[1-3] | ID7,ID2,ID7 | 3 | 3 | 2.0 | 3 | 0 | chr2R,chr3R,chr11R |
| chr12R|chr4R | ID2,ID2 | 3 | 3 | 3.5 | 0 | 0 | chr5L,chr5R,chr15R |
| chr13L[2-3] | ID2,ID7 | 2 | 2 | 1.5 | 0 | 0 | chr2L,chr12L |
| chr13L[2-3] | ID7,ID2 | 2 | 2 | 2.0 | 0 | 0 | chr2R,chr4L |
| chr13L[1-2] | ID2,ID7 | 2 | 2 | 2.0 | 0 | 0 | chr3R,chr14R |
| chr4R[2-5] | ID2,ID2,ID2,ID2 | 2 | 2 | 1.2 | 0 | 0 | chr5L,chr14R |
| self[2-4] | ID2,ID2,ID2 | 2 | 1 | 1.3 | 0 | 0 | chr12R |
| chr12R[2-5] | ID2,ID2,ID2,ID2 | 1 | 1 | 1.8 | 1 | 0 | chr8R |
| chr12R|chr13L|chr14L|chr4R | ID2 | 1 | 1 | 9.0 | 0 | 0 | chr12L |
| chr12R|chr13L|chr14L|chr15R|chr4R | ID2 | 1 | 1 | 3.0 | 0 | 0 | chr3R |
| chr13L[1-4] | ID7,ID2,ID7,ID2 | 1 | 1 | 1.2 | 1 | 0 | chr11L |
| chr14L[1-2] | ID2,ID5 | 1 | 1 | 1.5 | 0 | 0 | chr16L |
| chr12R|chr13L|chr14L | ID2 | 1 | 1 | 3.0 | 0 | 0 | chr5R |
| chr12R[2-3] | ID2,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr16R |
| chr14L[2-3] | ID3,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr5L |
| chr14L[2-3] | ID2,ID3 | 1 | 1 | 1.5 | 0 | 0 | chr7R |
| chr4R[1-5] | ID2,ID2,ID2,ID2,ID2 | 1 | 1 | 1.2 | 0 | 0 | chr16R |
| chr4R[1-4] | ID2,ID2,ID2,ID2 | 1 | 1 | 1.2 | 0 | 0 | chr7R |
| chr4R[3-4] | ID2,ID2 | 1 | 1 | 2.5 | 0 | 0 | chr8L |
| chr4R[3-5] | ID2,ID2,ID2 | 1 | 1 | 1.7 | 0 | 0 | chr13R |
| chr4R[3-7] | ID2,ID2,ID2,ID2,ID2 | 1 | 1 | 1.2 | 0 | 0 | chr2L |
| chr4R[2-4] | ID2,ID2,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr13R |
| chr4R[4-7] | ID2,ID2,ID2,ID2 | 1 | 1 | 1.2 | 0 | 0 | chr1L |
| self[1-2] | ID2,ID7 | 1 | 1 | 3.5 | 1 | 0 | chr13L |
| self[1-3] | ID7,ID7,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr13L |
| self[1-5] | ID2,ID2,ID2,ID2,ID2 | 1 | 1 | 1.2 | 0 | 0 | chr4R |
| self[2-4] | ID2,ID2,ID7 | 1 | 1 | 1.3 | 0 | 0 | chr13L |
| self[2-4] | ID2,ID3,ID3 | 1 | 1 | 1.3 | 0 | 0 | chr14L |

_(2 more rows in the TSV)_

### Composite paths (>= 2 segments): 71 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10L | chr14L[3]:ID3 > chr14L[3]:ID3 |
| chr10L | chr14L?[1]:ID5 > chr14L[3]:ID3 |
| chr10R | chr13L[1]:ID7 > chr13L[1]:ID7 |
| chr11L | chr13L[1-4]:ID7,ID2,ID7,ID2,ID7(circ x1.2 strong) > chr14L[3]:ID3,ID3,ID3,ID3,ID3,ID3(circ x6.0 strong) |
| chr11L | chr14L[1-3]:ID5,ID2,ID3 > chr14L[2-3]:ID2,ID3 |
| chr11L | chr2L|chr6L:ID4 > chr12R|chr13L|chr14L|chr15R|chr4R:ID2 |
| chr11R | chr13L[1-3]:ID7,ID2,ID7,ID7,ID2(circ x1.7 strong) > chr12R|chr13L|chr14L|chr15R|chr4R:ID2 |
| chr12L | chr14L[3]:ID3 > chr12R|chr4R:ID2,ID2 |
| chr12R | chr13L?[2]:ID2 > self[2-4]:ID2,ID2,ID2,ID2(circ x1.3 moderate | alt chr12R[2] + chr12R[2-4]) |
| chr13L | self[1]:ID7 > self[1-2]:ID2,ID7,ID2,ID7,ID2,ID7,ID2(circ x3.5 strong) > self[2-3]:ID2,ID7 |
| chr13L | self[1-2]:ID7,ID2 > self[2-4]:ID2,ID2,ID7,ID2(circ x1.3 moderate | alt chr13L[2] + chr13L[2-4]) |
| chr13L | self[1]:ID7 > self[1]:ID7 |


## 7302_day0_with_selection
117 gain-like reads, 19 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 11 |
| unresolved | 4 |
| self (tandem amplification) | 2 |
| unique donor by IDs+ITS | 2 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 5 |
| chr14L | 4 |
| chr4R | 3 |
| chr12R | 1 |

### Circles from the path (5 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr12R[2-3] | ID2,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr12L |
| chr12R[2-6] | ID2,ID2,ID2,ID2,ID2 | 1 | 1 | 1.2 | 0 | 0 | chr6L |
| chr13L[1-2] | ID7,ID2 | 1 | 1 | 2.0 | 0 | 0 | chr13R |
| chr14L[3-4] | ID3,ID3 | 1 | 1 | 1.5 | 0 | 0 | chr8L |
| chr14L[3] | ID3 | 1 | 1 | 5.0 | 0 | 0 | chr10R |

### Composite paths (>= 2 segments): 6 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10R | chr14L?[1]:ID5 > chr14L[3]:ID3,ID3,ID3,ID3,ID3(circ x5.0 moderate | alt chr14L[3-5] + chr14L[3-4]) |
| chr12L | chr14L[3]:ID3 > chr12R[2-3]:ID2,ID2,ID2(circ x1.5 moderate | alt chr12R[2-3] + chr12R[2]) |
| chr13L | chr10L|chr14R|chr5L|chr9L:ID6 > self[1-2]:ID7,ID2 |
| chr13R | chr10L|chr14R|chr5L|chr9L:ID6 > chr12R|chr13L|chr14L|chr15R|chr4R:ID2 |
| chr13R | chr10L|chr14R|chr5L|chr9L:ID6 > chr12R|chr13L|chr14L|chr15R|chr4R:ID2 |
| chr1L | chr10L|chr14R|chr5L|chr9L:ID6 > chr13L[2-4]:ID2,ID7,ID2 |

