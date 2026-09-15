# Y' path report -- verification/snapshot_v2path

## 7302_day4_with_selection
250 gain-like reads, 123 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 59 |
| unresolved | 25 |
| unique donor by IDs+ITS | 23 |
| self (tandem amplification) | 16 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 33 |
| chr14L | 28 |
| chr4R | 15 |
| chr12R | 6 |

### Circles from the path (29 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr13L[1-2] | ID2,ID1 | 4 | 4 | 2.5 | 3 | 0 | chr2R,chr4L,chr7R,chr10R |
| chr12R|chr13L | ID1 | 2 | 2 | 3.0 | 0 | 0 | chr5L,chr15R |
| chr13L[1-2] | ID1,ID2 | 2 | 2 | 2.5 | 1 | 0 | chr5R,chr8L |
| chr12R[2-4] | ID1,ID1,ID1 | 2 | 2 | 1.7 | 0 | 0 | chr7L,chr15R |
| self[3] | ID2 | 2 | 1 | 7.0 | 1 | 0 | chr14L |
| self[2] | ID1 | 2 | 1 | 4.0 | 0 | 0 | chr13L |
| chr14L[3-4] | ID2,ID2 | 2 | 2 | 2.0 | 0 | 0 | chr8R,chr14R |
| chr12R[2-3] | ID1,ID1 | 1 | 1 | 2.0 | 0 | 0 | chr15R |
| chr10L|chr9L | ID3 | 1 | 1 | 3.0 | 0 | 0 | chr14R |
| chr12R|chr14L|chr4R | ID1,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr7R |
| chr12R|chr4R | ID1,ID1,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr9L |
| chr12R|chr13L|chr14L|chr4R | ID1 | 1 | 1 | 6.0 | 0 | 0 | chr12L |
| chr14L[1-4] | ID2,ID2,ID1,ID1 | 1 | 1 | 1.2 | 0 | 0 | chr8L |
| chr13L[2-4] | ID1,ID1,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr3R |
| chr4R[1-4] | ID1,ID1,ID1,ID1 | 1 | 1 | 1.2 | 0 | 0 | chr7R |
| chr4R[3-4] | ID1,ID1 | 1 | 1 | 3.5 | 0 | 0 | chr12L |
| self[1-2] | ID1,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr14L |
| chr4R[3-5] | ID1,ID1,ID1 | 1 | 1 | 1.7 | 0 | 0 | chr4L |
| self[1] | ID1 | 1 | 1 | 6.0 | 1 | 0 | chr14L |
| self[1-3] | ID2,ID2,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr13L |

### Composite paths (>= 2 segments): 24 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10R | chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 > chr12R|chr14L|chr4R:ID1 > chr14L?[3]:ID2,ID2,ID2,ID2(circ x4.0 moderate | alt chr14L[3-5] + chr14L[3]) |
| chr11L | chr14L?[3]:ID2 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr12L | chr13L?[2]:ID1 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr12R | chr14L?[1]:ID1 > chr14L[3-4]:ID2,ID2 |
| chr13L | self[2-4]:ID1,ID2,ID1 > self[2]:ID1,ID1,ID1,ID1(circ x4.0 moderate) |
| chr13L | chr14L?[2]:ID1 > self[2]:ID1,ID1,ID1(circ x3.0 moderate) |
| chr15R | chr12R|chr14L|chr4R:ID1,ID1 > chr12R[2-4]:ID1,ID1,ID1,ID1,ID1(circ x1.7 moderate | alt chr12R[2-4] + chr12R[2-3]) |
| chr16R | chr14L?[3]:ID2 > chr12R|chr14L|chr4R:ID1,ID1 |
| chr16R | chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr13L[1-2]:ID2,ID1 |
| chr16R | chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr13L[1-2]:ID2,ID1 |
| chr1R | chr13L[2-3]:ID1,ID2 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr1R | chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 > chr14R[1]:ID6 |


## 7302_day5_with_selection
376 gain-like reads, 202 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 100 |
| unique donor by IDs+ITS | 40 |
| unresolved | 39 |
| self (tandem amplification) | 23 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 52 |
| chr14L | 47 |
| chr4R | 31 |
| chr12R | 9 |
| chr14R | 1 |

### Circles from the path (53 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr13L[1-3] | ID2,ID1,ID2 | 4 | 4 | 2.0 | 3 | 0 | chr2R,chr3R,chr7R,chr11R |
| chr12R|chr14L|chr4R | ID1,ID1 | 3 | 3 | 1.5 | 0 | 0 | chr1L,chr5L,chr15R |
| chr13L[1-2] | ID1,ID2 | 3 | 3 | 2.0 | 0 | 0 | chr1L,chr3R,chr14R |
| chr14L[2-3] | ID2,ID1 | 3 | 3 | 2.0 | 0 | 0 | chr2R,chr5L,chr5R |
| chr13L[1-2] | ID2,ID1 | 2 | 2 | 2.5 | 1 | 1 | chr1R,chr15R |
| chr12R|chr13L|chr14L|chr4R | ID1 | 2 | 2 | 9.0 | 0 | 0 | chr5R,chr12L |
| chr13L[2-3] | ID1,ID2 | 2 | 2 | 1.5 | 0 | 0 | chr2L,chr12L |
| chr12R[2-5] | ID1,ID1,ID1,ID1 | 2 | 2 | 1.8 | 1 | 0 | chr8R,chr15R |
| chr14L[2-3] | ID1,ID2 | 2 | 2 | 2.0 | 0 | 0 | chr7R,chr11L |
| self[1-5] | ID1,ID1,ID1,ID1,ID1 | 2 | 1 | 1.2 | 0 | 0 | chr4R |
| self[2-4] | ID1,ID1,ID1 | 2 | 1 | 1.3 | 0 | 0 | chr12R |
| chr4R[3-7] | ID1,ID1,ID1,ID1,ID1 | 2 | 2 | 1.4 | 0 | 0 | chr2L,chr11R |
| chr4R[2-7] | ID1,ID1,ID1,ID1,ID1,ID1 | 2 | 2 | 1.2 | 0 | 0 | chr1R,chr9L |
| chr4R[2-5] | ID1,ID1,ID1,ID1 | 2 | 2 | 1.2 | 0 | 0 | chr5L,chr14R |
| chr4R[1-5] | ID1,ID1,ID1,ID1,ID1 | 2 | 2 | 1.2 | 0 | 0 | chr13R,chr16R |
| chr4R[1-4] | ID1,ID1,ID1,ID1 | 2 | 2 | 1.2 | 0 | 0 | chr7R,chr13R |
| chr14L[3-4] | ID2,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr5R |
| chr14L[1-2] | ID1,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr3R |
| chr13L[2-4] | ID1,ID1,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr16R |
| chr13L[2-4] | ID1,ID2,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr3R |
| chr13L[2-3] | ID2,ID1 | 1 | 1 | 2.0 | 0 | 0 | chr2R |
| chr12R|chr4R | ID1,ID1,ID1 | 1 | 1 | 2.0 | 0 | 0 | chr5R |
| chr12L|chr13L|chr14L|chr8L|chr8R | ID2 | 1 | 1 | 3.0 | 0 | 0 | chr16R |
| chr4R[2-4] | ID1,ID1,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr13R |
| chr4R[3-4] | ID1,ID1 | 1 | 1 | 2.5 | 0 | 0 | chr8L |
| chr4R[2-6] | ID1,ID1,ID1,ID1,ID1 | 1 | 1 | 1.2 | 0 | 0 | chr1L |
| self[1-3] | ID2,ID2,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr13L |
| self[1-2] | ID1,ID2 | 1 | 1 | 3.5 | 1 | 0 | chr13L |
| self[2-4] | ID1,ID1,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr13L |
| self[2-4] | ID1,ID2,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr14L |

_(2 more rows in the TSV)_

### Composite paths (>= 2 segments): 61 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10L | chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr10R | chr13L?[1]:ID2 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr11L | chr13L[1-4]:ID2,ID1,ID2,ID1 > chr14L?[3]:ID2,ID2,ID2,ID2,ID2,ID2,ID2(circ x7.0 strong) |
| chr11L | chr14L?[2]:ID1 > chr14L[2-3]:ID1,ID2,ID1,ID2(circ x2.0 moderate | alt chr14L[2-3] + chr14L[2-3]) |
| chr11L | chr2L|chr6L:ID4 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr11R | chr13L[1-3]:ID2,ID1,ID2,ID2,ID1(circ x1.7 strong) > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr12R | chr13L?[2]:ID1 > self[2-4]:ID1,ID1,ID1,ID1(circ x1.3 moderate | alt chr12R[2] + chr12R[2-4]) |
| chr13L | chr14L?[3]:ID2 > self[1-2]:ID1,ID2,ID1,ID2,ID1,ID2,ID1(circ x3.5 strong) > self[2-3]:ID1,ID2 |
| chr13L | self[1-2]:ID2,ID1 > self[2-4]:ID1,ID1,ID2,ID1(circ x1.3 moderate | alt chr13L[2] + chr13L[2-4]) |
| chr13L | self[1]:ID2 > self[1]:ID2 |
| chr13L | self[1]:ID2 > self[1]:ID2 |
| chr13R | chr4R[2-5]:ID1,ID1,ID1,ID1 > chr4R[1-5]:ID1,ID1,ID1,ID1,ID1,ID1(circ x1.2 moderate | alt chr4R[3-7] + chr4R[1]) |


## 7302_day0_with_selection
66 gain-like reads, 19 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 10 |
| unique donor by IDs+ITS | 5 |
| self (tandem amplification) | 2 |
| unresolved | 2 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 5 |
| chr14L | 3 |
| chr4R | 3 |
| chr12R | 2 |
| chr14R | 2 |

### Circles from the path (4 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr12R[2-3] | ID1,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr12L |
| chr12R[2-6] | ID1,ID1,ID1,ID1,ID1 | 1 | 1 | 1.2 | 0 | 0 | chr6L |
| chr14L[2-3] | ID2,ID1 | 1 | 1 | 2.0 | 0 | 0 | chr13R |
| chr14L[3-4] | ID2,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr8L |

### Composite paths (>= 2 segments): 5 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10R | chr14L[2-5]:ID1,ID2,ID2,ID2 > chr14L[3-4]:ID2,ID2 |
| chr13L | chr10L|chr9L:ID3 > self[1-2]:ID2,ID1 |
| chr13R | chr14R[1]:ID6 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr13R | chr14R[1]:ID6 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr1L | chr14R[1]:ID6 > chr13L[2-4]:ID1,ID2,ID1 |

