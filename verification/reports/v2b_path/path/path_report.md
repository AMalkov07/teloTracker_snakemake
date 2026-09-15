# Y' path report -- verification/snapshot_v2b_path

## 7302_day2_with_selection
138 gain-like reads, 53 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 22 |
| self (tandem amplification) | 13 |
| unique donor by IDs+ITS | 11 |
| unresolved | 7 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 22 |
| chr4R | 6 |
| chr14L | 3 |
| chr12R | 2 |

### Circles from the path (8 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr13L[1-2] | ID2,ID1 | 2 | 1 | 2.5 | 2 | 0 | chr2R |
| chr12R|chr14L | ID1,ID1 | 1 | 1 | 2.0 | 0 | 0 | chr6R |
| chr13L[1-2] | ID1,ID2 | 1 | 1 | 2.0 | 0 | 0 | chr5L |
| chr14L[2-3] | ID1,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr8L |
| chr4R[5-7] | ID1,ID1,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr3R |
| self[2] | ID1 | 1 | 1 | 7.0 | 1 | 0 | chr12R |
| self[3-5] | ID2,ID2,ID2 | 1 | 1 | 1.7 | 0 | 0 | chr14L |

### Composite paths (>= 2 segments): 7 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr13L | self[1-3]:ID2,ID1,ID2 > chr2L|chr6L:ID4 |
| chr16R | chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr2L | chr13L?[1]:ID2 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr2R | chr13L?[1]:ID2 > chr13L[1-2]:ID2,ID1 |
| chr6L | chr13L?[1]:ID2 > chr13L[1-2]:ID2,ID1 |
| chr6L | chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr7L | chr13L[1-3]:ID2,ID1,ID2 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |


## 7302_day3_with_selection
309 gain-like reads, 129 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 60 |
| unique donor by IDs+ITS | 30 |
| self (tandem amplification) | 21 |
| unresolved | 18 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 40 |
| chr14L | 23 |
| chr4R | 19 |
| chr12R | 8 |

### Circles from the path (18 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr13L[1-2] | ID2,ID1 | 2 | 2 | 3.5 | 2 | 0 | chr5R,chr10R |
| chr12R|chr14L|chr4R | ID1,ID1 | 2 | 2 | 2.0 | 0 | 0 | chr6L,chr12L |
| chr12R[2-6] | ID1,ID1,ID1,ID1,ID1 | 1 | 1 | 1.4 | 1 | 0 | chr7L |
| chr12R[2-5] | ID1,ID1,ID1,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr14L |
| chr13L[2-3] | ID2,ID1 | 1 | 1 | 2.0 | 0 | 0 | chr3R |
| chr13L[2-4] | ID1,ID1,ID2 | 1 | 1 | 1.7 | 1 | 0 | chr2L |
| chr13L[2-4] | ID1,ID2,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr1L |
| chr14L[1-3] | ID1,ID2,ID1 | 1 | 1 | 1.7 | 0 | 0 | chr7R |
| chr14L[2-4] | ID2,ID2,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr16R |
| chr14L[3-4] | ID2,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr12L |
| chr4R[1-2] | ID1,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr7R |
| chr4R[3-4] | ID1,ID1 | 1 | 1 | 3.0 | 0 | 0 | chr13L |
| self[1-2] | ID1,ID2 | 1 | 1 | 2.0 | 0 | 0 | chr13L |
| self[2-4] | ID1,ID2,ID1 | 1 | 1 | 2.0 | 1 | 0 | chr13L |
| self[2] | ID1 | 1 | 1 | 6.0 | 0 | 0 | chr12R |
| self[3-4] | ID2,ID2 | 1 | 1 | 4.0 | 1 | 0 | chr14L |

### Composite paths (>= 2 segments): 17 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr11R | chr13L[1-3]:ID2,ID1,ID2 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr12R | self[2]:ID1 > chr10L|chr9L:ID3 |
| chr13L | self[1]:ID2 > self[2]:ID1 |
| chr13L | self[1]:ID2 > self[1]:ID2 |
| chr13R | chr14L[1-5]:ID1,ID1,ID2,ID2,ID2 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr13R | chr2L|chr6L:ID4 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr14R | chr10L|chr9L:ID3 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr15L | chr10L|chr9L:ID3 > chr14L[3-4]:ID2,ID2 |
| chr15L | chr14L?[3]:ID2 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr16R | chr14L[3-4]:ID2,ID2 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr1L | chr12L|chr14L|chr8L|chr8R:ID2 > chr13L[2-4]:ID1,ID2,ID1,ID1(circ x1.3 moderate | alt chr13L[2-4] + chr13L[2]) |
| chr2L | chr14L?[2]:ID1 > chr13L[2-4]:ID1,ID1,ID2,ID1,ID1(circ x1.7 strong) |


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


## 7302_day6_with_selection_repeat
301 gain-like reads, 69 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 50 |
| unresolved | 9 |
| unique donor by IDs+ITS | 9 |
| self (tandem amplification) | 1 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr14L | 27 |
| chr13L | 21 |
| chr4R | 6 |
| chr12R | 2 |
| chr5L | 1 |
| chr14R | 1 |
| chr16R | 1 |

### Circles from the path (6 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr12R|chr13L|chr14L|chr4R | ID1 | 1 | 1 | 7.0 | 0 | 0 | chr10R |
| chr12R|chr4R | ID1,ID1,ID1 | 1 | 1 | 1.7 | 0 | 0 | chr5L |
| chr13L|chr14L | ID2 | 1 | 1 | 3.0 | 0 | 0 | chr12L |
| chr14L[2-3] | ID1,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr5R |
| chr14L[3-4] | ID2,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr13R |
| chr2L|chr6L | ID4 | 1 | 1 | 7.0 | 0 | 0 | chr2R |

### Composite paths (>= 2 segments): 11 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10R | chr2L|chr6L:ID4 > chr12R[2-5]:ID1,ID1,ID1,ID1 > chr12R|chr13L|chr14L|chr4R:ID1,ID1,ID1,ID1,ID1,ID1,ID1(circ x7.0) |
| chr13R | chr2L|chr6L:ID4 > chr12R|chr14L|chr4R:ID1,ID1 > chr14L[3-4]:ID2,ID2,ID2(circ x1.5 moderate | alt chr14L[3] + chr14L[3-4]) > chr12R[2-5]:ID1,ID1,ID1,ID1 |
| chr14R | chr5L[1]:ID8 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr15L | chr2L|chr6L:ID4 > chr2L|chr6L:ID4 |
| chr15R | chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > self[1]:ID1 |
| chr2R | chr12L|chr14L|chr8L|chr8R:ID2 > chr12R[2-5]:ID1,ID1,ID1,ID1 |
| chr2R | chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr4R[1-3]:ID1,ID1,ID1 |
| chr2R | chr2L|chr6L:ID4 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr5R | chr14L?[3]:ID2 > self[1]:ID1 |
| chr7L | chr14R[1]:ID6 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr8L | chr16R[1]:ID7 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |


## 7302_day9_with_selection_repeat
257 gain-like reads, 64 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 52 |
| unique donor by IDs+ITS | 6 |
| unresolved | 5 |
| self (tandem amplification) | 1 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr14L | 40 |
| chr13L | 16 |
| chr12R | 1 |
| chr4R | 1 |

### Circles from the path (10 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr14L[1-3] | ID2,ID1,ID1 | 2 | 2 | 1.7 | 0 | 0 | chr1L,chr2R |
| chr14L[3-4] | ID2,ID2 | 2 | 2 | 2.0 | 0 | 0 | chr2R,chr5L |
| chr12R[2-5] | ID1,ID1,ID1,ID1 | 1 | 1 | 1.2 | 0 | 0 | chr14R |
| chr12L|chr13L|chr14L|chr8L|chr8R | ID2 | 1 | 1 | 6.0 | 0 | 0 | chr7L |
| chr14L[1-3] | ID1,ID2,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr5L |
| chr14L[2-3] | ID1,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr15R |
| chr14L[2-3] | ID2,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr2R |
| self[3-4] | ID2,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr14L |

### Composite paths (>= 2 segments): 14 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10L | chr14L?[3]:ID2 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr14R | chr12R[2-5]:ID1,ID1,ID1,ID1,ID1(circ x1.2 moderate | alt chr12R[2-5] + chr12R[2]) > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 > chr14L[3-4]:ID2,ID2 |
| chr14R | chr10L|chr9L:ID3 > chr14L[3-4]:ID2,ID2 |
| chr15L | chr16R[1]:ID7 > chr14L[3-4]:ID2,ID2 |
| chr1L | chr2L|chr6L:ID4 > chr10L|chr9L:ID3 |
| chr1R | chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 > chr14L?[3]:ID2 > chr12R|chr14L|chr4R:ID1,ID1 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr2R | chr14L[2-3]:ID2,ID1,ID2(circ x1.5 moderate | alt chr14L[3] + chr14L[2-3]) > chr14L[3-4]:ID2,ID2,ID2(circ x1.5 moderate | alt chr14L[3-4] + chr14L[3]) |
| chr2R | chr13L[1-2]:ID2,ID1 > chr12R|chr14L|chr4R:ID1,ID1 |
| chr5L | chr12L|chr14L|chr8L|chr8R:ID2 > chr14L[1-3]:ID1,ID2,ID1,ID1(circ x1.3 moderate | alt chr14L[2-3] + chr14L[1-2]) |
| chr5R | chr14L[2-3]:ID1,ID2 > chr14L[3-4]:ID2,ID2 |
| chr5R | chr10L|chr9L:ID3 > self[1]:ID1 |
| chr6L | chr14L?[2]:ID1 > chr14L[1-4]:ID1,ID1,ID2,ID2 |


## 7302_survivor_IT169
182 gain-like reads, 118 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unresolved | 45 |
| unique donor by IDs | 41 |
| unique donor by IDs+ITS | 32 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr13L | 27 |
| chr14L | 18 |
| chr12R | 14 |
| chr4R | 14 |

### Circles from the path (34 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr12R|chr13L|chr14L|chr4R | ID1 | 10 | 8 | 7.0 | 0 | 0 | chr1L,chr2R,chr5L,chr7L,chr7R,chr8L,chr13R,chr15L |
| chr12R|chr4R | ID1,ID1,ID1 | 3 | 2 | 3.0 | 0 | 0 | chr3L,chr14R |
| chr12R[2-6] | ID1,ID1,ID1,ID1,ID1 | 3 | 2 | 1.8 | 3 | 0 | chr7L,chr11L |
| chr12R[2-3] | ID1,ID1 | 2 | 2 | 2.5 | 0 | 0 | chr5L,chr13L |
| chr12R|chr13L|chr14L | ID1 | 2 | 2 | 6.0 | 0 | 0 | chr2R,chr3L |
| chr13L[1-2] | ID1,ID2 | 2 | 2 | 1.5 | 0 | 0 | chr7L,chr8L |
| chr12R|chr13L | ID1 | 2 | 2 | 9.0 | 0 | 0 | chr8L,chr15L |
| chr12R[2-4] | ID1,ID1,ID1 | 1 | 1 | 1.7 | 0 | 0 | chr5R |
| chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R | ID1 | 1 | 1 | 3.0 | 0 | 0 | chr7L |
| chr12R|chr4R | ID1,ID1 | 1 | 1 | 2.5 | 0 | 0 | chr6R |
| chr12R|chr14L|chr4R | ID1,ID1 | 1 | 1 | 2.0 | 0 | 0 | chr7L |
| chr14L[1-2] | ID1,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr9L |
| chr14L[1-3] | ID1,ID1,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr14R |
| chr14L[1-3] | ID2,ID1,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr16R |
| chr14L[1-4] | ID1,ID1,ID2,ID2 | 1 | 1 | 1.5 | 1 | 0 | chr14R |
| chr14L[2-3] | ID1,ID2 | 1 | 1 | 1.5 | 0 | 0 | chr8L |
| chr4R[2-6] | ID1,ID1,ID1,ID1,ID1 | 1 | 1 | 1.4 | 0 | 0 | chr7R |

### Composite paths (>= 2 segments): 56 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10L | chr13L?[2]:ID1 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr11L | chr10L|chr9L:ID3 > chr12R[2-6]:ID1,ID1,ID1,ID1,ID1,ID1,ID1(circ x1.4 strong) |
| chr11L | chr10L|chr9L:ID3 > chr12R[2-5]:ID1,ID1,ID1,ID1 |
| chr11L | chr10L|chr9L:ID3 > chr12R[2-4]:ID1,ID1,ID1 |
| chr11L | chr10L|chr9L:ID3 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr11L | chr10L|chr9L:ID3 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr11L | chr10L|chr9L:ID3 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr11L | chr10L|chr9L:ID3 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr11L | chr10L|chr9L:ID3 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr11L | chr10L|chr9L:ID3 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr11R | chr12R[3-6]:ID1,ID1,ID1,ID1 > chr14L[2-3]:ID1,ID2 |
| chr12L | chr12R[2-5]:ID1,ID1,ID1,ID1 > chr14L[3-4]:ID2,ID2 |


## 7302_survivor_IT171
125 gain-like reads, 71 with >= 2 gained Y'.
### Donor resolution of multi-Y' gains
| category | n_reads |
|---|---|
| unique donor by IDs | 33 |
| unresolved | 19 |
| unique donor by IDs+ITS | 19 |

### Donors (unique calls)
| donor | n_reads |
|---|---|
| chr14L | 29 |
| chr16R | 11 |
| chr13L | 8 |
| chr12R | 4 |

### Circles from the path (13 reads): donor x repeat unit
| donor | unit | n_reads | n_ends | max_repeats | n_strong | n_weak | ends |
|---|---|---|---|---|---|---|---|
| chr12R|chr13L | ID1 | 3 | 2 | 4.0 | 0 | 0 | chr3R,chr11L |
| chr14L[2-3] | ID1,ID2 | 2 | 2 | 2.0 | 0 | 0 | chr1L,chr16L |
| chr14L[1-3] | ID1,ID1,ID2 | 2 | 1 | 1.3 | 0 | 0 | chr13R |
| chr12R[2-5] | ID1,ID1,ID1,ID1 | 1 | 1 | 1.2 | 0 | 0 | chr1L |
| chr12R|chr13L|chr14L | ID1 | 1 | 1 | 3.0 | 0 | 0 | chr9L |
| chr12R|chr4R | ID1,ID1,ID1 | 1 | 1 | 1.3 | 0 | 0 | chr9L |
| chr14L[2-3] | ID2,ID1 | 1 | 1 | 1.5 | 0 | 0 | chr8R |
| chr14L[2-4] | ID2,ID1,ID2 | 1 | 1 | 1.3 | 0 | 0 | chr7R |
| chr16R[1] | ID7 | 1 | 1 | 3.0 | 1 | 0 | chr6R |

### Composite paths (>= 2 segments): 48 reads; examples
| chr_end | y_prime_path |
|---|---|
| chr10L | chr13L?[2]:ID1 > chr12R|chr13L|chr14L|chr15R|chr16L|chr4R|chr5R|chr7R:ID1 |
| chr10R | chr2L|chr6L:ID4 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 |
| chr11L | chr2L|chr6L:ID4 > chr14L?[3]:ID2 > chr12R[2-5]:ID1,ID1,ID1,ID1 > chr2L|chr6L:ID4 |
| chr11L | chr13L[1-2]:ID2,ID1 > chr12R|chr13L:ID1,ID1,ID1,ID1(circ x4.0) |
| chr11L | chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr14L?[2]:ID1 > chr12R|chr13L:ID1,ID1,ID1,ID1(circ x4.0) |
| chr12L | self[1]:ID2 > chr13L[1-2]:ID2,ID1 |
| chr12L | chr13L?[1]:ID2 > self[1]:ID2 |
| chr13L | chr12R[2-4]:ID1,ID1,ID1 > chr14L?[3]:ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2(circ x8.0 strong) |
| chr13L | chr14L?[2]:ID1 > chr14L[1-4]:ID1,ID1,ID2,ID2 |
| chr13R | chr14L?[2]:ID1 > chr14L[1-3]:ID1,ID1,ID2 |
| chr14R | chr16R[1]:ID7 > chr14L?[3]:ID2,ID2,ID2,ID2,ID2,ID2,ID2(circ x7.0 strong) |
| chr14R | chr16R[1]:ID7 > chr16R[1]:ID7 > chr14L[2-3]:ID1,ID2 |


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

