# Circle candidates -- snapshot verification/snapshot_v2c_path

## 7302_day4_with_selection
250 gain-like reads; **52 circle candidates** ({'homopolymer_run': 38, 'repeated_unit_from_donor': 14}); ends with >= 3 candidate reads: chr2R, chr3R, chr4L, chr5R, chr12L, chr14L, chr14R, chr15R

### Recurrent (end, class, donor, repeat unit) with >= 2 reads
| chr_end | circle_class | donor | unit | n_reads | max_copies | pct_of_end |
|---|---|---|---|---|---|---|
| chr15R | homopolymer_run |  | ID1 | 4 | 4.0 | 2.4 |
| chr2R | repeated_unit_from_donor | chr13L | ID2,ID1 | 4 | 2.5 | 1.6 |
| chr3R | repeated_unit_from_donor | chr13L | ID2,ID1 | 3 | 2.0 | 0.9 |
| chr14L | homopolymer_run |  | ID1 | 3 | 3.0 | 3.4 |
| chr14L | homopolymer_run |  | ID2 | 3 | 7.0 | 3.4 |
| chr7R | homopolymer_run |  | ID1 | 2 | 5.0 | 0.9 |
| chr12R | homopolymer_run |  | ID1 | 2 | 3.0 | 4.2 |
| chr13L | homopolymer_run |  | ID1 | 2 | 4.0 | 1.6 |
| chr10R | repeated_unit_from_donor | chr13L | ID2,ID1 | 2 | 2.5 | 1.2 |
| chr14R | homopolymer_run | chr14L | ID2 | 2 | 4.0 | 0.8 |
| chr4L | repeated_unit_from_donor | chr13L | ID2,ID1 | 2 | 2.5 | 1.7 |
| chr12L | homopolymer_run | chr4R | ID1 | 2 | 7.0 | 0.8 |

### Donor ends of circle candidates
| donor | n |
|---|---|
| chr13L | 14 |
| chr14L | 7 |
| chr4R | 4 |


## 7302_day5_with_selection
376 gain-like reads; **76 circle candidates** ({'homopolymer_run': 53, 'repeated_unit_from_donor': 18, 'tandem_amplification_own_array': 2, 'repeated_unit_donor_ambiguous': 3}); ends with >= 3 candidate reads: chr2R, chr3R, chr4R, chr5L, chr7R, chr11R, chr12L, chr12R, chr13R, chr14R, chr15R, chr16L

### Recurrent (end, class, donor, repeat unit) with >= 2 reads
| chr_end | circle_class | donor | unit | n_reads | max_copies | pct_of_end |
|---|---|---|---|---|---|---|
| chr2R | repeated_unit_from_donor | chr13L | ID2,ID1 | 5 | 2.0 | 2.1 |
| chr15R | homopolymer_run |  | ID1 | 5 | 5.0 | 2.5 |
| chr12L | homopolymer_run |  | ID1 | 4 | 3.0 | 1.4 |
| chr12R | homopolymer_run |  | ID1 | 4 | 5.0 | 5.3 |
| chr4R | homopolymer_run |  | ID1 | 3 | 10.0 | 4.3 |
| chr16L | homopolymer_run | chr14L | ID2 | 3 | 3.0 | 1.5 |
| chr7R | homopolymer_run |  | ID1 | 3 | 5.0 | 1.3 |
| chr14R | homopolymer_run |  | ID1 | 2 | 5.0 | 0.8 |
| chr5L | homopolymer_run |  | ID1 | 2 | 3.0 | 0.9 |
| chr13R | homopolymer_run | chr4R | ID1 | 2 | 10.0 | 1.0 |
| chr2L | homopolymer_run |  | ID1 | 2 | 4.0 | 0.9 |
| chr4L | repeated_unit_from_donor | chr13L | ID2,ID1 | 2 | 2.0 | 1.2 |
| chr3R | repeated_unit_from_donor | chr13L | ID2,ID1 | 2 | 2.0 | 0.7 |

### Donor ends of circle candidates
| donor | n |
|---|---|
| chr13L | 18 |
| chr4R | 11 |
| chr14L | 7 |

