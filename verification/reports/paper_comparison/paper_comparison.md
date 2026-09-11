# Summary

| strain | PD | paper_reads | found_in_ours | chr_end_agree | yprime_count_agree | paper_Y | caught | missed | extra_switches | agreement_pct |
|---|---|---|---|---|---|---|---|---|---|---|
| 7172.0 | 28.0 | 60.0 | 59.0 | 59.0 | 53.0 | 9.0 | 9.0 | 0.0 | 15.0 | 74.6 |
| 7172.0 | 33.0 | 125.0 | 119.0 | 119.0 | 114.0 | 24.0 | 22.0 | 0.0 | 17.0 | 85.7 |
| 7302.0 | 28.0 | 103.0 | 103.0 | 103.0 | 102.0 | 9.0 | 7.0 | 2.0 | 7.0 | 91.3 |
| 7302.0 | 33.0 | 86.0 | 83.0 | 83.0 | 80.0 | 7.0 | 6.0 | 1.0 | 11.0 | 85.5 |


# Supplementary Data 6 vs our curated-library calls

snapshot: `verification/snapshot_curated_all`

## strain 7172, PD 28  ->  7172_day3_with_selection (60)
* 59 of 60 paper reads found in our output; chr_end agrees 59/59; Y' copy count agrees 53/59
* template switching: paper flags 9; we catch 9, miss 0, and flag 15 reads the paper calls N (74.6 % agreement)

### the paper's switching reads, as our pipeline called them
| paper_chr_end | paper_ids | our_ids | our_blocks | our_source | our_mechanism |
|---|---|---|---|---|---|
| chr6L | ID5,ID3,ID3 | ID3_Orange,ID3_Orange,ID5_Blue-Dark | chr14L > chr16L|chr7R | ambiguous | unmatched_array |
| chr1R | ID3,ID1 | ID1_Gray,ID3_Orange | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr14L | chr1L | subtelomere_switch |
| chr15L | ID4,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID4_Green-Light | chr16L|chr7R > chr14L > chr2L|chr6L | chr14L | subtelomere_switch |
| chr11R | ID3,ID3,ID3,ID3,ID3,ID4 | ID4_Green-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr2L|chr6L > chr14L | chr14L | donor_transfer |
| chr11L | ID3,ID3,ID3,ID3,ID6 | ID6_Purple-Neutral,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr14R > chr14L | chr11R | subtelomere_switch |
| chr16R | ID6,ID2,ID2,ID1 | ID1_Gray,ID2_Red,ID2_Red,ID6_Purple-Dark | chr12R|chr4R > chr10L|chr9L | chr12R | donor_transfer_candidates:3 |
| chr9L | ID5,ID3,ID3,ID5,ID6 | ID6_Purple-Dark,ID5_Blue-Dark,ID3_Orange,ID3_Orange,ID5_Blue-Light | chr16L > chr14L | ambiguous | unmatched_array |
| chr1L | ID5,ID5,ID5 | ID5_Blue-Dark,ID5_Blue-Light,ID5_Blue-Light | chr16L > chr14L | ambiguous | unmatched_array |
| chr13R | ID6,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID1 | ID1_Gray,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID6_Purple-Dark | chr12R > chr10L|chr9L | chr12R | donor_transfer |

### reads we call a switch and the paper does not
| paper_chr_end | paper_ids | our_ids | our_blocks |
|---|---|---|---|
| chr16L | ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red | chr12R > chr16L |
| chr13R | ID2,ID2,ID2,ID5 | ID5_Blue-Light,ID2_Red,ID2_Red,ID2_Red | chr14L > chr12R|chr4R |
| chr14L | ID2,ID2,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID2_Red,ID2_Red | chr7R > chr14L |
| chr14L | ID2,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID3_Orange,ID3_Orange,ID2_Red | chr16L > chr14L |
| chr15L | ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange | chr16L|chr7R > chr14L |
| chr14L | ID3,ID2,ID3,ID3,ID2,ID3,ID3,ID2,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID2_Red,ID3_Orange,ID3_Orange,ID2_Red,ID3_Orange,ID3_Orange,ID2_Red,ID3_Orange | chr7R > chr14L |
| chr13R | ID3,ID2,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID3_Orange,ID2_Red,ID3_Orange | chr16L > chr14L |
| chr14L | ID2,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID3_Orange,ID2_Red | chr7R > chr14L |
| chr15R | ID2,ID3,ID2 | ID2_Red,ID3_Orange,ID2_Red | chr14L > chr15R |
| chr6L | ID2,ID2,ID2,ID2,ID2,ID2,ID4 | ID4_Green-Light,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red | chr14L > chr12R |
| chr2R | ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange | chr7R > chr14L |
| chr10R | ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange | chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID5,ID2 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID3_Orange,ID5_Blue-Light,ID2_Red | chr7R > chr14L |
| chr14L | ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red | chr16L > chr4R > chr12R |
| chr14L | ID2,ID3,ID2,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID2_Red,ID3_Orange,ID2_Red | chr7R > chr14L |


## strain 7172, PD 33  ->  7172_day4_with_selection_repeat (84), 7172_day4_with_selection (41)
* 119 of 125 paper reads found in our output; chr_end agrees 119/119; Y' copy count agrees 114/119
* template switching: paper flags 24; we catch 22, miss 0, and flag 17 reads the paper calls N (85.7 % agreement)

### the paper's switching reads, as our pipeline called them
| paper_chr_end | paper_ids | our_ids | our_blocks | our_source | our_mechanism |
|---|---|---|---|---|---|
| chr14L | ID4,ID2,ID4,ID2,ID4,ID4,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID4_Green-Light,ID4_Green-Light,ID2_Red,ID4_Green-Light,ID2_Red,ID4_Green-Light | chr16L > chr2L|chr6L > chr14L > chr2L|chr6L > chr14L > chr2L|chr6L | chr16L | donor_transfer |
| chr11L | ID3,ID5 | ID5_Blue-Dark,ID3_Orange | chr16L > chr14L | chr11R | subtelomere_switch |
| chr15R | ID4,ID2,ID2 | ID2_Red,ID2_Red,ID4_Green-Light | chr12R|chr14L|chr4R > chr2L|chr6L | chr7 | unmatched_array |
| chr12R | ID4,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID1 | ID1_Gray,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID4_Green-Light | chr12R > chr2L|chr6L | chr12R | tandem_amplification_same_end |
| chr10R | ID4,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID4_Green-Light | chr16L > chr2L|chr6L | chr16L | donor_transfer |
| chr3R | ID2,ID4,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID4_Green-Light,ID2_Red | chr16L > chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr13L | subtelomere_switch |
| chr14R | ID4,ID2,ID4,ID2,ID4,ID2,ID4,ID2,ID4,ID6 | ID6_Purple-Neutral,ID4_Green-Light,ID2_Red,ID4_Green-Light,ID2_Red,ID4_Green-Light,ID2_Red,ID4_Green-Light,ID2_Red,ID4_Green-Light | chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L | chr16 | unmatched_array |
| chr8L | ID1,ID2,ID4,ID2 | ID1_Gray,ID2_Red,ID4_Green-Light,ID2_Red | chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | ambiguous | unmatched_array |
| chr7R | ID4,ID3,ID4,ID5 | ID5_Blue-Dark,ID4_Green-Light,ID3_Orange,ID4_Green-Light | chr2L|chr6L > chr14L > chr2L|chr6L | ambiguous | unmatched_array |
| chr11R | ID4,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID4_Green-Light | chr16L > chr2L|chr6L | chr16L | donor_transfer |
| chr4L | ID4,ID2,ID5 | ID5_Blue-Light,ID2_Red,ID4_Green-Light | chr14L > chr2L|chr6L | chr14L | donor_transfer |
| chr13L | ID2,ID2,ID4 | ID4_Green-Light,ID2_Red,ID2_Red | chr2L|chr6L > chr16L|chr4R | chr12R | donor_transfer_candidates:3 |
| chr1L | ID6,ID5 | ID5_Blue-Dark,ID6_Purple-Dark | chr16L|chr7R > chr10L|chr9L | chr14 | unmatched_array |
| chr11L | ID2,ID2,ID2,ID6 | ID6_Purple-Dark,ID2_Red,ID2_Red,ID2_Red | chr10L|chr9L > chr12R | chr10L | subtelomere_switch |
| chr9R | ID3,ID2,ID6 | ID6_Purple-Neutral,ID2_Red,ID3_Orange | chr14R > chr14L | chr14L | donor_transfer |
| chr12L | ID4,ID2,ID4,ID2,ID1 | ID1_Gray,ID2_Red,ID4_Green-Light,ID2_Red,ID4_Green-Light | chr12R|chr14L|chr4R > chr2L|chr6L > chr12R|chr14L|chr4R > chr2L|chr6L | ambiguous | unmatched_array |
| chr7R | ID2,ID4,ID2,ID5 |  |  |  |  |
| chr5L | ID6,ID4,ID2,ID6 |  |  |  |  |
| chr15L | ID5,ID3,ID4 | ID5_Blue-Dark,ID3_Orange,ID4_Green-Light | chr16L > chr14L > chr2L|chr6L | ambiguous | unmatched_array |
| chr2L | ID6,ID2,ID2,ID2,ID2,ID2,ID4 | ID4_Green-Light,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID6_Purple-Dark | chr12R|chr4R > chr10L|chr9L | chr12R | donor_transfer_candidates:2 |
| chr7R | ID4,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID4_Green-Light | chr12R|chr14L|chr4R > chr2L|chr6L | ambiguous | unmatched_array |
| chr5R | ID4,ID2,ID1 | ID1_Gray,ID2_Red,ID4_Green-Light | chr12R|chr14L|chr4R > chr2L|chr6L | ambiguous | unmatched_array |
| chr1L | ID2,ID4 | ID4_Green-Light,ID2_Red | chr2L|chr6L > chr12R|chr14L|chr15R|chr16L|chr4R | chr14 | unmatched_array |
| chr11R | ID4,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID4_Green-Light | chr16L > chr2L|chr6L | chr16L | donor_transfer |

### reads we call a switch and the paper does not
| paper_chr_end | paper_ids | our_ids | our_blocks |
|---|---|---|---|
| chr2L | ID4,ID2,ID2,ID3 | ID4_Green-Light,ID2_Red,ID2_Red,ID3_Orange | chr4R > chr14L |
| chr2R | ID1,ID2,ID2 | ID1_Gray,ID2_Red,ID2_Red | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr5R | ID1,ID3,ID3,ID2,ID2 | ID1_Gray,ID3_Orange,ID3_Orange,ID2_Red,ID2_Red | chr14L > chr16L|chr4R |
| chr10R | ID1,ID2,ID2 | ID1_Gray,ID2_Red,ID2_Red | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr14L | ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID3_Orange,ID3_Orange | chr16L > chr14L |
| chr11L | ID5,ID2,ID3 | ID5_Blue-Dark,ID2_Red,ID3_Orange | chr7R > chr14L |
| chr4R | ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red | chr12R > chr4R |
| chr4L | ID2,ID2,ID2,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID2_Red,ID2_Red,ID2_Red | chr16L|chr7R > chr12R |
| chr13R | ID2,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID2_Red | chr7R > chr12R|chr4R |
| chr6L | ID2,ID2,ID2,ID1 | ID1_Gray,ID2_Red,ID2_Red,ID2_Red | chr12L|chr13L|chr16R|chr5R|chr8L|chr8R > chr4R |
| chr2L | ID2,ID2,ID3 | ID3_Orange,ID2_Red,ID2_Red | chr14L > chr12R |
| chr7L | ID5,ID2,ID2,ID2,ID2,ID2 | ID5_Blue-Dark,ID2_Red,ID2_Red,ID2_Red,ID2_Red,ID2_Red | chr7R > chr12R |
| chr14L | ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID3_Orange,ID3_Orange | chr16L > chr14L |
| chr11L | ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange | chr7R > chr14L |
| chr6L | ID2,ID2,ID1 | ID1_Gray,ID2_Red,ID2_Red | chr12L|chr12R|chr13L|chr16R|chr5R|chr8L|chr8R > chr16L|chr4R |
| chr14L | ID3,ID2,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red,ID3_Orange,ID3_Orange,ID2_Red,ID3_Orange | chr7R > chr14L |
| chr2R | ID2,ID2,ID2,ID3 | ID3_Orange,ID2_Red,ID2_Red,ID2_Red | chr14L > chr4R |


## strain 7302, PD 28  ->  7302_day3_with_selection (103)
* 103 of 103 paper reads found in our output; chr_end agrees 103/103; Y' copy count agrees 102/103
* template switching: paper flags 9; we catch 7, miss 2, and flag 7 reads the paper calls N (91.3 % agreement)

### the paper's switching reads, as our pipeline called them
| paper_chr_end | paper_ids | our_ids | our_blocks | our_source | our_mechanism |
|---|---|---|---|---|---|
| chr13L | ID7,ID2,ID8,ID2,ID8,ID2,ID8 | ID7_Yellow,ID2_Red-Dark,ID2_Red-Light,ID2_Red-Dark,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light | chr15R > chr12R|chr14L|chr4R > chr15R > chr13L | chr13L | tandem_amplification_same_end |
| chr13R | ID4,ID2 | ID4_Green-Light,ID2_Red-Light | chr2L > chr12R|chr13L|chr14L|chr4R | ambiguous | unmatched_array |
| chr14R | ID3,ID3,ID3,ID5 | ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange | chr14L | chr14L | donor_transfer |
| chr13R | ID7,ID5 | ID5_Blue-Dark,ID7_Yellow | chr16L|chr7R > chr13L | chr15R | subtelomere_switch |
| chr7R | ID5,ID5,ID3,ID5,ID5,ID3 | ID5_Blue-Dark,ID5_Blue-Dark,ID3_Orange,ID5_Blue-Dark,ID5_Blue-Dark,ID3_Orange | chr7R > chr14L > chr7R > chr14L | ambiguous | unmatched_array |
| chr13L | ID2,ID7,ID8,ID2,ID7,ID8,ID2,ID7,ID8,ID7 | ID7_Yellow,ID2_Red-Light,ID7_Yellow,ID2_Red-Light,ID2_Red-Light,ID7_Yellow,ID2_Red-Light,ID2_Red-Light,ID7_Yellow,ID2_Red-Light | chr13L | chr13L | tandem_amplification_same_end |
| chr15L | ID6,ID3,ID3 | ID6_Purple-Dark,ID3_Orange,ID3_Orange | chr10L|chr9L > chr14L | chr14L | donor_transfer |
| chr16R | ID2,ID3,ID3,ID1 | ID1_Gray,ID3_Orange,ID3_Orange,ID2_Red-Dark | chr14L > chr15R | chr14L | donor_transfer |
| chr14L | ID2,ID3,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID2_Red-Dark | chr16L|chr7R > chr14L > chr15R | chr14L | tandem_amplification_same_end |

### reads we call a switch and the paper does not
| paper_chr_end | paper_ids | our_ids | our_blocks |
|---|---|---|---|
| chr14L | ID3,ID3,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID2,ID2,ID2,ID2,ID2,ID2 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light | chr16L|chr7R > chr14L > chr12R |
| chr14L | ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr16L|chr7R > chr14L |
| chr12L | ID2,ID2,ID2,ID3 | ID3_Orange,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light | chr14L > chr12R|chr4R |
| chr14L | ID3,ID3,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr16L|chr7R > chr14L |
| chr16R | ID1,ID2,ID3 | ID1_Gray,ID2_Red-Light,ID3_Orange | chr13L > chr14L |
| chr13R | ID3,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr16L|chr7R > chr14L |


## strain 7302, PD 33  ->  7302_day4_with_selection (86)
* 83 of 86 paper reads found in our output; chr_end agrees 83/83; Y' copy count agrees 80/83
* template switching: paper flags 7; we catch 6, miss 1, and flag 11 reads the paper calls N (85.5 % agreement)

### the paper's switching reads, as our pipeline called them
| paper_chr_end | paper_ids | our_ids | our_blocks | our_source | our_mechanism |
|---|---|---|---|---|---|
| chr3R | ID6,ID2,ID3,ID3 | ID6_Purple-Neutral,ID2_Red-Light,ID3_Orange,ID3_Orange | chr14R > chr14L | chr14L | donor_transfer |
| chr13L | ID2,ID2,ID2,ID2,ID2,ID7,ID8,ID8,ID7,ID8,ID7 | ID7_Yellow,ID2_Red-Light,ID7_Yellow,ID2_Red-Light,ID2_Red-Light,ID7_Yellow,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light | chr13L | chr13L | tandem_amplification_same_end |
| chr16R | ID8,ID7,ID3 | ID3_Orange,ID7_Yellow,ID2_Red-Light | chr14L > chr13L | chr13L | donor_transfer |
| chr15R | ID2,ID2,ID8,ID8 | ID2_Red-Dark,ID2_Red-Dark,ID2_Red-Light,ID2_Red-Light | chr15R > chr13L | chr12R | donor_transfer_candidates:2 |
| chr16R | ID3,ID7,ID8 | ID3_Orange,ID7_Yellow,ID2_Red-Light | chr14L > chr13L | chr13L | donor_transfer |
| chr1R | ID5,ID7,ID3 | ID5_Blue-Dark,ID7_Yellow,ID3_Orange | chr16L|chr7R > chr13L > chr14L | chr14 | unmatched_array |
| chr3R | ID5,ID8,ID7,ID8 | ID5_Blue-Dark,ID2_Red-Light,ID7_Yellow,ID2_Red-Light | chr16L|chr7R > chr13L | chr7R | subtelomere_switch |

### reads we call a switch and the paper does not
| paper_chr_end | paper_ids | our_ids | our_blocks |
|---|---|---|---|
| chr2L | ID2,ID2,ID3 | ID3_Orange,ID2_Red-Light,ID2_Red-Light | chr14L > chr12R|chr4R |
| chr14L | ID2,ID3,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID2_Red-Light | chr16L|chr7R > chr14L |
| chr14L | ID2,ID5,ID5,ID2,ID5,ID5,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID5_Blue-Light,ID5_Blue-Light,ID2_Red-Light,ID5_Blue-Light,ID5_Blue-Light,ID2_Red-Light | chr16L|chr7R > chr14L |
| chr14L | ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr16L|chr7R > chr14L |
| chr14L | ID3,ID3,ID3,ID3,ID3,ID3,ID2,ID5 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr16L|chr7R > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID2,ID3 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID2_Red-Light,ID3_Orange | chr16L|chr7R > chr14L |
| chr12R | ID3,ID3,ID2,ID2,ID2,ID2,ID2,ID2,ID1 | ID1_Gray,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID3_Orange,ID3_Orange | chr13L > chr14L |
| chr14L | ID5,ID2,ID3,ID3,ID3,ID3,ID3,ID3,ID3,ID3 | ID5_Blue-Dark,ID2_Red-Light,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange,ID3_Orange | chr16L|chr7R > chr14L |
| chr2R | ID2,ID2,ID2,ID2,ID3 | ID3_Orange,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light,ID2_Red-Light | chr14L > chr12R |
| chr15R | ID2,ID2,ID2 | ID2_Red-Dark,ID2_Red-Dark,ID2_Red-Light | chr15R > chr12R|chr13L|chr14L|chr4R |
| chr8L | ID8,ID8,ID7 | ID7_Yellow,ID2_Red-Light,ID2_Red-Light | chr13L > chr12R|chr4R |

