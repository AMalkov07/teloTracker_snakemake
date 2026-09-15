
## B3 positive control: chr11L -> chr11R

Criterion: >= 90.0% of chr11L reads attributed to chr11R.

| sample | chr_end | n_reads | n_expected_source | pct_expected_source | mean_confidence_expected | pct_complex | sources | verdict |
|---|---|---|---|---|---|---|---|---|
| 7172_day4_with_selection | chr11L | 301 | 292 | 97.0 | 0.933 | 1.0 | chr11R:292,chr10L:1,chr14R:1,chr2L:1,chr13R:1 | PASS |
| 7172_day6_with_selection | chr11L | 133 | 123 | 92.5 | 0.924 | 5.3 | chr11R:123,chr10L:4,chr14R:3,chr4L:2,chr1R:1 | PASS |
| 7172_day9_with_selection | chr11L | 105 | 101 | 96.2 | 0.88 | 7.6 | chr11R:101,chr10L:1,chr14R:1,chr4R:1,chr1R:1 | PASS |


## B5 truth set: alternating ID2/ID1 fingerprint of chr13L

Reference array at chr13L: `ID2,ID1,ID2,ID1`. Longest alternating window at any other end: 2; ends whose own array holds both IDs are excluded (chr14L).
Truth set: reads at other ends with a gain-like Y' status carrying an alternating window >= 3. **Strict set (window >= 4, a full period): 38 reads, 92.1% attributed to chr13L (mean confidence of hits 0.51) -> PASS** (criterion >= 80.0%). All windows >= 3: 63 reads, 88.9%.

### By window length
| window_len | n | n_hit |
|---|---|---|
| 3 | 25 | 21 |
| 4 | 30 | 27 |
| 5 | 8 | 8 |

### Per sample
| sample | n | n_hit | pct_hit |
|---|---|---|---|
| 7302_day4_with_selection | 21 | 19 | 90.5 |
| 7302_day5_with_selection | 40 | 35 | 87.5 |
| 7302_day6_with_selection_repeat | 2 | 2 | 100.0 |

### Attributed sources
| source | n |
|---|---|
| chr13L | 56 |
| chr14L | 6 |
| chr7R | 1 |

### Reads
| sample | read_id | chr_end | observed_array | window_len | recombination_source | overall_confidence | mechanism |
|---|---|---|---|---|---|---|---|
| 7302_day5_with_selection | SRR33298447.322782 | chr10R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.490537 | chr11L | ID2,ID1,ID2,ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | 5 | chr13L | -0.0026 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.149279 | chr11L | ID1,ID1,ID2,ID1,ID2 | 4 | chr13L | 0.0036 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.146060 | chr11L | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.268817 | chr11R | ID2,ID1,ID2,ID2,ID1,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.12915 | chr11R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.28086 | chr12L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.170876 | chr13R | ID2,ID1,ID2,ID1,ID1,ID2,ID1 | 4 | chr13L | 0.0008 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.85642 | chr14R | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.151294 | chr15R | ID1,ID2,ID1,ID2,ID1,ID1 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.78734 | chr15R | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.335886 | chr16L | ID1,ID2,ID1,ID1 | 3 | chr13L | 0.33 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.108657 | chr16L | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.502469 | chr16R | ID7,ID1,ID1,ID1,ID2,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.328635 | chr1L | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.148608 | chr1R | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.118629 | chr1R | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.151907 | chr2L | ID4,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.318819 | chr2R | ID2,ID1,ID2,ID1,ID1,ID1 | 4 | chr13L | 0.008 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.149253 | chr2R | ID2,ID1,ID2,ID2,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.129927 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.74716 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.460402 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.33469 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.9362 | chr3R | ID3,ID1,ID2,ID1,ID1,ID2,ID1,ID2 | 4 | chr13L | 0.0032 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.180074 | chr3R | ID2,ID1,ID2,ID2,ID1,ID2 | 3 | chr13L | 0.33 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.490300 | chr3R | ID2,ID1,ID2,ID1,ID1 | 4 | chr13L | -0.0007 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.207200 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.168241 | chr4L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.434153 | chr4L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.295395 | chr5L | ID8,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.16244 | chr5L | ID8,ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.318447 | chr5R | ID1,ID2,ID1,ID2,ID2,ID2,ID2,ID2,ID4 | 4 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.325287 | chr5R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.498649 | chr5R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.22803 | chr7R | ID1,ID2,ID1,ID2,ID1 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.508772 | chr7R | ID1,ID2,ID1,ID2,ID2 | 4 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.150844 | chr7R | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.329763 | chr7R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.18466 | chr7R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day6_with_selection_repeat | SRR33298446.373861 | chr5R | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day6_with_selection_repeat | SRR33298446.252236 | chr5R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.327615 | chr10R | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.274261 | chr10R | ID2,ID1,ID2,ID1 | 4 | chr14L | 0.3559 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.368231 | chr14R | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.497101 | chr15L | ID2,ID1,ID2,ID1 | 4 | chr13L | -0.0009 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.293351 | chr16L | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.250891 | chr16L | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.444426 | chr2R | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.196678 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.507268 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.0002 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.564793 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.567566 | chr3R | ID1,ID1,ID2,ID1 | 3 | chr7R | 0.7 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.486235 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.43071 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.277608 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.4646 | chr3R | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.23058 | chr4L | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.141373 | chr4L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.131857 | chr5R | ID1,ID1,ID2,ID1,ID2,ID1 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.463689 | chr5R | ID1,ID2,ID1,ID1 | 3 | chr13L | 0.33 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.347234 | chr7R | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.17944 | chr8L | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |


## B1 null (day-0 self-run): 7302_day0_with_selection

Criterion: non-Loss recombination <= 1.0% per end. Y' Loss is reported as the per-end baseline (subclonal copy-number heterogeneity), not as a failure.
**7 of 32 ends FAIL**; total reads 14202; arm-less sources: 4; ambiguous: 75.

### Failing ends
| chr_end | n_reads | n_recomb | pct_recomb | n_recomb_nonloss | pct_recomb_nonloss | n_loss | pct_loss | n_loss_confirmed_end | n_gain | n_first_change | n_yp_recomb | n_spacer_switch | n_x_switch | n_ambiguous | n_armless_source | nonloss_sources | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr1R | 213 | 4 | 1.9 | 4 | 1.9 | 0 | 0.0 | 0 | 3 | 0 | 0 | 0 | 1 | 0 | 0 | chr14R:3,chr1L:1 | FAIL |
| chr4L | 273 | 3 | 1.1 | 3 | 1.1 | 0 | 0.0 | 0 | 2 | 0 | 0 | 0 | 3 | 0 | 0 | chr10L:2,chr10R:1 | FAIL |
| chr6L | 595 | 27 | 4.5 | 14 | 2.4 | 13 | 2.2 | 9 | 7 | 7 | 0 | 0 | 1 | 18 | 3 | ambiguous:6,chr13:3,chr4R:1,chr14L:1 | FAIL |
| chr10R | 395 | 5 | 1.3 | 5 | 1.3 | 0 | 0.0 | 0 | 1 | 0 | 0 | 0 | 4 | 0 | 0 | chr4L:4,chr14L:1 | FAIL |
| chr11R | 519 | 7 | 1.3 | 7 | 1.3 | 0 | 0.0 | 0 | 2 | 0 | 0 | 3 | 5 | 0 | 1 | chr3L:3,chr4R:1,chr8:1,chr1L:1 | FAIL |
| chr13L | 292 | 15 | 5.1 | 3 | 1.0 | 12 | 4.1 | 7 | 2 | 0 | 1 | 0 | 0 | 12 | 0 | chr13L:3 | FAIL |
| chr13R | 432 | 11 | 2.5 | 11 | 2.5 | 0 | 0.0 | 0 | 10 | 0 | 0 | 2 | 3 | 5 | 0 | ambiguous:5,chr15R:2,chr2L:2,chr13L:1 | FAIL |

### Ends with Y' Loss baseline > 0
| chr_end | n_reads | n_loss | pct_loss | n_loss_confirmed_end |
|---|---|---|---|---|
| chr2L | 424 | 1 | 0.2 | 0 |
| chr4R | 113 | 12 | 10.6 | 10 |
| chr5R | 582 | 1 | 0.2 | 0 |
| chr6L | 595 | 13 | 2.2 | 9 |
| chr8R | 228 | 1 | 0.4 | 0 |
| chr12L | 513 | 1 | 0.2 | 1 |
| chr12R | 123 | 8 | 6.5 | 8 |
| chr13L | 292 | 12 | 4.1 | 7 |
| chr14L | 133 | 13 | 9.8 | 13 |
| chr14R | 536 | 1 | 0.2 | 1 |


## B5 truth set: alternating ID2/ID1 fingerprint of chr13L

Reference array at chr13L: `ID2,ID1,ID2,ID1`. Longest alternating window at any other end: 2; ends whose own array holds both IDs are excluded (chr14L).
Truth set: reads at other ends with a gain-like Y' status carrying an alternating window >= 3. **Strict set (window >= 4, a full period): 38 reads, 92.1% attributed to chr13L (mean confidence of hits 0.51) -> PASS** (criterion >= 80.0%). All windows >= 3: 63 reads, 88.9%.

### By window length
| window_len | n | n_hit |
|---|---|---|
| 3 | 25 | 21 |
| 4 | 30 | 27 |
| 5 | 8 | 8 |

### Per sample
| sample | n | n_hit | pct_hit |
|---|---|---|---|
| 7302_day4_with_selection | 21 | 19 | 90.5 |
| 7302_day5_with_selection | 40 | 35 | 87.5 |
| 7302_day6_with_selection_repeat | 2 | 2 | 100.0 |

### Attributed sources
| source | n |
|---|---|
| chr13L | 56 |
| chr14L | 6 |
| chr7R | 1 |

### Reads
| sample | read_id | chr_end | observed_array | window_len | recombination_source | overall_confidence | mechanism |
|---|---|---|---|---|---|---|---|
| 7302_day5_with_selection | SRR33298447.322782 | chr10R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.490537 | chr11L | ID2,ID1,ID2,ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | 5 | chr13L | -0.0026 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.149279 | chr11L | ID1,ID1,ID2,ID1,ID2 | 4 | chr13L | 0.0036 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.146060 | chr11L | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.268817 | chr11R | ID2,ID1,ID2,ID2,ID1,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.12915 | chr11R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.28086 | chr12L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.170876 | chr13R | ID2,ID1,ID2,ID1,ID1,ID2,ID1 | 4 | chr13L | 0.0008 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.85642 | chr14R | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.151294 | chr15R | ID1,ID2,ID1,ID2,ID1,ID1 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.78734 | chr15R | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.335886 | chr16L | ID1,ID2,ID1,ID1 | 3 | chr13L | 0.33 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.108657 | chr16L | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.502469 | chr16R | ID7,ID1,ID1,ID1,ID2,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.328635 | chr1L | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.148608 | chr1R | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.118629 | chr1R | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.151907 | chr2L | ID4,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.318819 | chr2R | ID2,ID1,ID2,ID1,ID1,ID1 | 4 | chr13L | 0.008 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.149253 | chr2R | ID2,ID1,ID2,ID2,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.129927 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.74716 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.460402 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.33469 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.9362 | chr3R | ID3,ID1,ID2,ID1,ID1,ID2,ID1,ID2 | 4 | chr13L | 0.0032 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.180074 | chr3R | ID2,ID1,ID2,ID2,ID1,ID2 | 3 | chr13L | 0.33 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.490300 | chr3R | ID2,ID1,ID2,ID1,ID1 | 4 | chr13L | -0.0007 | subtelomere_switch |
| 7302_day5_with_selection | SRR33298447.207200 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.168241 | chr4L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.434153 | chr4L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.295395 | chr5L | ID8,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.16244 | chr5L | ID8,ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.318447 | chr5R | ID1,ID2,ID1,ID2,ID2,ID2,ID2,ID2,ID4 | 4 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.325287 | chr5R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.498649 | chr5R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.22803 | chr7R | ID1,ID2,ID1,ID2,ID1 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.508772 | chr7R | ID1,ID2,ID1,ID2,ID2 | 4 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.150844 | chr7R | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.329763 | chr7R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.18466 | chr7R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day6_with_selection_repeat | SRR33298446.373861 | chr5R | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day6_with_selection_repeat | SRR33298446.252236 | chr5R | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.327615 | chr10R | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.274261 | chr10R | ID2,ID1,ID2,ID1 | 4 | chr14L | 0.3559 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.368231 | chr14R | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.497101 | chr15L | ID2,ID1,ID2,ID1 | 4 | chr13L | -0.0009 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.293351 | chr16L | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.250891 | chr16L | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.444426 | chr2R | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.196678 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.507268 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.0002 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.564793 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.567566 | chr3R | ID1,ID1,ID2,ID1 | 3 | chr7R | 0.7 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.486235 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.43071 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.277608 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.4646 | chr3R | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.23058 | chr4L | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.141373 | chr4L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.131857 | chr5R | ID1,ID1,ID2,ID1,ID2,ID1 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.463689 | chr5R | ID1,ID2,ID1,ID1 | 3 | chr13L | 0.33 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.347234 | chr7R | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.17944 | chr8L | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |

