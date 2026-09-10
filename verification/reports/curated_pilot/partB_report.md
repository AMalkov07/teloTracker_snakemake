
## B4 Y' Loss classification: 7302_day0_with_selection

real_contraction = telomere end confirmed (adapter after telomere, repeat >= 30 bp), probe count agrees with RepeatMasker, no long non-Y' tail; rm_miss = probe count > RepeatMasker count (library/RepeatMasker miss); unconfirmed_end = read may be truncated; other = >500 bp non-Y' sequence between Y' and telomere.

| chr_end | n_reads | ref_yprimes | n_loss | pct_loss | real_contraction | rm_miss | unconfirmed_end | other | pct_real | pct_rm_miss | rm_miss_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chr2L | 424 | 1 | 1 | 0.2 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr4R | 113 | 7 | 13 | 11.5 | 11 | 0 | 2 | 0 | 84.6 | 0.0 |  |
| chr5R | 582 | 1 | 1 | 0.2 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr6L | 595 | 1 | 13 | 2.2 | 9 | 0 | 4 | 0 | 69.2 | 0.0 |  |
| chr8R | 228 | 1 | 1 | 0.4 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr12L | 513 | 1 | 1 | 0.2 | 1 | 0 | 0 | 0 | 100.0 | 0.0 |  |
| chr12R | 123 | 6 | 8 | 6.5 | 8 | 0 | 0 | 0 | 100.0 | 0.0 |  |
| chr13L | 292 | 4 | 12 | 4.1 | 7 | 0 | 5 | 0 | 58.3 | 0.0 |  |
| chr14L | 133 | 5 | 13 | 9.8 | 13 | 0 | 0 | 0 | 100.0 | 0.0 |  |
| chr14R | 536 | 1 | 1 | 0.2 | 1 | 0 | 0 | 0 | 100.0 | 0.0 |  |


## B4 Y' Loss classification: 7302_day4_with_selection

real_contraction = telomere end confirmed (adapter after telomere, repeat >= 30 bp), probe count agrees with RepeatMasker, no long non-Y' tail; rm_miss = probe count > RepeatMasker count (library/RepeatMasker miss); unconfirmed_end = read may be truncated; other = >500 bp non-Y' sequence between Y' and telomere.

| chr_end | n_reads | ref_yprimes | n_loss | pct_loss | real_contraction | rm_miss | unconfirmed_end | other | pct_real | pct_rm_miss | rm_miss_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chr2L | 203 | 1 | 1 | 0.5 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr4R | 41 | 7 | 6 | 14.6 | 3 | 0 | 3 | 0 | 50.0 | 0.0 |  |
| chr8R | 118 | 1 | 1 | 0.8 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr12L | 250 | 1 | 2 | 0.8 | 2 | 0 | 0 | 0 | 100.0 | 0.0 |  |
| chr12R | 48 | 6 | 9 | 18.8 | 3 | 0 | 6 | 0 | 33.3 | 0.0 |  |
| chr13L | 123 | 4 | 15 | 12.2 | 13 | 0 | 2 | 0 | 86.7 | 0.0 |  |
| chr14L | 88 | 5 | 17 | 19.3 | 12 | 0 | 5 | 0 | 70.6 | 0.0 |  |


## B4 Y' Loss classification: 7302_day5_with_selection

real_contraction = telomere end confirmed (adapter after telomere, repeat >= 30 bp), probe count agrees with RepeatMasker, no long non-Y' tail; rm_miss = probe count > RepeatMasker count (library/RepeatMasker miss); unconfirmed_end = read may be truncated; other = >500 bp non-Y' sequence between Y' and telomere.

| chr_end | n_reads | ref_yprimes | n_loss | pct_loss | real_contraction | rm_miss | unconfirmed_end | other | pct_real | pct_rm_miss | rm_miss_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chr4R | 69 | 7 | 17 | 24.6 | 14 | 0 | 3 | 0 | 82.4 | 0.0 |  |
| chr6L | 266 | 1 | 3 | 1.1 | 1 | 0 | 2 | 0 | 33.3 | 0.0 |  |
| chr8L | 259 | 1 | 2 | 0.8 | 0 | 0 | 2 | 0 | 0.0 | 0.0 |  |
| chr8R | 125 | 1 | 2 | 1.6 | 0 | 0 | 2 | 0 | 0.0 | 0.0 |  |
| chr12R | 76 | 6 | 24 | 31.6 | 20 | 1 | 3 | 0 | 83.3 | 4.2 |  |
| chr13L | 151 | 4 | 22 | 14.6 | 13 | 0 | 9 | 0 | 59.1 | 0.0 |  |
| chr14L | 84 | 5 | 11 | 13.1 | 8 | 0 | 3 | 0 | 72.7 | 0.0 |  |
| chr14R | 248 | 1 | 2 | 0.8 | 2 | 0 | 0 | 0 | 100.0 | 0.0 |  |


## B1 null (day-0 self-run): 7302_day0_with_selection

Criterion: non-Loss recombination <= 1.0% per end. Y' Loss is reported as the per-end baseline (subclonal copy-number heterogeneity), not as a failure.
**8 of 32 ends FAIL**; total reads 14202; arm-less sources: 3; ambiguous: 69.

### Failing ends
| chr_end | n_reads | n_recomb | pct_recomb | n_recomb_nonloss | pct_recomb_nonloss | n_loss | pct_loss | n_loss_confirmed_end | n_gain | n_first_change | n_yp_recomb | n_spacer_switch | n_x_switch | n_ambiguous | n_armless_source | nonloss_sources | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr1R | 213 | 4 | 1.9 | 4 | 1.9 | 0 | 0.0 | 0 | 3 | 0 | 0 | 0 | 1 | 0 | 0 | chr14R:3,chr1L:1 | FAIL |
| chr4L | 273 | 3 | 1.1 | 3 | 1.1 | 0 | 0.0 | 0 | 2 | 0 | 0 | 0 | 3 | 0 | 0 | chr10L:2,chr10R:1 | FAIL |
| chr6L | 595 | 27 | 4.5 | 14 | 2.4 | 13 | 2.2 | 9 | 7 | 7 | 0 | 0 | 1 | 13 | 1 | chr14L:10,chr4R:1,chr12R:1,chr14:1 | FAIL |
| chr10R | 395 | 5 | 1.3 | 5 | 1.3 | 0 | 0.0 | 0 | 1 | 0 | 0 | 0 | 4 | 0 | 0 | chr4L:4,chr14L:1 | FAIL |
| chr11R | 519 | 7 | 1.3 | 7 | 1.3 | 0 | 0.0 | 0 | 2 | 0 | 0 | 2 | 5 | 0 | 1 | chr3L:3,chr14L:1,chr8:1,chr1L:1 | FAIL |
| chr13L | 292 | 15 | 5.1 | 3 | 1.0 | 12 | 4.1 | 7 | 2 | 0 | 1 | 0 | 0 | 12 | 0 | chr13L:3 | FAIL |
| chr13R | 432 | 11 | 2.5 | 11 | 2.5 | 0 | 0.0 | 0 | 10 | 0 | 0 | 1 | 3 | 3 | 0 | chr13L:3,ambiguous:3,chr2L:2,chr16L:1 | FAIL |
| chr15R | 363 | 56 | 15.4 | 56 | 15.4 | 0 | 0.0 | 0 | 3 | 53 | 0 | 0 | 0 | 0 | 0 | chr14L:53,chr15R:2,chr13L:1 | FAIL |

### Ends with Y' Loss baseline > 0
| chr_end | n_reads | n_loss | pct_loss | n_loss_confirmed_end |
|---|---|---|---|---|
| chr2L | 424 | 1 | 0.2 | 0 |
| chr4R | 113 | 13 | 11.5 | 11 |
| chr5R | 582 | 1 | 0.2 | 0 |
| chr6L | 595 | 13 | 2.2 | 9 |
| chr8R | 228 | 1 | 0.4 | 0 |
| chr12L | 513 | 1 | 0.2 | 1 |
| chr12R | 123 | 8 | 6.5 | 8 |
| chr13L | 292 | 12 | 4.1 | 7 |
| chr14L | 133 | 13 | 9.8 | 13 |
| chr14R | 536 | 1 | 0.2 | 1 |


## B5 truth set: alternating ID7/ID2 fingerprint of chr13L

Reference array at chr13L: `ID7,ID2,ID7,ID2`. Longest alternating window at any other end: 1; ends whose own array holds both IDs are excluded (none).
Truth set: reads at other ends with a gain-like Y' status carrying an alternating window >= 3. **Strict set (window >= 4, a full period): 32 reads, 96.9% attributed to chr13L (mean confidence of hits 0.66) -> PASS** (criterion >= 80.0%). All windows >= 3: 49 reads, 95.9%.

### By window length
| window_len | n | n_hit |
|---|---|---|
| 3 | 17 | 16 |
| 4 | 25 | 24 |
| 5 | 7 | 7 |

### Per sample
| sample | n | n_hit | pct_hit |
|---|---|---|---|
| 7302_day4_with_selection | 19 | 17 | 89.5 |
| 7302_day5_with_selection | 30 | 30 | 100.0 |

### Attributed sources
| source | n |
|---|---|
| chr13L | 47 |
| chr14L | 1 |
| chr7R | 1 |

### Reads
| sample | read_id | chr_end | observed_array | window_len | recombination_source | overall_confidence | mechanism |
|---|---|---|---|---|---|---|---|
| 7302_day4_with_selection | SRR33298449.327615 | chr10R | ID7,ID2,ID7,ID2,ID7 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.274261 | chr10R | ID7,ID2,ID7,ID2 | 4 | chr14L | 0.3559 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.368231 | chr14R | ID7,ID2,ID7 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.497101 | chr15L | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.293351 | chr16L | ID5,ID2,ID7,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.444426 | chr2R | ID7,ID2,ID7,ID2,ID7 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.196678 | chr2R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.507268 | chr2R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.564793 | chr2R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.567566 | chr3R | ID5,ID2,ID7,ID2 | 3 | chr7R | 0.7 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.486235 | chr3R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.43071 | chr3R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.277608 | chr3R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.4646 | chr3R | ID7,ID2,ID7 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.23058 | chr4L | ID7,ID2,ID7,ID2,ID7 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.141373 | chr4L | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.131857 | chr5R | ID1,ID2,ID7,ID2,ID7,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.347234 | chr7R | ID5,ID7,ID2,ID7 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.17944 | chr8L | ID2,ID7,ID2,ID7 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.322782 | chr10R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.490537 | chr11L | ID7,ID2,ID7,ID2,ID7,ID3,ID3,ID3,ID3,ID3,ID3 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.146060 | chr11L | ID7,ID2,ID7 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.268817 | chr11R | ID7,ID2,ID7,ID7,ID2,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.12915 | chr11R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.28086 | chr12L | ID1,ID2,ID7,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.170876 | chr13R | ID7,ID2,ID7,ID2,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.85642 | chr14R | ID2,ID7,ID2,ID7 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.151294 | chr15R | ID2,ID7,ID2,ID7,ID2,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.78734 | chr15R | ID2,ID2,ID7,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.328635 | chr1L | ID5,ID7,ID2,ID7 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.148608 | chr1R | ID7,ID2,ID7,ID2,ID7 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.118629 | chr1R | ID7,ID2,ID7 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.151907 | chr2L | ID4,ID2,ID7,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.318819 | chr2R | ID7,ID2,ID7,ID2,ID2,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.149253 | chr2R | ID7,ID2,ID7,ID7,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.129927 | chr2R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.74716 | chr2R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.460402 | chr2R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.33469 | chr2R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.9362 | chr3R | ID6,ID2,ID7,ID2,ID2,ID7,ID2,ID7 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.180074 | chr3R | ID7,ID2,ID7,ID7,ID2,ID7 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.490300 | chr3R | ID7,ID2,ID7,ID2,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.207200 | chr3R | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.168241 | chr4L | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.434153 | chr4L | ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.295395 | chr5L | ID6,ID2,ID7,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.318447 | chr5R | ID1,ID7,ID2,ID7,ID7,ID7,ID7,ID7,ID4 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.22803 | chr7R | ID5,ID7,ID2,ID7,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.508772 | chr7R | ID5,ID7,ID2,ID7,ID3 | 3 | chr13L | 0.66 | donor_transfer |

