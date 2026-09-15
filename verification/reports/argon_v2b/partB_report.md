# Part B recombination verification — snapshot: verification/snapshot_argon_v2b

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
| chr11R | 519 | 7 | 1.3 | 7 | 1.3 | 0 | 0.0 | 0 | 2 | 0 | 0 | 2 | 5 | 0 | 1 | chr3L:3,chr4R:1,chr8:1,chr1L:1 | FAIL |
| chr13L | 292 | 15 | 5.1 | 3 | 1.0 | 12 | 4.1 | 7 | 2 | 0 | 1 | 0 | 0 | 12 | 0 | chr13L:3 | FAIL |
| chr13R | 432 | 11 | 2.5 | 11 | 2.5 | 0 | 0.0 | 0 | 10 | 0 | 0 | 0 | 3 | 5 | 0 | ambiguous:5,chr15R:2,chr2L:2,chr13L:1 | FAIL |

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


## B1 null (day-0 self-run): 7172_day0_with_selection

Criterion: non-Loss recombination <= 1.0% per end. Y' Loss is reported as the per-end baseline (subclonal copy-number heterogeneity), not as a failure.
**12 of 32 ends FAIL**; total reads 9317; arm-less sources: 8; ambiguous: 39.

### Failing ends
| chr_end | n_reads | n_recomb | pct_recomb | n_recomb_nonloss | pct_recomb_nonloss | n_loss | pct_loss | n_loss_confirmed_end | n_gain | n_first_change | n_yp_recomb | n_spacer_switch | n_x_switch | n_ambiguous | n_armless_source | nonloss_sources | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr2R | 327 | 4 | 1.2 | 4 | 1.2 | 0 | 0.0 | 0 | 2 | 0 | 0 | 0 | 3 | 0 | 1 | chr7R:1,chr7:1,chr3R:1,chr10R:1 | FAIL |
| chr4R | 64 | 9 | 14.1 | 4 | 6.2 | 5 | 7.8 | 4 | 4 | 0 | 0 | 0 | 0 | 5 | 0 | chr4R:4 | FAIL |
| chr5R | 345 | 4 | 1.2 | 4 | 1.2 | 0 | 0.0 | 0 | 2 | 1 | 0 | 0 | 1 | 1 | 0 | chr16L:1,chr5R:1,chr3L:1,ambiguous:1 | FAIL |
| chr7R | 233 | 3 | 1.3 | 3 | 1.3 | 0 | 0.0 | 0 | 3 | 0 | 0 | 0 | 0 | 1 | 0 | chr12R:2,ambiguous:1 | FAIL |
| chr8L | 393 | 5 | 1.3 | 4 | 1.0 | 1 | 0.3 | 1 | 4 | 0 | 0 | 1 | 1 | 3 | 0 | ambiguous:3,chr12R:1 | FAIL |
| chr8R | 126 | 5 | 4.0 | 2 | 1.6 | 3 | 2.4 | 1 | 2 | 0 | 0 | 0 | 0 | 3 | 1 | chr12R:1,chr16:1 | FAIL |
| chr10R | 271 | 4 | 1.5 | 4 | 1.5 | 0 | 0.0 | 0 | 2 | 0 | 0 | 0 | 2 | 0 | 0 | chr4L:2,chr7R:1,chr13L:1 | FAIL |
| chr12R | 51 | 51 | 100.0 | 45 | 88.2 | 6 | 11.8 | 4 | 0 | 0 | 0 | 51 | 0 | 0 | 0 | chr3R:45 | FAIL |
| chr14L | 144 | 13 | 9.0 | 6 | 4.2 | 7 | 4.9 | 6 | 5 | 0 | 1 | 0 | 0 | 7 | 1 | chr14L:5,chr12:1 | FAIL |
| chr15R | 226 | 7 | 3.1 | 7 | 3.1 | 0 | 0.0 | 0 | 7 | 0 | 0 | 0 | 0 | 0 | 0 | chr12R:3,chr16L:1,chr7R:1,chr15R:1 | FAIL |
| chr16L | 139 | 13 | 9.4 | 4 | 2.9 | 9 | 6.5 | 5 | 4 | 0 | 0 | 0 | 0 | 8 | 1 | chr14L:2,chr12R:1,chr16L:1 | FAIL |
| chr16R | 345 | 4 | 1.2 | 4 | 1.2 | 0 | 0.0 | 0 | 2 | 2 | 0 | 0 | 0 | 2 | 0 | ambiguous:2,chr14L:1,chr12R:1 | FAIL |

### Ends with Y' Loss baseline > 0
| chr_end | n_reads | n_loss | pct_loss | n_loss_confirmed_end |
|---|---|---|---|---|
| chr4R | 64 | 5 | 7.8 | 4 |
| chr5L | 370 | 1 | 0.3 | 0 |
| chr6L | 414 | 5 | 1.2 | 2 |
| chr8L | 393 | 1 | 0.3 | 1 |
| chr8R | 126 | 3 | 2.4 | 1 |
| chr12R | 51 | 6 | 11.8 | 4 |
| chr13L | 339 | 1 | 0.3 | 0 |
| chr14L | 144 | 7 | 4.9 | 6 |
| chr14R | 305 | 3 | 1.0 | 1 |
| chr16L | 139 | 9 | 6.5 | 5 |


## B2 replicate concordance: 7172_day4_with_selection vs 7172_day4_with_selection_repeat

Criterion: |delta pct_recombination| <= 5.0 points per end.
**9 of 32 ends flagged**.

| chr_end | n_a | n_b | pct_recomb_a | pct_recomb_b | delta | pct_gain_a | pct_gain_b | pct_loss_a | pct_loss_b | top_source_a | top_source_b | source_agrees | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr1L | 300 | 125 | 16.3 | 9.6 | -6.7 | 7.3 | 6.4 | 0.0 | 0.0 | chr1R | chr1R | True | FLAG |
| chr1R | 101 | 56 | 8.9 | 14.3 | 5.4 | 4.0 | 7.1 | 0.0 | 0.0 | chr1L | chr1L | True | FLAG |
| chr6L | 345 | 260 | 6.7 | 100.0 | 93.3 | 1.2 | 0.0 | 2.3 | 0.8 | chr14 | chr5R | False | FLAG |
| chr7L | 337 | 190 | 1.2 | 6.3 | 5.1 | 1.2 | 5.8 | 0.0 | 0.0 | chr8 | chr10L | False | FLAG |
| chr8L | 364 | 231 | 0.3 | 7.4 | 7.1 | 0.3 | 3.0 | 0.0 | 0.9 | chr5 | chr2L | False | FLAG |
| chr10R | 254 | 176 | 2.0 | 9.7 | 7.7 | 1.6 | 5.7 | 0.0 | 0.0 | chr13L | chr4L | False | FLAG |
| chr11L | 301 | 215 | 100.0 | 3.3 | -96.7 | 5.6 | 3.3 | 0.0 | 0.0 | chr11R | chr12L | False | FLAG |
| chr13R | 262 | 152 | 11.5 | 3.9 | -7.5 | 10.7 | 3.9 | 0.0 | 0.0 | chr15R | chr15R | True | FLAG |
| chr14L | 94 | 94 | 31.9 | 97.9 | 66.0 | 7.4 | 0.0 | 20.2 | 28.7 | chr14L | chr2L | False | FLAG |


## B2 replicate concordance: 6991_day0_with_selection vs 6991_day0_with_selection_repeat

Criterion: |delta pct_recombination| <= 5.0 points per end.
**1 of 32 ends flagged**.

| chr_end | n_a | n_b | pct_recomb_a | pct_recomb_b | delta | pct_gain_a | pct_gain_b | pct_loss_a | pct_loss_b | top_source_a | top_source_b | source_agrees | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr12R | 309 | 121 | 26.9 | 38.8 | 12.0 | 12.3 | 11.6 | 14.2 | 27.3 | chr4 | chr4 | True | FLAG |


## B2 replicate concordance: 6991_day0_with_selection_repeat vs 6991_day0_with_selection_repeat2

Criterion: |delta pct_recombination| <= 5.0 points per end.
**1 of 32 ends flagged**.

| chr_end | n_a | n_b | pct_recomb_a | pct_recomb_b | delta | pct_gain_a | pct_gain_b | pct_loss_a | pct_loss_b | top_source_a | top_source_b | source_agrees | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr12R | 121 | 28 | 38.8 | 28.6 | -10.3 | 11.6 | 7.1 | 27.3 | 21.4 | chr4 | chr3R | False | FLAG |


## B3 positive control: chr11L -> chr11R

Criterion: >= 90.0% of chr11L reads attributed to chr11R.

| sample | chr_end | n_reads | n_expected_source | pct_expected_source | mean_confidence_expected | pct_complex | sources | verdict |
|---|---|---|---|---|---|---|---|---|
| 7172_day4_with_selection | chr11L | 301 | 292 | 97.0 | 0.933 | 1.0 | chr11R:292,chr10L:1,chr14R:1,chr2L:1,chr13R:1 | PASS |
| 7172_day6_with_selection | chr11L | 133 | 123 | 92.5 | 0.932 | 3.8 | chr11R:123,chr10L:4,chr14R:3,chr4L:2,chr1R:1 | PASS |
| 7172_day9_with_selection | chr11L | 105 | 101 | 96.2 | 0.9 | 2.9 | chr11R:101,chr10L:1,chr14R:1,chr4R:1,chr1R:1 | PASS |


## B4 Y' Loss classification: 7302_day0_with_selection

real_contraction = telomere end confirmed (adapter after telomere, repeat >= 30 bp), probe count agrees with RepeatMasker, no long non-Y' tail; rm_miss = probe count > RepeatMasker count (library/RepeatMasker miss); unconfirmed_end = read may be truncated; other = >500 bp non-Y' sequence between Y' and telomere.

| chr_end | n_reads | ref_yprimes | n_loss | pct_loss | real_contraction | rm_miss | unconfirmed_end | other | pct_real | pct_rm_miss | rm_miss_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chr2L | 424 | 1 | 1 | 0.2 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr4R | 113 | 7 | 12 | 10.6 | 10 | 0 | 2 | 0 | 83.3 | 0.0 |  |
| chr5R | 582 | 1 | 1 | 0.2 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr6L | 595 | 1 | 13 | 2.2 | 9 | 0 | 4 | 0 | 69.2 | 0.0 |  |
| chr8R | 228 | 1 | 1 | 0.4 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr12L | 513 | 1 | 1 | 0.2 | 1 | 0 | 0 | 0 | 100.0 | 0.0 |  |
| chr12R | 123 | 6 | 8 | 6.5 | 8 | 0 | 0 | 0 | 100.0 | 0.0 |  |
| chr13L | 292 | 4 | 12 | 4.1 | 7 | 0 | 5 | 0 | 58.3 | 0.0 |  |
| chr14L | 133 | 5 | 13 | 9.8 | 13 | 0 | 0 | 0 | 100.0 | 0.0 |  |
| chr14R | 536 | 1 | 1 | 0.2 | 1 | 0 | 0 | 0 | 100.0 | 0.0 |  |


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


## B4 Y' Loss classification: 7302_day6_with_selection_repeat

real_contraction = telomere end confirmed (adapter after telomere, repeat >= 30 bp), probe count agrees with RepeatMasker, no long non-Y' tail; rm_miss = probe count > RepeatMasker count (library/RepeatMasker miss); unconfirmed_end = read may be truncated; other = >500 bp non-Y' sequence between Y' and telomere.

| chr_end | n_reads | ref_yprimes | n_loss | pct_loss | real_contraction | rm_miss | unconfirmed_end | other | pct_real | pct_rm_miss | rm_miss_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chr4R | 43 | 7 | 43 | 100.0 | 38 | 0 | 5 | 0 | 88.4 | 0.0 |  |
| chr5L | 161 | 1 | 1 | 0.6 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr6L | 171 | 1 | 2 | 1.2 | 1 | 0 | 1 | 0 | 50.0 | 0.0 |  |
| chr10L | 98 | 1 | 1 | 1.0 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr12R | 16 | 6 | 2 | 12.5 | 1 | 0 | 1 | 0 | 50.0 | 0.0 |  |
| chr13L | 98 | 4 | 95 | 96.9 | 85 | 0 | 10 | 0 | 89.5 | 0.0 |  |
| chr14L | 29 | 5 | 5 | 17.2 | 3 | 0 | 2 | 0 | 60.0 | 0.0 |  |
| chr14R | 144 | 1 | 2 | 1.4 | 1 | 0 | 1 | 0 | 50.0 | 0.0 |  |


## B4 Y' Loss classification: 7172_day0_with_selection

real_contraction = telomere end confirmed (adapter after telomere, repeat >= 30 bp), probe count agrees with RepeatMasker, no long non-Y' tail; rm_miss = probe count > RepeatMasker count (library/RepeatMasker miss); unconfirmed_end = read may be truncated; other = >500 bp non-Y' sequence between Y' and telomere.

| chr_end | n_reads | ref_yprimes | n_loss | pct_loss | real_contraction | rm_miss | unconfirmed_end | other | pct_real | pct_rm_miss | rm_miss_flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chr4R | 64 | 7 | 5 | 7.8 | 4 | 0 | 1 | 0 | 80.0 | 0.0 |  |
| chr5L | 370 | 1 | 1 | 0.3 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr6L | 414 | 1 | 5 | 1.2 | 2 | 0 | 3 | 0 | 40.0 | 0.0 |  |
| chr8L | 393 | 1 | 1 | 0.3 | 1 | 0 | 0 | 0 | 100.0 | 0.0 |  |
| chr8R | 126 | 1 | 3 | 2.4 | 1 | 0 | 2 | 0 | 33.3 | 0.0 |  |
| chr12R | 51 | 8 | 6 | 11.8 | 4 | 0 | 2 | 0 | 66.7 | 0.0 |  |
| chr13L | 339 | 1 | 1 | 0.3 | 0 | 0 | 1 | 0 | 0.0 | 0.0 |  |
| chr14L | 144 | 3 | 7 | 4.9 | 6 | 0 | 1 | 0 | 85.7 | 0.0 |  |
| chr14R | 305 | 1 | 3 | 1.0 | 1 | 0 | 2 | 0 | 33.3 | 0.0 |  |
| chr16L | 139 | 3 | 9 | 6.5 | 5 | 0 | 4 | 0 | 55.6 | 0.0 |  |


## B5 truth set: alternating ID2/ID1 fingerprint of chr13L

Reference array at chr13L: `ID2,ID1,ID2,ID1`. Longest alternating window at any other end: 2; ends whose own array holds both IDs are excluded (chr14L).
Truth set: reads at other ends with a gain-like Y' status carrying an alternating window >= 3. **Strict set (window >= 4, a full period): 39 reads, 92.3% attributed to chr13L (mean confidence of hits 0.64) -> PASS** (criterion >= 80.0%). All windows >= 3: 68 reads, 85.3%.

### By window length
| window_len | n | n_hit |
|---|---|---|
| 3 | 29 | 22 |
| 4 | 31 | 28 |
| 5 | 8 | 8 |

### Per sample
| sample | n | n_hit | pct_hit |
|---|---|---|---|
| 7302_day4_with_selection | 21 | 19 | 90.5 |
| 7302_day5_with_selection | 40 | 35 | 87.5 |
| 7302_day6_with_selection_repeat | 2 | 2 | 100.0 |
| 7302_day9_with_selection_repeat | 5 | 2 | 40.0 |

### Attributed sources
| source | n |
|---|---|
| chr13L | 58 |
| chr14L | 8 |
| chr7R | 1 |
| chr1L | 1 |

### Reads
| sample | read_id | chr_end | observed_array | window_len | recombination_source | overall_confidence | mechanism |
|---|---|---|---|---|---|---|---|
| 7302_day5_with_selection | SRR33298447.322782 | chr10R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.490537 | chr11L | ID2,ID1,ID2,ID1,ID2,ID2,ID2,ID2,ID2,ID2,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.149279 | chr11L | ID1,ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.146060 | chr11L | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.268817 | chr11R | ID2,ID1,ID2,ID2,ID1,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.12915 | chr11R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.28086 | chr12L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.170876 | chr13R | ID2,ID1,ID2,ID1,ID1,ID2,ID1 | 4 | chr13L | 0.33 | donor_transfer_candidates:2 |
| 7302_day5_with_selection | SRR33298447.85642 | chr14R | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.151294 | chr15R | ID1,ID2,ID1,ID2,ID1,ID1 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.78734 | chr15R | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.335886 | chr16L | ID1,ID2,ID1,ID1 | 3 | chr13L | 0.33 | donor_transfer_candidates:2 |
| 7302_day5_with_selection | SRR33298447.108657 | chr16L | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.502469 | chr16R | ID7,ID1,ID1,ID1,ID2,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.328635 | chr1L | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.148608 | chr1R | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.118629 | chr1R | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.151907 | chr2L | ID4,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.318819 | chr2R | ID2,ID1,ID2,ID1,ID1,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.149253 | chr2R | ID2,ID1,ID2,ID2,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.129927 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.74716 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.460402 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.33469 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day5_with_selection | SRR33298447.9362 | chr3R | ID3,ID1,ID2,ID1,ID1,ID2,ID1,ID2 | 4 | chr13L | 0.33 | donor_transfer_candidates:2 |
| 7302_day5_with_selection | SRR33298447.180074 | chr3R | ID2,ID1,ID2,ID2,ID1,ID2 | 3 | chr13L | 0.33 | donor_transfer_candidates:2 |
| 7302_day5_with_selection | SRR33298447.490300 | chr3R | ID2,ID1,ID2,ID1,ID1 | 4 | chr13L | 0.66 | donor_transfer |
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
| 7302_day4_with_selection | SRR33298449.497101 | chr15L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.293351 | chr16L | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.250891 | chr16L | ID1,ID2,ID1 | 3 | chr13L | 0.44 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.444426 | chr2R | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.196678 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.507268 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.564793 | chr2R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.567566 | chr3R | ID1,ID1,ID2,ID1 | 3 | chr7R | 0.7 | subtelomere_switch |
| 7302_day4_with_selection | SRR33298449.486235 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.43071 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.277608 | chr3R | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.4646 | chr3R | ID2,ID1,ID2 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.23058 | chr4L | ID2,ID1,ID2,ID1,ID2 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.141373 | chr4L | ID2,ID1,ID2,ID1 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.131857 | chr5R | ID1,ID1,ID2,ID1,ID2,ID1 | 5 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.463689 | chr5R | ID1,ID2,ID1,ID1 | 3 | chr13L | 0.33 | donor_transfer_candidates:2 |
| 7302_day4_with_selection | SRR33298449.347234 | chr7R | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day4_with_selection | SRR33298449.17944 | chr8L | ID1,ID2,ID1,ID2 | 4 | chr13L | 0.66 | donor_transfer |
| 7302_day9_with_selection_repeat | SRR33298445.216902 | chr15R | ID1,ID1,ID2,ID1 | 3 | chr13L | 0.66 | donor_transfer |
| 7302_day9_with_selection_repeat | SRR33298445.87017 | chr1R | ID1,ID2,ID1,ID1,ID2 | 3 | chr1L | 0.517 | subtelomere_switch |
| 7302_day9_with_selection_repeat | SRR33298445.127200 | chr2R | ID2,ID1,ID2,ID2,ID2,ID2 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day9_with_selection_repeat | SRR33298445.36153 | chr2R | ID2,ID1,ID1,ID2,ID1 | 3 | chr14L | 0.66 | donor_transfer |
| 7302_day9_with_selection_repeat | SRR33298445.90156 | chr5L | ID8,ID2,ID1,ID2,ID1,ID1 | 4 | chr13L | 0.66 | donor_transfer |

