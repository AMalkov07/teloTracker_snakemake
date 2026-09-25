# Recombination confidence: old vs new

Same reads, one pass: `verification/snapshot_v3conf` (28 samples, 234,899 reads). `overall_confidence` is the old score (still written, deprecated); `recombination_confidence` and `donor_confidence` are the new ones.

## 0. Nothing else changed

234,899 reads matched to `verification/snapshot_v2c_path`: old score identical for **100.000%**, donor call identical for **100.000%**.

## 1. Shape

Hardcoded values: 0.30 / 0.33 (old recombinant floor), 0.95 / 0.90 (old unchanged constants). The new unchanged-read score still takes two values by design: 0.95 when the read reaches the telomere, 0.70 when it stops early -- which, unlike before, is information.

| score | reads | n | distinct values | % on a hardcoded value | mean | median |
|---|---|---|---|---|---|---|
| old overall_confidence | recombinant | 33987 | 3372 | 59.326 | 0.459 | 0.330 |
| new recombination_confidence | recombinant | 33987 | 2731 | 0.000 | 0.720 | 0.769 |
| new donor_confidence | recombinant | 33987 | 3148 | 0.000 | 0.289 | 0.238 |
| old overall_confidence | unchanged | 200912 | 2 | 100.000 | 0.950 | 0.950 |
| new recombination_confidence | unchanged | 200912 | 2 | 77.191 | 0.893 | 0.950 |

## 2. Per-end summary columns

R^2 of each per-end mean against the recombination rate, over 883 (sample, end) pairs. High = the column mostly restates the rate printed next to it.

| column | R^2 |
|---|---|
| mean_confidence (old) | 0.824 |
| mean_recombination_confidence (new) | 0.180 |
| mean_donor_confidence (new) | 0.000 |

## 3. Does the donor score separate right from wrong donors?

| truth set | correct | wrong | AUC old | AUC new (donor) | new: mean correct | new: mean wrong |
|---|---|---|---|---|---|---|
| positive control 7172 chr11L -> chr11R | 516 | 23 | 0.829 | 0.810 | 0.551 | 0.481 |
| chr13L truth set, strict | 36 | 3 | 0.630 | 0.769 | 0.941 | 0.616 |

The chr13L strict set has very few wrong calls, so its AUC is noisy either way.

## 4. Does the call score separate real events from day-0 calls?

Day-0 calls include genuine standing variation, so no score can reach 1.0 here.

| score | AUC later vs day-0 | median day-0 | median later |
|---|---|---|---|
| old overall_confidence | 0.633 | 0.330 | 0.330 |
| new recombination_confidence | 0.657 | 0.500 | 0.800 |

## 5. Worked examples: before and after

### Y'-only donor transfer that sat on the old 0.30 floor

`7302_day3_with_selection` chr7L read `SRR33298450.273100`

- Y' array: `ID2,ID1,ID1,ID1,ID1,ID1,ID1,ID1` (Y' Gain; telomere reached: True)
- spacer: no_change (conf 0.95); X element: no_change (conf 0.6736)
- Y' fingerprint: `-` (length 3); path: `chr14L?[3]:ID2 > chr12R[2-6]:ID1,ID1,ID1,ID1,ID1,ID1,ID1(circ x1.4 strong)`
- donor named: **chr12R** (donor_transfer)

| | old | new |
|---|---|---|
| score | `overall_confidence` = **0.3** | `recombination_confidence` = **0.9961**, `donor_confidence` = **0.9375** |
| why | `base=0.30, complexity=1.00` | `call[yprime=1.00] donor[support=0.94;margin=1.00]` |

### Spacer switch + Y' gain, where a weak spacer switch hid the Y' evidence

`6991_day0_with_selection` chr9R read `SRR33298395.1615337`

- Y' array: `ID1,ID1,ID1` (Y' Gain; telomere reached: False)
- spacer: switch_detected (conf -0.0007); X element: no_change (conf 0.8756)
- Y' fingerprint: `-` (length 3); path: `chr4R[5-7]:ID1,ID1,ID1`
- donor named: **chr4R** (subtelomere_switch)

| | old | new |
|---|---|---|
| score | `overall_confidence` = **-0.0008** | `recombination_confidence` = **0.7**, `donor_confidence` = **1.0** |
| why | `base=-0.00, complexity=1.10` | `call[spacer=0.00;yprime=0.70] donor[support=1.00;margin=1.00]` |

### Positive control: chr11L -> chr11R

`7172_day6_with_selection` chr11L read `SRR33298426.293534`

- Y' array: `` (No Change; telomere reached: False)
- spacer: no_change (conf 0.95); X element: full_switch (conf 1.0)
- Y' fingerprint: `-` (length 0); path: `-`
- donor named: **chr11R** (subtelomere_switch)

| | old | new |
|---|---|---|
| score | `overall_confidence` = **1.0** | `recombination_confidence` = **1.0**, `donor_confidence` = **0.6** |
| why | `base=1.00, complexity=1.00` | `call[x=1.00] donor[support=0.60;margin=1.00]` |

### Y' Loss on a read that never reaches the telomere (day 0)

`6991_day0_with_selection` chr2L read `SRR33298395.226515`

- Y' array: `` (Y' Loss; telomere reached: False)
- spacer: no_change (conf 0.95); X element: no_change (conf 0.8914)
- Y' fingerprint: `-` (length 0); path: `-`
- donor named: **ambiguous** (array_contraction_unconfirmed)

| | old | new |
|---|---|---|
| score | `overall_confidence` = **0.15** | `recombination_confidence` = **0.4**, `donor_confidence` = **0.0** |
| why | `base=0.30, complexity=1.00` | `call[yprime=0.40] donor[no_donor]` |

### Donor left ambiguous

`7172_survivor_IT155` chr16R read `SRR33298416.419352`

- Y' array: `ID7,ID4,ID2,ID4,ID2,ID4,ID2,ID4,ID2,ID4,ID2,ID2,ID4,ID2,ID4,ID4` (Y' Gain; telomere reached: True)
- spacer: no_change (conf 0.95); X element: no_change (conf 0.7542)
- Y' fingerprint: `-` (length 0); path: `chr2L|chr6L:ID4 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr2L|chr6L:ID4 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr2L|chr6L:ID4 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr2L|chr6L:ID4 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr2L|chr6L:ID4 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr2L|chr6L:ID4 > chr12L|chr13L|chr14L|chr8L|chr8R:ID2 > chr2L|chr6L:ID4 > chr2L|chr6L:ID4`
- donor named: **ambiguous** (unmatched_array)

| | old | new |
|---|---|---|
| score | `overall_confidence` = **0.3** | `recombination_confidence` = **1.0**, `donor_confidence` = **0.0** |
| why | `base=0.30, complexity=1.00` | `call[yprime=1.00] donor[no_donor]` |

### Tandem amplification of the read's own end

`6991_day0_with_selection` chr14L read `SRR33298395.1415479`

- Y' array: `ID1,ID1,ID2,ID2,ID2,ID1,ID2,ID2,ID2` (Y' Gain; telomere reached: True)
- spacer: no_change (conf 0.95); X element: no_change (conf 0.9072)
- Y' fingerprint: `-` (length 4); path: `self[2-5]:ID1,ID2,ID2,ID2`
- donor named: **chr14L** (tandem_amplification_same_end)

| | old | new |
|---|---|---|
| score | `overall_confidence` = **0.33** | `recombination_confidence` = **0.9375**, `donor_confidence` = **1.0** |
| why | `base=0.30, complexity=1.10` | `call[yprime=0.94] donor[support=1.00;margin=1.00]` |

### Unchanged read that stops before the telomere

`6991_day0_with_selection` chr2L read `SRR33298395.682704` -- Y' array `ID4` matches the reference as far as the read goes.

| | old | new |
|---|---|---|
| score | `overall_confidence` = **0.95** | `recombination_confidence` = **0.7** |
| why | `no_recombination` | `call[no_change;end_unconfirmed] donor[n/a]` |

