# Recombination confidence audit

Snapshot: `verification/snapshot_v2c_path` -- 28 samples, 234,899 reads. Read-only: this measures the current `overall_confidence`; nothing is changed.

## 1. What values the score takes, by kind of read

Of 33,987 recombinant reads, **59.3% sit exactly on the 0.30 / 0.33 floor** (base 0.3, x1.1 when the Y' array agrees). Every non-recombinant read is a constant.

| event_class | reads | mean | median | distinct values | % at 0.95 | % at 0.30/0.33 |
|---|---|---|---|---|---|---|
| no recombination | 200912 | 0.950 | 0.950 | 2 | 99.892 | 0.000 |
| Y' | 24147 | 0.344 | 0.330 | 5 | 0.000 | 82.963 |
| X+Y' | 4511 | 0.789 | 0.809 | 1798 | 0.842 | 0.399 |
| X | 3047 | 0.822 | 0.914 | 1281 | 18.149 | 0.000 |
| X+Y' (complex) | 1542 | 0.539 | 0.520 | 863 | 0.000 | 0.000 |
| spacer+Y' | 259 | 0.188 | 0.058 | 123 | 0.000 | 42.085 |
| spacer+X+Y' | 159 | 0.838 | 0.900 | 38 | 0.629 | 0.629 |
| spacer+X+Y' (complex) | 159 | 0.593 | 0.641 | 70 | 0.000 | 1.258 |
| spacer+X | 88 | 0.903 | 0.900 | 27 | 1.136 | 0.000 |
| spacer | 51 | 0.719 | 0.900 | 18 | 0.000 | 0.000 |
| spacer+X (complex) | 17 | 0.617 | 0.639 | 9 | 0.000 | 0.000 |
| spacer+Y' (complex) | 7 | 0.270 | 0.000 | 5 | 0.000 | 0.000 |

Most common `confidence_factors` among recombinant reads:

- `base=0.30, complexity=1.00`: 10,678
- `base=0.30, complexity=1.10`: 10,416
- `base=0.60, complexity=1.10`: 1,944
- `base=0.40, complexity=1.10`: 1,353
- `base=0.95, complexity=1.00`: 954
- `base=0.74, complexity=1.10`: 449
- `base=0.73, complexity=1.10`: 350
- `base=0.90, complexity=1.00`: 251

## 2. Is the per-end mean just the recombination rate restated?

Across 883 (sample, end) pairs, fitting `mean_confidence = a*(1-p) + b*p` (p = fraction recombinant) gives a = 0.943, b = 0.450, **R^2 = 0.824** (Pearson r = -0.908). a ~ 0.95 is the constant given to every non-recombinant read. An R^2 this high means the per-end mean carries almost no information beyond the recombination rate already reported next to it.

## 3. Does confidence separate correct from wrong donor calls?

AUC = probability a correctly attributed read scores higher than a wrongly attributed one. 0.5 = no information, 1.0 = perfect.

| truth set | correct donor | wrong donor | mean conf (correct) | mean conf (wrong) | AUC |
|---|---|---|---|---|---|
| positive control (7172 chr11L -> chr11R) | 516 | 23 | 0.922 | 0.773 | 0.829 |
| chr13L truth set, strict (7302) | 36 | 3 | 0.642 | 0.559 | 0.630 |
| both pooled | 552 | 26 | 0.904 | 0.748 | 0.804 |

Wrong calls went to: `chr10L` 6, `chr14R` 5, `chr4L` 3, `chr14L` 3, `chr4R` 2, `chr1R` 2

## 4. Does confidence separate day-0 calls from real events?

At day 0 there should be no recombination, so day-0 recombinant calls are at best standing variation and at worst false. If the score measured "did recombination happen", they should score clearly lower than calls at later timepoints.

| recombinant reads in | n | mean conf | median conf |
|---|---|---|---|
| day-0 self-runs | 1602 | 0.362 | 0.330 |
| later timepoints | 9445 | 0.517 | 0.330 |
| survivors | 22940 | 0.442 | 0.330 |

**AUC (later timepoint vs day-0) = 0.633.**

