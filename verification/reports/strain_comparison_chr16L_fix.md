# chr16L boundary fix across 6991, 7172, 7302

## 1. The defect, measured independently on each reference

`fix_yprime_boundary.py` BLASTs `chr16L_Y_Prime_1` against every element it is ≥99 %
identical to and ≥90 % of the longer length, and takes the median implied trim across
donors that agree within 5 bp. Applied with no strain-specific input:

| reference | old length | trim | new length | donors agreeing |
|---|---|---|---|---|
| 6991_day0 | 6732 | 77 bp | 6655 | 17 |
| 6991_day0_TeloTag | 6731 | 77 bp | 6654 | 10 |
| 6991_day0_TeloTag_with_selection | 6732 | 77 bp | 6655 | 17 |
| 6991_day0_reference | 6731 | 77 bp | 6654 | 16 |
| 6991_day0_reference_promethion | 6732 | 77 bp | 6655 | 17 |
| 6991_day0_with_selection | 6731 | 77 bp | 6654 | 16 |
| 6991_day0_with_selection_repeat | 6732 | 77 bp | 6655 | 17 |
| 6991_day0_with_selection_repeat2 | 6732 | 77 bp | 6655 | 17 |
| **7172_day0_with_selection** | 6732 | **77 bp** | 6655 | 19 |
| **7302_day0_with_selection** | 6732 | **77 bp** | 6655 | 18 |

All ten references independently converge on the same 77 bp, with 10–19 elements agreeing
on it each time. That consistency is itself the evidence this is a property of the
labelling step (or of the curated library it inherits boundaries from), not per-assembly
noise: ten unrelated BLAST measurements do not agree to the base pair by chance.

## 2. False recombination-call rate on day-0 reads, before vs after

Reads that should show no recombination (anchored copy count matches the reference array)
whose Y' array nonetheless disagrees with their own anchor's reference array.

| scheme | 6991 (8 refs, 79,590 reads) | 7172 (4,471 reads) | 7302 (7,112 reads) |
|---|---|---|---|
| **element** (raw match, no grouping) | 8.95 → **4.00** % | 13.31 → **8.16** % | 14.75 → **8.93** % |
| **condensed** | 7.21 → **0.75** % | 8.21 → **0.20** % | 7.75 → **0.28** % |
| **cut99** | 7.12 → **0.67** % | 8.10 → **0.09** % | 7.68 → **0.20** % |
| cut97 | 0.66 → 0.66 % | 0.09 → 0.09 % | 0.20 → 0.20 % |
| **curated_variant** | 1.44 → **0.81** % | 2.97 → **0.09** % | 1.88 → **0.27** % |
| silhouette | 0.16 → 0.17 % | 0.09 → 0.09 % | 0.20 → 0.20 % |
| curated_family | 0.27 → 0.17 % | 0.13 → 0.09 % | 0.28 → 0.20 % |

## 3. Group count, before → after

| scheme | 6991 | 7172 | 7302 |
|---|---|---|---|
| condensed | 18 → 17 | 19 → 18 | 18 → 17 |
| cut99 | 13 → 12 | 13 → 12 | 13 → 12 |
| cut97 | 10 → 10 | 10 → 10 | 10 → 10 |
| curated_variant | 10 → 10 | 9 → 9 | 12 → 12 |
| silhouette | 8 → 8 | 8 → 8 | 8 → 8 |
| curated_family | 6 → 6 | 6 → 6 | 7 → 7 |

Every strain loses exactly one `condensed`/`cut99` group (chr7R/chr14L/chr16L merging into
one), and no scheme's group count changes anywhere else — the fix is surgical.

## 4. What's consistent across all three, and what isn't

**Consistent everywhere:**
- The trim itself: 77 bp, every reference.
- `cut97` and `silhouette` are untouched (they already grouped the three together).
- `condensed` and `cut99` both fall by roughly an order of magnitude, landing in the
  same 0.09–0.75 % band as the untouched schemes.
- `element`-level (raw matching, no grouping applied) improves substantially but stays
  the highest of any scheme on all three strains — it is the floor of what boundary
  correction alone can fix; the residual error there is a different, unrelated set of
  confusions (see `verification/reports/grouping_benchmark/cut97_errors_*` for the kind
  of thing that remains).

**Differs by strain — `curated_variant`:**
Both 7172 and 7302 recover to the same floor as every other scheme (0.09 % / 0.27 %).
6991 only partially recovers (0.81 %), because its curated library's `ID5_Blue-Dark`
(chr7R1/chr16L1) vs `ID5_Blue-Dark`/`ID5_Blue-Light` (chr14L1) split — unrelated to this
boundary and unfixed by it — happens to catch far more anchored reads in 6991's
population than in 7172 or 7302's. All three curated libraries carry the identical
label split; only how often it's hit differs.

Source data: `verification/grouping_benchmark/fixed_beds/trim_report.tsv`,
`verification/reports/grouping_benchmark_{uncorrected8,fixed}/grouping_benchmark.md` (6991),
`verification/reports/grouping_benchmark_7172_7302_{uncorrected,fixed}/grouping_benchmark.md`
(7172/7302).
