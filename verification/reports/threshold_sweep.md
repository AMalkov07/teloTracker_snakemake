# Mismatch rate vs grouping threshold, per reference

Built with `verification/threshold_sweep.py`. For each day-0 reference the Y' grouping is cut
at a series of identity thresholds (the pipeline's `--stop-mode threshold` path) and every read
whose Y' copy count matches its anchor end's reference array is scored. A copy is a mismatch
when its matched element's group differs from the group of the element that positionally
belongs there. A read counts once if any of its copies mismatched.

**These are mismatch rates, not error rates.** At day 0 a mismatch is not automatically wrong —
some are genuine recombinants and some are standing variation already in the culture. Treat
them as an upper bound.

## Read-level mismatch %

| reference | reads | 95 | 96 | 97 | 98 | 99 | 99.5 | silhouette | curated variant | curated family |
|---|---|---|---|---|---|---|---|---|---|---|
| 6991_day0 | 3,489 | 0.115 | 0.115 | 0.115 | 0.115 | 0.115 | 0.115 | 0.115 | **1.204** | 0.143 |
| 6991_day0_TeloTag | 7,428 | 0.054 | 0.054 | 0.067 | 0.067 | 0.081 | 0.350 | 0.054 | **0.539** | 0.094 |
| 6991_day0_TeloTag_with_selection | 8,765 | 0.125 | 0.125 | 0.125 | 0.125 | 0.125 | 0.137 | 0.125 | 0.194 | 0.148 |
| 6991_day0_reference | 7,511 | 0.359 | 0.359 | 0.359 | 0.359 | 0.359 | 0.373 | 0.346 | **0.533** | 0.320 |
| 6991_day0_reference_promethion | 7,180 | 0.153 | 0.153 | 0.153 | 0.153 | 0.167 | 0.167 | 0.139 | 0.251 | 0.139 |
| 6991_day0_with_selection | 27,409 | 1.594 | 1.598 | 1.598 | 1.598 | 1.602 | 1.605 | *0.179* | 1.649 | *0.182* |
| 6991_day0_with_selection_repeat | 12,420 | 0.145 | 0.145 | 0.145 | 0.145 | 0.145 | 0.145 | 0.137 | 0.209 | 0.121 |
| 6991_day0_with_selection_repeat2 | 5,391 | 0.241 | 0.241 | 0.241 | 0.241 | 0.241 | 0.241 | 0.204 | 0.241 | 0.185 |
| 7172_day0_with_selection | 4,471 | 0.089 | 0.089 | 0.089 | 0.089 | 0.089 | 0.089 | 0.089 | 0.089 | 0.089 |
| 7302_day0_with_selection | 7,112 | 0.197 | 0.197 | 0.197 | 0.197 | 0.197 | 0.211 | 0.197 | 0.267 | 0.197 |

## Number of groups

| reference | elements | 95 | 96 | 97 | 98 | 99 | 99.5 | silhouette | curated variant | curated family |
|---|---|---|---|---|---|---|---|---|---|---|
| 6991_day0 | 34 | 9 | 10 | 10 | 11 | 12 | 14 | 8 | 10 | 6 |
| 6991_day0_TeloTag | 27 | 9 | 10 | 11 | 11 | 13 | 15 | 8 | 10 | 6 |
| 6991_day0_TeloTag_with_selection | 34 | 9 | 10 | 10 | 11 | 12 | 14 | 8 | 10 | 6 |
| 6991_day0_reference | 33 | 9 | 10 | 10 | 11 | 12 | 14 | 8 | 10 | 6 |
| 6991_day0_reference_promethion | 34 | 9 | 10 | 10 | 11 | 12 | 14 | 8 | 10 | 6 |
| 6991_day0_with_selection | 34 | 10 | 11 | 11 | 12 | 13 | 15 | 8 | 10 | 6 |
| 6991_day0_with_selection_repeat | 34 | 9 | 10 | 10 | 11 | 12 | 14 | 8 | 10 | 6 |
| 6991_day0_with_selection_repeat2 | 34 | 9 | 10 | 10 | 11 | 12 | 14 | 8 | 10 | 6 |
| 7172_day0_with_selection | 35 | 9 | 10 | 10 | 11 | 12 | 14 | 8 | 9 | 6 |
| 7302_day0_with_selection | 36 | 9 | 10 | 10 | 11 | 12 | 14 | 8 | 12 | 7 |

## Curated variant is beaten by a threshold cut at the same resolution

The fair comparison is against the threshold cut producing the **same number of groups**, since
resolution and mismatch rate trade off:

| reference | curated variant | groups | equivalent cut | that cut's rate |
|---|---|---|---|---|
| 6991_day0 | **1.204 %** | 10 | 96 % | **0.115 %** |
| 6991_day0_TeloTag | **0.539 %** | 10 | 96 % | **0.054 %** |
| 6991_day0_reference | 0.533 % | 10 | 96 % | 0.359 % |
| 6991_day0_reference_promethion | 0.251 % | 10 | 96 % | 0.153 % |
| 6991_day0_with_selection_repeat | 0.209 % | 10 | 96 % | 0.145 % |
| 6991_day0_TeloTag_with_selection | 0.194 % | 10 | 96 % | 0.125 % |
| 6991_day0_with_selection_repeat2 | 0.241 % | 10 | 96 % | 0.241 % |
| 7302_day0_with_selection | 0.267 % | 12 | 99 % | 0.197 % |
| 7172_day0_with_selection | 0.089 % | 9 | 95 % | 0.089 % |

**Curated variant is worse than, or equal to, the equivalent threshold cut on all ten
references** — never better. On `6991_day0` it is 10x worse (1.204 % vs 0.115 %) for the same
10 groups. The curated scheme buys no resolution that a threshold cut does not also buy, and
costs accuracy to do it.

The likely cause is the one already documented: the curated scheme separates `ID5_Blue-Light`
(chr14L-1) from `ID5_Blue-Dark` (chr7R-1, chr16L-1) on a distinction resting on 1-2 bp across
6.6 kb — below ONT read accuracy. A sequence-derived cut never makes that split because the
sequences do not support it.

## Curated family is competitive, and matches silhouette

| | curated family (6-7 groups) | silhouette (8 groups) |
|---|---|---|
| better on | 4 references | 4 references |
| tied | 2 references | |

Ranges 0.089-0.320 % against silhouette's 0.089-0.346 %, excluding the defective reference.
At this coarse resolution the two are interchangeable; curated family gets there with one or
two fewer groups.

## The headline: the rate is now flat, and that is the point


Eight of the ten references vary by **less than 0.05 percentage points across the whole 95-99
range**, while the group count rises from 9 to 12. Going from 95 % to 99 % buys **three extra
groups for essentially nothing**; against the silhouette default it is **four extra groups
(8 → 12)** at the same mismatch rate.

This is a change from the earlier benchmark, where cut99 scored 6.92 % against cut97's 0.17 %.
That gap was entirely the chr16L boundary defect: the over-extended element split into its own
group at the finer cut and then mis-sorted every chr7R and chr14L read. With the boundary
corrected, chr16L merges with chr7R and chr14L at every threshold and the sensitivity is gone.
**The threshold was never the problem; the reference was.**

## The two references that still move

**`6991_day0_with_selection` — 1.59-1.61 % at every threshold but 0.179 % under silhouette.**
Its `chr14L-1` is mis-assembled at **5,720 bp** against 6,654 in every clean reference. Being
900 bp short it splits into its own group at any threshold cut, and mis-sorts every chr14L read.
Silhouette's coarser 8-group partition absorbs it and hides the defect — the low silhouette
number here is concealment, not accuracy.

**`6991_day0_TeloTag` — 0.054 % up to 99, then 0.350 % at 99.5.** The only reference where the
finest cut genuinely costs something, and the only one with 27 elements rather than 33-36 (the
known missing chr4R array). 99.5 % is the point where this reference starts splitting groups it
should not.

## Practical reading

99 % looks like the right operating point: four more groups than silhouette at no measurable
cost on nine of ten references. 99.5 % is where the first reference starts to break, so it is
the edge of what the data supports. And a flat curve is a reference-quality check in itself —
if one reference shows threshold sensitivity while its siblings do not, suspect the reference
before the threshold.
