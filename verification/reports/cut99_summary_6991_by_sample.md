# 6991 day-0, split by reference

Each of the eight 6991 day-0 populations is assembled independently and gets its own Y'
element library and its own cut99 grouping, so they are eight separate analyses rather than
eight replicates of one. Pooling them hides both the reference defect and the spread in rate.

| reference (sample) | scored reads | elements | groups | chr14L-1 len | mismatches | strong | weak | FAILS | no junction | ref defect | supported | rate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 6991_day0_with_selection | 27,409 | 34 | **13** | **5,720** | 439 | 20 | 5 | 12 | 14 | **388** | 25 | 0.091 % |
| 6991_day0_with_selection_repeat | 12,420 | 34 | 12 | 6,654 | 18 | 5 | 2 | 2 | 9 | 0 | 7 | 0.056 % |
| 6991_day0_TeloTag_with_selection | 8,765 | 34 | 12 | 6,654 | 11 | 5 | 0 | 4 | 2 | 0 | 5 | 0.057 % |
| 6991_day0_reference | 7,511 | 33 | 12 | 6,654 | 27 | 10 | 6 | 8 | 3 | 0 | 16 | **0.213 %** |
| 6991_day0_TeloTag | 7,428 | **27** | **13** | 6,653 | 6 | 2 | 1 | 0 | 3 | 0 | 3 | 0.040 % |
| 6991_day0_reference_promethion | 7,180 | 34 | 12 | 6,654 | 12 | 3 | 3 | 2 | 4 | 0 | 6 | 0.084 % |
| 6991_day0_with_selection_repeat2 | 5,391 | 34 | 12 | 6,654 | 13 | 5 | 4 | 1 | 3 | 0 | 9 | **0.167 %** |
| 6991_day0 | 3,489 | 34 | 12 | 6,653 | 4 | 2 | 1 | 1 | 0 | 0 | 3 | 0.086 % |
| **TOTAL** | **79,593** | | | | **530** | **52** | **22** | **30** | **38** | **388** | **74** | **0.093 %** |

"supported" = strong + weak, i.e. reads where a mid-Y' junction is evidenced by both halves.

## The defective reference distorts its own grouping

`6991_day0_with_selection` carries a mis-assembled `chr14L-1` (5,720 bp against 6,654 in every
other sample). The consequence is not only that its reads mismatch -- the defect changes the
grouping itself:

| sample | chr14L-1 sits in | with |
|---|---|---|
| 6991_day0_with_selection | G10 | **nothing — a singleton group** |
| 6991_day0_with_selection_repeat | G8 | chr7R-1, chr16L-1, chr4R-1..7, chr12R-2..7, chr14L-2, chr15R-1 (17 others) |

Because the truncated element is 900 bp short, it no longer clusters with the 17 elements it
belongs with, and splits off alone. That phantom singleton is what produces 385 of the 388
"reference defect" calls, and it is why this sample reports 13 groups where its siblings
report 12.

This is the same failure mode as the chr16L boundary defect found earlier: a mis-bounded
element becomes its own group, and that group then mis-sorts every read at the affected end.
`6991_day0_TeloTag` also reports 13 groups, but for an unrelated reason — it has only 27
elements (the known missing chr4R array), so its clustering differs structurally.

## The rate is not stable across preparations

Excluding the reference defect, supported events range from **0.040 % to 0.213 %** — a
five-fold spread across eight preparations of the same strain. The two extremes
(`6991_day0_TeloTag` at 0.040 %, `6991_day0_reference` at 0.213 %) differ by more than
Poisson noise on these counts would comfortably explain, so the pooled 0.093 % should be
treated as an average over heterogeneous preparations rather than a property of the strain.

Note also that `6991_day0_with_selection` still yields 25 supported events at 0.091 % once its
defect rows are set aside, so the defect inflates its mismatch count without much changing its
underlying recombination rate.
