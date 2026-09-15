# The chr6L "mislabelled" reads are recombinant Y' elements, not matching failures

## Correction to an earlier conclusion

An earlier pass over these reads (`verification/reports/grouping_benchmark_7172_7302_fixed/cut99_errors_analysis.md`)
concluded they were **not** recombination, on the grounds that none of the confusions
involved a near-identical element -- every "wrong" match sat at 28-72% similarity to the
true element, where a genuine template switch should produce a high-identity donor match.

**That reasoning was backwards.** A read carrying a recombinant Y' -- anchor-proximal half
from one element, telomere-distal half from another -- will by construction match *no single
library entry* end to end. Scoring it against the library one element at a time is guaranteed
to return a mediocre best hit. The very signature I treated as evidence against recombination
is what recombination looks like under that test.

## The correct test, and what it shows

Splice the two candidate references at the homology block and ask whether the read matches
the splice product in one clean full-length alignment (`verification/test_recombinant_hybrid.py`).

For `SRR33298452.400002` (7302 day-0, anchored chr6L, 23,890 bp):

| reference | result |
|---|---|
| chr6L1 alone (5,975 bp, its positional truth) | fragments into 3,860 bp @ 98.6% + 2,040 bp @ 96.7% |
| chr7R1 alone (6,655 bp) | fragments similarly |
| **chr6L1[1-3433] + chr7R1[4619-6655]** | **ONE block, 5,473 bp @ 99.762%** |

Across all 8 chr6L misreads in 7302/7172, with 12 correctly-assigned chr6L reads as controls:

| set | n | verdict |
|---|---|---|
| misreads | 5 of 8 | **RECOMBINANT** (hybrid bitscore > 1.15x either single reference) |
| misreads | 3 of 8 | single-element (different junction or donor; not resolved here) |
| controls | 12 of 12 | single-element -- each matches chr6L1 in one full-length block, 94-99.6% |

The controls are the load-bearing part: they confirm the chr6L1 reference is sound. Its own
reads match it end to end. So the hybrid reads are genuinely structurally different, not an
artifact of a bad reference.

## Structure of the recombinant

Read orientation telomere -> anchor:

```
[telomere][~1.6 kb Long-Y'-type @ 99.94%][425 bp shared cassette][~3.4 kb chr6L1-type @ 99.65%][chr6L anchor]
                                          ^ junction lies in here
```

The junction is bounded to roughly chr6L1 positions 3433-3858 -- inside a 425 bp cassette that
is 99.8-100% identical across 25 of 36 elements in the library. A junction falling in the
region of maximal homology is what the mechanism predicts, and it is also why the junction
cannot be localised more precisely than that window.

**The donor cannot be pinpointed.** chr7R1, chr16L1, chr5R1, chr13L2, chr13L4 and chr14L1 all
match the telomere-distal flank at 99.94-100%. chr7R1 is used above as a representative, not
an identification.

Five of the six hybrid reads share the same junction bounds (chr6L1 block ending ~3858,
Long-Y' block starting ~2999), consistent either with one clonal event that expanded or with
a hotspot at the homology-block boundary. The data here do not separate those.

## Frequency

8 of 595 chr6L reads in 7302 day-0 (1.3%); 568 (95.5%) assign cleanly to chr6L1. A low-level
recombination background at day 0 of this order is what was originally hypothesised.

## Caveat not excluded

These are day-0 "with selection" samples, and a chimera formed during library prep would also
tend to join at a homologous region. Nothing here distinguishes an in vivo recombinant from a
library chimera; the ONT prep used (ligation vs amplification) and whether the same junction
recurs in independent preps would be the way to settle it.
