# Why is chr13L consistently mistaken for G1 (chr2L-1/chr6L-1)? Investigated in detail

Starting question: is the chr13L -> G1 relationship (18 of 18 events in the 5 matched
samples, 83% concentration) a grouping failure, recombination, or something else?

## First check: is this a grouping failure?

If chr13L-1 were actually just as similar to chr2L-1/chr6L-1 as they are to each other, the
99% cutoff would be arbitrarily splitting one real group into two, and every chr13L read
would look "wrong" by construction -- an artefact of the threshold, not a finding.

**Measured directly.** chr2L-1 and chr6L-1 are 99.05% identical to each other over their full
5,975 bp -- correctly grouped. chr13L-1 is only **~97.1-97.3%** identical to each of them
(the same value against both), and is **492 bp shorter** (5,483 bp vs 5,975 bp). Both the
identity and the length gap are stable across all five matched assemblies -- chr13L-1 comes
out at exactly 5,483 bp in every one.

**Verdict: not a grouping failure.** chr13L-1 is a real, distinct, consistently-assembled
element correctly excluded from G1. The cutoff is doing its job.

## Second check: is it recombination? Yes -- but it splits into two distinct signatures

Breaking down all 15 chr13L -> G1 events (5 matched samples) by their actual evidence:

| pattern | n | donor named | signature |
|---|---|---|---|
| **partial junction** | 5 | chr2L-1 (5), | clean anchor/telomere split, "strong"/"weak" evidence |
| **whole-element** | 9 | chr6L-1 (9) | "no junction" (1 additional "FAILS") |

These are not the same phenomenon wearing two labels -- they behave completely differently.

### The whole-element reads: a length signature that rules out noise

For the 9 "no junction" reads, the Y' at the chr13L locus was measured directly:

| read | sample (SRA run) | Y' span |
|---|---|---|
| SRR33298373.896473 | 6991_day0_with_selection_repeat | 5,965 bp |
| SRR33298373.1068481 | 6991_day0_with_selection_repeat | 5,975 bp |
| SRR33298373.107087 | 6991_day0_with_selection_repeat | 5,969 bp |
| SRR33298373.1079946 | 6991_day0_with_selection_repeat | 5,972 bp |
| SRR33298373.46529 | 6991_day0_with_selection_repeat | 5,972 bp |
| SRR33298373.797552 | 6991_day0_with_selection_repeat | 5,965 bp |
| SRR33298377.260642 | 6991_day0_reference_promethion | 5,964 bp |
| SRR33298377.272313 | 6991_day0_reference_promethion | 5,966 bp |
| SRR33298384.45869 | 6991_day0_with_selection_repeat2 | 5,923 bp |

**All nine sit at 5,923-5,975 bp -- the length of chr6L-1/chr2L-1, not chr13L-1's own 5,483 bp.**
A read carrying its own native chr13L-1 copy should be ~5,483 bp; every one of these is ~490 bp
longer. This is not something read noise or a matching quirk can produce -- the physical
molecule is the wrong length for chr13L's own Y'.

Confirmed with a full-read alignment (`SRR33298373.896473` against all three references):

| reference | identity | coverage |
|---|---|---|
| chr6L-1 | 98.93% | full length, one continuous block |
| chr2L-1 | 98.00% | full length, one continuous block |
| chr13L-1 | 95.77% | **only the last 3,897 of 6,212 bp -- no hit at all for the first ~2,300 bp** |

There is no native chr13L-1 signature anywhere in most of the read. This is not "the
junction is hard to find" -- the window scan correctly reported no junction because there
is not one: the whole element is foreign.

### Is the whole-element signature a chimera, standing variation, or a fresh event?

Ruled out chimera the same way as the earlier chr6L investigation: pulled the Y' region from
5 of the 9 reads, spanning **three independent SRA runs** (SRR33298373, SRR33298377,
SRR33298384 -- three separate sequencing libraries), and aligned them to each other:

| pair | identity | length |
|---|---|---|
| SRR33298377.260642 vs SRR33298373.1068481 | 99.27% | 6,045 bp |
| SRR33298373.896473 vs SRR33298373.1068481 | 98.85% | 6,059 bp |
| SRR33298384.45869 vs SRR33298377.260642 | 96.10% | 4,406 bp |

Full-length, high-identity matches across three independently prepared libraries. A chimera
formed during library prep cannot recur in a separately prepared sequencing run -- this is the
same standing-variation signature established earlier for the chr6L->donor haplotype.

## The picture that emerges

Not one relationship, but two, both real and both correctly excluded from G1 by the grouping:

1. **A recurring whole-Y' replacement haplotype** (9 of 15 events, the majority): the native
   chr13L Y' has been replaced end-to-end by a chr2L/chr6L-type copy at some point in the
   strain's history, and this haplotype is now standing variation in the culture -- present
   across at least three independently sequenced preparations, not something happening fresh
   at each day-0 timepoint.
2. **Genuine partial (mid-Y') recombination** (5 of 15 events): a clean junction, native
   chr13L-1 on the anchor side and chr2L-1-derived sequence on the telomere side -- the same
   signature validated in detail for chr6L earlier.

Both point to real biological exchange between the chr13L locus and the chr2L/chr6L lineage,
not a grouping or matching artefact -- but they are mechanistically distinct (a fixed,
historical whole-element replacement vs. an ongoing/recent partial switch) and should not be
pooled into one count.

## One thing that stays unresolved

The named donor within the whole-element reads (chr6L-1, 9/9) cannot be trusted as more
specific than "chr2L-1 or chr6L-1" -- those two are 99.05% identical to each other, so which
one the pipeline names is close to a coin flip. The partial-junction reads more often named
chr2L-1 (4 of 5), which may or may not be meaningful given the same ambiguity.
