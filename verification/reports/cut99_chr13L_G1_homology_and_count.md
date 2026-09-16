# Why G1 specifically? Homology doesn't explain it -- but the "15 events" aren't 15 events

## Part 1: local homology tracts don't favour G1 either

The earlier ranking used overall identity. HR strand invasion cares about a contiguous
near-perfect tract, not a genome-wide average, so the same candidates were re-checked for
their single best unbroken high-identity block against chr13L-1:

| candidate | best block identity | length | >=99% block | >=99.5% block |
|---|---|---|---|---|
| chr14L-3/4/5 (chr13L-1's own group -- invisible to this method) | 99.60% | 5,482 bp (full) | 5,482 bp | 5,482 bp |
| chr12L-1 | 99.17% | 5,188 bp | 5,188 bp | 0 |
| chr8R-1 | 98.56% | 5,495 bp | 0 | 0 |
| chr8L-1 | 98.14% | 5,059 bp | 0 | 0 |
| chr16R-1 | 98.61% | 3,588 bp | 0 | 0 |
| **chr6L-1 (G1)** | **97.27%** | 3,884 bp | **0** | **0** |
| **chr2L-1 (G1)** | **97.12%** | 3,883 bp | **0** | **0** |

Not one single alignment block between chr13L-1 and either G1 member ever reaches 99%
identity anywhere along its length. chr12L-1 has a 5,188 bp block at >=99%; G1 has nothing
at that stringency at all. **The local-tract version of the homology hypothesis fails even
more clearly than the overall-identity version did.**

Also checked the X element (the other candidate homology source, since Y' elements are not
the only shared sequence at a chromosome end): chr13L's X element is 88-94% identical to
every candidate's X element with no differential signal -- chr16R is nominally highest (94.2%,
but only 34% coverage) and G1's members are unremarkable (88.9%, 90.4%). The X element does
not explain the preference either.

**Sequence homology, at both the Y' and the X element, rules out "G1 is simply the best
available match" as the explanation.**

## Part 2: the reframe -- 15 events is not 15 independent occurrences

Both signature classes were tested for cross-read identity, the same test that established
these are standing variation rather than chimeras:

**The 9 whole-element reads** (already reported): pairwise 94-99% identity, full length,
across 3 independent SRA runs -- one shared haplotype.

**The 5 partial-junction reads (donor chr2L-1), now checked the same way:**

| pair | identity | length | independent samples? |
|---|---|---|---|
| SRR33298373.8939 vs SRR33298373.1520960 | 98.98% | 6,051 bp | same sample |
| SRR33298377.125819 vs SRR33298373.8939 | 96.50% | 6,117 bp | **different** (reference_promethion vs with_selection_repeat) |
| SRR33298377.267122 vs SRR33298373.1520960 | 96.00% | 6,146 bp | **different** |
| SRR33298461.44252 vs SRR33298373.1520960 | 95.53% | 4,251 bp | **different** (6991_day0) |

Full-length, high-identity matches across the same three independent SRA runs
(SRR33298373, SRR33298377, SRR33298461). **This is also one shared haplotype**, not five
independent partial-recombination events.

## The resolution

All 15 chr13L -> G1 "events" reduce to **two distinct ancestral molecules**:

1. one whole-Y'-replacement haplotype (9 reads, 3 independent preps)
2. one partial-junction haplotype (5-6 reads, 3 independent preps)

Both are resampled repeatedly because each independently prepared library draws reads from
the same underlying population of cells, and both haplotypes are apparently present at
appreciable frequency in that population. The count of "15" was never 15 independent
recombination events choosing G1 over ten other options -- it is closer to **2**.

This resolves the original puzzle without requiring G1 to be an exceptional homology match.
A single historical recombination event does not need to have used the best available
partner at the time it happened; it only had to happen once, successfully, and then persist.
Two such events landing on the same donor pair is a much weaker coincidence to explain than
fifteen would be, though it is not zero -- something (possibly spatial/nuclear proximity of
these three subtelomeres, unmeasured here, or simply chance among a small number of historical
events) still put both of chr13L's known exchange partners in the same place.

## What this means for the numbers reported earlier

The per-dataset and combined tables that reported "18 events, 15 to G1" should be read as
**"18 reads sampled, reducible to a small number of distinct underlying recombination
haplotypes, at least 2 of which are chr13L<->G1"** -- the tables were correct as read-counts,
but read-count is not event-count when preparations are resampling the same population.
