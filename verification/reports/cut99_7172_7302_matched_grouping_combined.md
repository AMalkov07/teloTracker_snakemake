# 7172 and 7302: cut99 Y' grouping mismatches, same method as the 6991 analysis

Unlike 6991 (5 identically-grouped day0 variant sequencing preps), 7172 and 7302 each have
**only one** day0 self-check run (`7172_day0_with_selection`, `7302_day0_with_selection`) --
there is no cross-sample grouping-identity step to do first, and far fewer mismatches to work
with: **4 for 7172, 14 for 7302** (vs. 58 pooled across 5 samples for 6991). Both strains
cluster into the same 12 groups as 6991, with the same group-to-chr_end architecture (G1 =
chr2L/chr6L, G2 = chr13L/chr14L family, G3 = chr8R solo, G8 = the big tandem-array group,
etc.) -- the underlying element library is the same species, just different sequenced isolates.

Method, identical to the corrected 6991 analysis: for every mismatched read, compute a true
global alignment (`edlib`, infix `HW` mode, both strands checked) of the read's Y' region
against its own expected reference and against the best member of the donor group. Gap =
donor% - own%. Thresholds established on 6991's native controls (consistently -4.8% to -10.8%
for correctly-assigned reads) carry over directly: gap > +5% = whole-element swap; -2% to +5% =
"mid-Y' partial junction, unconfirmed breakpoint" (real signal, not yet shown to be a discrete
crossover); gap < -2% = not supported. **No fresh 7172/7302-specific control reads were pulled**
-- the per-chr-end recombination feature files that the mismatch tables were built from are no
longer present on Argon in an easily-relocatable form, so the already-established 6991 baseline
range is reused rather than re-derived. This is a reasonable assumption (same library, same
species) but is a real limitation worth flagging.

## Group key

Both strains: 12 groups, same composition pattern as 6991 (sizes differ slightly because the
tandem arrays at chr12R/chr4R/chr13L/chr14L/chr16L have different copy counts per strain).

**7172** (35 elements):

| group | size | members |
|---|---|---|
| G8 | 21 | chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr12R-8, chr14L-1, chr14L-2, chr15R-1, chr16L-1, chr16L-2, chr16L-3, chr4R-1, chr4R-2, chr4R-3, chr4R-4, chr4R-5, chr4R-6, chr4R-7, chr7R-1 |
| G7 | 2 | chr10L-1, chr9L-1 |
| G2 | 2 | chr13L-1, chr14L-3 |
| G1 | 2 | chr2L-1, chr6L-1 |
| G4, G11, G10, G6, G12, G9, G5, G3 | 1 each | chr12L-1, chr12R-1, chr14R-1, chr16R-1, chr5L-1, chr5R-1, chr8L-1, chr8R-1 |

**7302** (36 elements):

| group | size | members |
|---|---|---|
| G8 | 19 | chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr13L-2, chr13L-4, chr14L-1, chr14L-2, chr15R-1, chr16L-1, chr4R-1, chr4R-2, chr4R-3, chr4R-4, chr4R-5, chr4R-6, chr4R-7, chr7R-1 |
| G2 | 5 | chr13L-1, chr13L-3, chr14L-3, chr14L-4, chr14L-5 |
| G7 | 2 | chr10L-1, chr9L-1 |
| G1 | 2 | chr2L-1, chr6L-1 |
| G4, G11, G10, G6, G12, G9, G5, G3 | 1 each | chr12L-1, chr12R-1, chr14R-1, chr16R-1, chr5L-1, chr5R-1, chr8L-1, chr8R-1 |

7302's chr13L has a 4-copy array (2 in G2, 2 in G8), unlike 6991/7172 where chr13L is a single
solo copy in G2 -- this is a genuine array-length difference between strains, not a grouping
error.

## Recipient/donor breakdown

**7172 (4 mismatches, 1 sample):**

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr16R | 2 | G2 (chr14L-3) | 2 | 100% |
| chr5R | 1 | G2 (chr13L-1) | 1 | 100% |
| chr6L | 1 | G2 (chr14L-3) | 1 | 100% |

**All 4 of 4 mismatches in 7172 go to G2.** Tiny N, but a striking signal in the same direction
as 7302 below.

**7302 (14 mismatches, 1 sample):**

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr6L | 7 | G2 (chr14L-5) | 4 | 57% |
| | | G8 (chr16L-1) | 1 | 14% |
| | | G10 (chr14R-1) | 1 | 14% |
| | | G3 (chr8R-1) | 1 | 14% |
| chr14R | 2 | G8 (chr14L-1) | 2 | 100% |
| chr5R | 2 | G7 (chr10L-1) | 1 | 50% |
| | | G10 (chr14R-1) | 1 | 50% |
| chr10L | 1 | G10 (chr14R-1) | 1 | 100% |
| chr2L | 1 | G3 (chr8R-1) | 1 | 100% |
| chr8R | 1 | G1 (chr2L-1) | 1 | 100% |

**chr6L accounts for half of all 7302 mismatches (7/14), and its dominant donor is again G2**
(4 of 7) -- the same group that dominates every 7172 mismatch and that dominated chr13L's
mismatches in 6991. chr2L -> G3 (chr8R-1) reproduces exactly the pattern already established
for 6991. chr8R -> G1 (chr2L-1) is the mirror image of that -- chr8R and chr2L/chr6L mistake
each other in both directions, consistent with genuine mutual sequence homology rather than a
one-way artefact.

## Is chr6L -> G2 homology-driven, like chr13L -> G1 was NOT?

The same homology-ranking test used for the four 6991 recipients (BLAST the recipient's own
reference against every other library element, excluding same-group members which are
structurally invisible to mismatch detection): chr6L-1 (7302) vs. the rest of the library:

| rank | element | group | %identity | bitscore |
|---|---|---|---|---|
| 1 | chr2L-1 | G1 (own group -- invisible) | 99.0% | 10,717 |
| 2 | chr8R-1 | **G3** | 97.7% | 6,623 |
| 3-7 | chr14L-3/4/5, chr13L-1/3 | **G2 (observed dominant donor)** | 97.3% | 6,549-6,565 |
| 8 | chr16R-1 | G6 | 97.6% | 6,307 |
| 9-35 | everything else | -- | <=97.5% | <=5,991 |

**G2 is the #2 detectable group by homology** (right behind G3, which only wins 1 of chr6L's 7
mismatches vs. G2's 4) -- most likely because G2 has 5 members to G3's 1, giving it more
chances to be the best-scoring match on any given read, the same "bigger group wins more often"
effect already seen for 6991's chr10L -> G8. **This reproduces the homology-driven pattern
found for chr16R, chr10L, and chr2L in 6991 -- chr6L -> G2 is not a chr13L-style outlier; it is
exactly what ordinary homology-driven strand invasion predicts**, and it recurs identically in
7172 (where G2 is the *only* donor group observed across all 4 mismatches).

## Master table: all 18 mismatched reads, global alignment vs. own reference and the full donor group

| strain | chr_end | read_id | own elem | own % | donor grp | donor elem | donor % | gap | verdict |
|---|---|---|---|---|---|---|---|---|---|
| 7302 | chr6L | SRR33298452.419935 | chr6L-1 | 73.71% | G10 | chr14R-1 | 98.47% | +24.76% | **whole-element swap** |
| 7302 | chr6L | SRR33298452.485661 | chr6L-1 | 75.92% | G8 | chr16L-1 | 99.17% | +23.26% | **whole-element swap** |
| 7172 | chr5R | SRR33298432.328531 | chr5R-1 | 80.20% | G2 | chr13L-1 | 96.53% | +16.33% | **whole-element swap** |
| 7302 | chr8R | SRR33298452.555243 | chr8R-1 | 89.29% | G1 | chr2L-1 | 98.91% | +9.63% | **whole-element swap** |
| 7302 | chr2L | SRR33298452.59861 | chr2L-1 | 90.93% | G3 | chr8R-1 | 99.76% | +8.83% | **whole-element swap** |
| 7302 | chr6L | SRR33298452.400002 | chr6L-1 | 90.96% | G2 | chr14L-3 | 98.78% | +7.82% | **whole-element swap** |
| 7302 | chr6L | SRR33298452.220014 | chr6L-1 | 89.09% | G3 | chr8R-1 | 96.58% | +7.49% | **whole-element swap** |
| 7302 | chr6L | SRR33298452.400093 | chr6L-1 | 87.00% | G2 | chr14L-3 | 94.31% | +7.32% | **whole-element swap** |
| 7302 | chr6L | SRR33298452.272227 | chr6L-1 | 90.24% | G2 | chr14L-3 | 97.56% | +7.31% | **whole-element swap** |
| 7172 | chr6L | SRR33298432.77931 | chr6L-1 | 88.92% | G2 | chr14L-3 | 95.79% | +6.87% | **whole-element swap** |
| 7302 | chr6L | SRR33298452.427945 | chr6L-1 | 90.54% | G2 | chr14L-3 | 97.34% | +6.79% | **whole-element swap** |
| 7302 | chr10L | SRR33298452.573441 | chr10L-1 | 94.07% | G10 | chr14R-1 | 98.47% | +4.39% | mid-Y' partial junction |
| 7302 | chr5R | SRR33298452.485827 | chr5R-1 | 91.18% | G7 | chr10L-1 | 93.35% | +2.16% | mid-Y' partial junction |
| 7302 | chr14R | SRR33298452.217955 | chr14R-1 | 96.71% | G8 | chr14L-1 | 98.35% | +1.64% | mid-Y' partial junction |
| 7302 | chr5R | SRR33298452.241882 | chr5R-1 | 96.65% | G10 | chr14R-1 | 98.18% | +1.52% | mid-Y' partial junction |
| 7302 | chr14R | SRR33298452.246202 | chr14R-1 | 97.04% | G8 | chr14L-1 | 98.18% | +1.14% | mid-Y' partial junction |
| 7172 | chr16R | SRR33298432.177252 | chr16R-1 | 97.55% | G2 | chr14L-3 | 96.68% | -0.87% | mid-Y' partial junction |
| 7172 | chr16R | SRR33298432.315476 | chr16R-1 | 97.66% | G2 | chr14L-3 | 96.10% | -1.56% | mid-Y' partial junction |

**11 of 18 (61%) are clean whole-element swaps** (gap > +5%, including several very large gaps
at chr6L -- +23% to +25% -- where the own reference barely explains the read at all: worth a
closer look at whether chr6L-1 is itself a good reference, separate from the recombination
question). **7 of 18 (39%) fall in the mid-Y' partial-junction band. 0 of 18 are unsupported** --
every single mismatched read in both strains shows a real, non-trivial signal.

## Position-enforced half-split test on the 7 mid-Y' candidates

Same rigorous test as the 6991 correction: search for the read's best-fitting split point,
then compare each half **only** to the proportionally corresponding half of each reference
(`edlib` global/NW alignment of matched-length slices -- no half can match an unrelated region
of the reference). A clean single breakpoint requires both halves decisive (>=5 points) in
**opposite** directions.

| chr_end | read_id | 1st half: own | 1st half: donor | margin | 2nd half: own | 2nd half: donor | margin | clean crossover? |
|---|---|---|---|---|---|---|---|---|
| chr16R | SRR33298432.177252 | 48.6% | 61.0% | +12.4 | 90.0% | 93.6% | +3.6 | no |
| chr16R | SRR33298432.315476 | 49.3% | 47.2% | -2.2 | 48.7% | 48.2% | -0.5 | no |
| chr10L | SRR33298452.573441 | 48.5% | 49.6% | +1.1 | 47.4% | 47.4% | +0.0 | no |
| chr14R | SRR33298452.217955 | 49.5% | 49.0% | -0.5 | 48.5% | 48.7% | +0.1 | no |
| chr14R | SRR33298452.246202 | 93.2% | 92.4% | -0.8 | 83.1% | 91.2% | +8.1 | no |
| chr5R | SRR33298452.485827 | 47.8% | 51.5% | +3.8 | 49.5% | 49.5% | +0.0 | no |
| chr5R | SRR33298452.241882 | 48.3% | 49.1% | +0.9 | 48.8% | 48.4% | -0.4 | no |

**0 of 7 pass.** This exactly reproduces the 6991 finding, including the specific pattern:
`SRR33298432.177252` (7172, chr16R -> G2) shows the same signature as 6991's chr16R reads --
own near-baseline (~49%) in the first half, donor moderately better (61%) but not clean, both
references high in the second half. That this specific asymmetric signature recurs
independently in a completely different strain and sequencing run is itself notable, but it is
still not a confirmed breakpoint under this test -- same caveat as 6991: this may reflect
non-Y' padding sequence pulled into the compared window rather than real internal structure,
and would need the padding narrowed and the test rerun to resolve.

## Cross-strain synthesis

| pattern | 6991 | 7172 | 7302 |
|---|---|---|---|
| chr2L -> G3 (chr8R) | 4 of 5 chr2L mismatches | not observed (no chr2L mismatches) | 1 of 1 |
| chr8R -> G1 (chr2L/chr6L) | not observed as a recipient pattern | not observed | 1 of 1 |
| chr16R -> G2 (chr13L/chr14L family) | 7 of 7 | 2 of 2 | not observed (no chr16R mismatches) |
| chr6L -> G2 (chr13L/chr14L family) | 1 of 1 | 1 of 1 | 4 of 7 (dominant, other 3 split across G3/G8/G10) |
| chr10L -> G8 (large tandem group) | 6 of 6 | not observed | not observed (7302's one chr10L mismatch went to G10/chr14R-1, unrelated) |
| chr5R -> various | not scored (no mismatch) | 1 of 1 -> G2 | 2 of 2 -> G7, G10 (no consistent donor) |
| chr13L -> G1 (chr2L/chr6L) | 15 of 18, dominant, **NOT homology-ranked #1** (the one confirmed outlier) | not observed | not observed |

**The two strongest, most consistent cross-strain patterns are chr2L<->chr8R (G1<->G3, mutual)
and chr6L/chr16R -> G2 (chr13L/chr14L family)** -- both recur in at least two of the three
strains with the same donor group every time, and both are explained by straightforward
sequence homology (chr6L's ranking above; chr16R and chr2L already established in the 6991
analysis). **chr13L -> G1 remains unique to 6991** among the patterns large enough to assess --
it does not recur in 7172 or 7302 (neither strain shows a chr13L mismatch at all in this
single-sample data), consistent with the earlier finding that it is driven by standing
population variation specific to that assembly's history rather than a general homology rule
that should be expected to recur everywhere.

## Files

- `verification/global_test_73xx.py`, `verification/half_split_73xx.py` -- the two scripts used
  (identical methodology to `global_identity_test.py` / `half_split_positional.py`, adapted for
  the two-strain, single-sample input shape).
- `verification/reports/cut99_7172_7302_data/` -- read fastas, element libraries, group JSONs,
  and the raw `master_*.tsv` / `half_split_*.tsv` results tables.
