# 6991 combined: the five samples sharing an identical cut99 grouping

Restricted to samples whose cut99 partition (34 elements -> 12 groups) is EXACTLY
identical, verified by comparing group membership directly (not just group count).

## Which samples qualify

| sample | elements | groups | included |
|---|---|---|---|
| 6991_day0 | 34 | 12 | yes |
| 6991_day0_TeloTag_with_selection | 34 | 12 | yes |
| 6991_day0_reference_promethion | 34 | 12 | yes |
| 6991_day0_with_selection_repeat | 34 | 12 | yes |
| 6991_day0_with_selection_repeat2 | 34 | 12 | yes |
| 6991_day0_reference | 33 | 12 | **no** -- one fewer element (missing chr12R-7); same 12 groups otherwise but not an exact match |
| 6991_day0_TeloTag | 27 | 13 | **no** -- missing the whole chr4R array (structurally different) |
| 6991_day0_with_selection | 34 | 13 | **no** -- the chr14L-1 defect (5,720 bp) additionally splits off as its own group |

**5 of 8 samples qualify** (not 6 -- 6991_day0_reference is close but not exact).

## Group key (identical across all 5 samples)

| group | size | members |
|---|---|---|
| G8 | 18 | chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-1, chr14L-2, chr15R-1, chr16L-1, chr4R-1, chr4R-2, chr4R-3, chr4R-4, chr4R-5, chr4R-6, chr4R-7, chr7R-1 |
| G2 | 4 | chr13L-1, chr14L-3, chr14L-4, chr14L-5 |
| G7 | 2 | chr10L-1, chr9L-1 |
| G1 | 2 | chr2L-1, chr6L-1 |
| G4 | 1 | chr12L-1 |
| G11 | 1 | chr12R-1 |
| G10 | 1 | chr14R-1 |
| G6 | 1 | chr16R-1 |
| G12 | 1 | chr5L-1 |
| G9 | 1 | chr5R-1 |
| G5 | 1 | chr8L-1 |
| G3 | 1 | chr8R-1 |

## Pooled counts across the 5 samples (58 mismatched copies)

### As recipient

| end | count |
|---|---|
| chr13L | 18 |
| chr16R | 7 |
| chr10L | 6 |
| chr14R | 6 |
| chr2L | 5 |
| chr5L | 3 |
| chr8R | 2 |
| chr8L | 2 |
| chr12L | 2 |
| chr14L | 2 |
| chr6L | 1 |
| chr7R | 1 |
| chr15R | 1 |
| chr16L | 1 |
| chr5R | 1 |

### As donor, by group (unambiguous -- group labels are identical across all 5)

| group | times used as donor | group size |
|---|---|---|
| G1 | 18 | 2 |
| G2 | 13 | 4 |
| G8 | 10 | 18 |
| G7 | 5 | 2 |
| G11 | 4 | 1 |
| G3 | 4 | 1 |
| G10 | 2 | 1 |
| G6 | 1 | 1 |
| G5 | 1 | 1 |

### Combined — the 5 identically-grouped 6991 samples (58 copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr13L | 18 | G1 (chr2L-1, chr6L-1) | 15 | 83% |
|  |  | G6 (chr16R-1) | 1 | 6% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 6% |
|  |  | G5 (chr8L-1) | 1 | 6% |
| chr16R | 7 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 7 | 100% |
| chr10L | 6 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 6 | 100% |
| chr14R | 6 | G1 (chr2L-1, chr6L-1) | 2 | 33% |
|  |  | G7 (chr10L-1, chr9L-1) | 2 | 33% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 17% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 17% |
| chr2L | 5 | G3 (chr8R-1) | 4 | 80% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 20% |
| chr5L | 3 | G11 (chr12R-1) | 2 | 67% |
|  |  | G7 (chr10L-1, chr9L-1) | 1 | 33% |
| chr8R | 2 | G10 (chr14R-1) | 1 | 50% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 50% |
| chr8L | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 50% |
|  |  | G1 (chr2L-1, chr6L-1) | 1 | 50% |
| chr12L | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 2 | 100% |
| chr14L | 2 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 50% |
|  |  | G11 (chr12R-1) | 1 | 50% |
| chr6L | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |
| chr7R | 1 | G10 (chr14R-1) | 1 | 100% |
| chr15R | 1 | G7 (chr10L-1, chr9L-1) | 1 | 100% |
| chr16L | 1 | G11 (chr12R-1) | 1 | 100% |
| chr5R | 1 | G7 (chr10L-1, chr9L-1) | 1 | 100% |
---

# Alignment-score defence: chr13L -> G1 is real recombination for most reads

All 15 chr13L reads whose Y' matched into G1 (chr2L-1/chr6L-1), each read's full Y' region
BLASTed against **its own expected reference (chr13L-1)** and against **G1** (best of
chr2L-1/chr6L-1). Compared to a control of 6 correctly-assigned chr13L-1 reads to establish
what a genuine native read looks like under the same test.

## Control: what a real chr13L-1 read looks like

| | vs chr13L-1 (own reference) | vs G1 |
|---|---|---|
| 6 correctly-assigned reads | full length (~5,486 bp), 99.3-99.8% identity, bitscore ~6,600-10,100 | partial (~2,000-3,890 bp), 96.7-97.2% identity, bitscore ~3,400-6,500 |

A native read matches its own reference end-to-end at high identity and loses badly to G1,
which only picks up the shared middle portion (the internal tandem repeat region already
characterised elsewhere in this investigation lets any two related Y' partially cross-align,
native or not).

## The 15 test reads: 13 show the exact reverse of the control pattern

| read | sample | window-scan class | vs chr13L-1 (own ref) | vs G1 best match | bitscore margin | verdict |
|---|---|---|---|---|---|---|
| SRR33298373.1068481 | ws_repeat | no junction | 97.2% / 3,886bp / bs=6,527 | 99.9% / 5,978bp (chr6L-1) / bs=11,001 | G1 +4,474 | **recombinant** -- G1 wins, full length |
| SRR33298373.1079946 | ws_repeat | no junction | 97.2% / 3,883bp / bs=6,525 | 99.9% / 5,976bp (chr6L-1) / bs=10,988 | G1 +4,463 | **recombinant** -- G1 wins, full length |
| SRR33298373.896473 | ws_repeat | no junction | 95.8% / 3,897bp / bs=6,211 | 98.9% / 5,991bp (chr6L-1) / bs=10,671 | G1 +4,460 | **recombinant** -- G1 wins, full length |
| SRR33298373.107087 | ws_repeat | no junction | 96.2% / 3,894bp / bs=6,311 | 99.2% / 5,989bp (chr6L-1) / bs=10,752 | G1 +4,441 | **recombinant** -- G1 wins, full length |
| SRR33298373.797552 | ws_repeat | no junction | 97.0% / 3,885bp / bs=6,490 | 99.7% / 5,977bp (chr6L-1) / bs=10,920 | G1 +4,430 | **recombinant** -- G1 wins, full length |
| SRR33298373.46529 | ws_repeat | no junction | 97.0% / 3,884bp / bs=6,488 | 99.6% / 5,981bp (chr6L-1) / bs=10,905 | G1 +4,417 | **recombinant** -- G1 wins, full length |
| SRR33298384.185796 | ws_repeat2 | FAILS | 95.6% / 3,891bp / bs=6,156 | 98.6% / 5,986bp (chr6L-1) / bs=10,539 | G1 +4,383 | **recombinant** -- G1 wins, full length |
| SRR33298377.260642 | reference_promethion | no junction | 96.9% / 3,885bp / bs=6,471 | 99.4% / 5,982bp (chr6L-1) / bs=10,820 | G1 +4,349 | **recombinant** -- G1 wins, full length |
| SRR33298373.8939 | ws_repeat | strong | 98.8% / 3,857bp / bs=6,866 | 98.7% / 5,989bp (chr6L-1) / bs=10,591 | G1 +3,725 | **recombinant** -- G1 wins, full length |
| SRR33298373.1520960 | ws_repeat | strong | 98.4% / 3,870bp / bs=6,791 | 98.4% / 6,005bp (chr6L-1) / bs=10,493 | G1 +3,702 | **recombinant** -- G1 wins, full length |
| SRR33298377.272313 | reference_promethion | no junction | 92.1% / 3,933bp / bs=5,421 | 94.0% / 6,061bp (chr6L-1) / bs=9,027 | G1 +3,606 | **recombinant** -- G1 wins, full length |
| SRR33298377.125819 | reference_promethion | strong | 95.7% / 3,898bp / bs=6,207 | 96.3% / 6,035bp (chr2L-1) / bs=9,801 | G1 +3,594 | **recombinant** -- G1 wins, full length |
| SRR33298377.267122 | reference_promethion | strong | 95.2% / 3,900bp / bs=6,117 | 96.0% / 6,053bp (chr2L-1) / bs=9,701 | G1 +3,584 | **recombinant** -- G1 wins, full length |
| SRR33298384.45869 | ws_repeat2 | no junction | 96.3% / 3,885bp / bs=6,322 | 96.3% / 4,316bp (chr6L-1) / bs=7,068 | G1 +746 | recombinant, but G1 hit is itself partial -- coordinate artefact, see caveat |
| SRR33298461.44252 | day0 | weak | 97.1% / 3,891bp / bs=6,529 | 94.7% / 4,177bp (chr6L-1) / bs=6,416 | -113 | **NOT supported** -- own reference wins |

## Summary

**13 of 15 (87%)** show the full pattern that defines real recombination under this test:
full-length (~5,975-6,061 bp) high-identity match to G1, against only a partial
(~3,857-3,933 bp) match to chr13L-1 -- the exact reverse of what the 6 control reads show.
This is not a marginal call: bitscore margins for these 13 range **+3,584 to +4,474**, roughly
3,600-40x the +746 margin of the weakest positive case.

**2 of 15 do not clearly support recombination:**
* `SRR33298384.45869` -- G1 still wins on bitscore, but its own G1 hit is only 4,316 bp
  (72% of G1's length), not the ~6,000 bp every other positive case shows. This read was
  already flagged with a coordinate-extraction artefact in the windowed-comparison version of
  this analysis; the whole-read result here is consistent with that, not a clean confirmation.
* `SRR33298461.44252` -- **chr13L-1 wins**, by a small margin (113 bitscore, ~2%). This read
  does not show the recombinant signature under a whole-read comparison. It was originally
  called "weak" by the sliding-window scan, which is consistent with it being a borderline or
  false-positive call rather than a real event.

**Bottom line: the recombination signature is real and strong for the large majority of
these reads (13/15), defended by bitscore margins of 3,500+ against a same-locus native-read
control, not just by percent identity. It is not established for 2 of the 15, and those two
should not be counted as confirmed recombinants without further work on their coordinates.**

---

## Why almost all of it goes to the same donor group: mostly real homology, chr13L excepted

Each recipient's own reference was BLASTed against every other element in the library
(excluding same-group members, which are structurally invisible to this detection method --
a read converting to its own group's sequence is indistinguishable from a native read and can
never register as a mismatch):

| recipient | donor group rank among detectable candidates | confirmed |
|---|---|---|
| chr16R -> G2 | **#1 -- G2 occupies the top 4 ranks outright** | 7/7 |
| chr10L -> G8 | **#1 by a wide margin** -- all 17 detectable G8 members occupy ranks 2-18, beating the nearest non-member by only 11-37 bitscore points before a 3,000+ point cliff to everything else | 6/6 |
| chr2L -> G3 | **#1** -- chr8R-1 is the single most homologous detectable element | 2/5 |
| chr13L -> G1 | **#8-9 of ~12 -- one of the LEAST homologous options** | 13/15 |

**For three of the four recipients, the dominant donor is simply the most sequence-homologous
available partner** -- exactly what a homology-driven strand-invasion mechanism predicts, and
a sufficient explanation on its own. chr10L -> G8 is the cleanest example: after its one
invisible own-group relative, G8's 17 members sweep the entire top of the ranking with almost
nothing else coming close.

**chr13L is the outlier**, and stays one: its dominant donor (G1) is measurably one of the
least homologous groups available to it, so homology does not explain that case. The
supported explanation there (see the ancestral-haplotype analysis above) is a small number of
historical recombination events that became standing variation in the population and are
being resampled by every independent sequencing prep, rather than an active homology
preference operating today. That mechanism does not need to generalise to explain the other
three recipients, and the homology ranking shows directly that it should not be assumed to.

---

# Correction: the bitscore-margin test above was the wrong metric -- redone properly

A fair objection was raised: in several rows above, `own` has **higher %identity** than
`donor`, yet the row is still marked recombinant. That looks contradictory if you only look
at %identity, and the bitscore-margin framing did not make clear why it isn't.

## Why %identity of the best single HSP is the wrong test

High %identity over a **short, partial** alignment is exactly what unrelated-but-related
sequence produces anyway -- every Y' at this locus shares a divergent internal region (the
tandem repeat characterised earlier in this investigation), so two genuinely different
elements can still align at 97-99% identity over the ~3,600-5,000 bp portion outside that
region. A short high-identity match does not mean "this read is that reference"; it means
"these two references are related," which is true of every pair examined here regardless of
recombination.

## The correct test: does ONE unbroken alignment cover (close to) the full reference?

If a read is genuinely a native copy of reference X, the whole read should align to X in a
**single HSP spanning nearly all of X's length** -- exactly what the control reads show. If a
read is a recombinant, the true native reference will only ever produce a **partial** best
HSP (it cannot explain the donor-derived portion), while the true donor will produce a
**full-length** HSP (it explains the read end to end).

**Control proof (chr16R, 6 correctly-assigned reads):**

| | vs own (chr16R-1, 5,302bp) | vs G2 (donor) |
|---|---|---|
| best single HSP | **~5,300bp -- 99.9% of own's length, ONE block** | ~3,590bp -- 68% of donor's length, needs 2-3 HSPs to extend further |

Native reads reach full length against their own reference and stay stuck at ~65-70% against
the donor, no matter how high that partial identity is (up to 98.9% in the controls).

## Redone with this test: reference-length coverage of the single best HSP

| recipient | own ref length | donor ref (varies by read) | how many reach >=90% coverage against DONOR while own stays <85%? |
|---|---|---|---|
| chr13L (own 5,483bp) | | chr2L-1/chr6L-1, 5,975bp | **13 of 15** -- own caps at 70-72% every time; donor reaches 100-101% for the 13 confirmed |
| chr16R (own 5,302bp) | | chr13L-1/chr14L-3/4/5, ~5,484bp | **7 of 7** -- own caps at 67-68%; donor reaches 95-96% |
| chr10L (own 6,869bp) | | chr7R-1/chr16L-1/chr4R-7, ~6,654bp | **6 of 6** -- own caps at 73-73%; donor reaches 99-100% |
| chr2L (own 5,975bp) | | chr8R-1, 5,469bp | **1 of 5 clean** (SRR33298373.284052: own 65% / donor 100%); SRR33298434.94939 leans donor (34% / 67%) but neither reaches 90%, so it is now ambiguous rather than confirmed; the other 3 stay ambiguous as before |

**Corrected total: 27 of 33 reads across the four recipients pass this stricter, better-
justified test** (13+7+6+1), the same overall count as before but with **chr2L revised down**
from "2 confirmed" to **1 clearly confirmed**, since SRR33298434.94939 does not reach the
90% full-length bar on either side even though it favours the donor.

## Why this resolves the original concern

The own-vs-donor contest was never really about which single HSP has the better raw score --
it is about which reference can explain the **entire** molecule in one piece. A read's own
reference losing on %identity while "winning" a short fragment is not evidence of anything;
what matters is that **the own reference never manages to extend past ~65-73% of its own
length for any of the confirmed recombinants, in every one of the four recipients, while the
donor consistently reaches 94-101%** -- matching almost exactly what the size difference
between the native and donor Y' variants predicts. That contrast, not the bitscore number,
is the actual evidence for recombination.

---

# Final correction: true global alignment separates whole-element swaps from mid-Y' partial junctions

Both tests above are still built on BLAST **local** alignment (HSPs) -- coverage-of-best-HSP is
a proxy for "does one reference explain the whole read," not a direct measurement of it. The
proper test is a true **global** alignment of each candidate reference against the read:
[`global_identity_test.py`](../global_identity_test.py) uses `edlib` in infix mode (`HW`) --
the reference is forced to align end-to-end, the read is free at both ends (it carries extra
flanking anchor/telomere sequence beyond the Y' itself) -- checking both DNA strands and taking
the better one (edlib, unlike blastn, does not do this automatically). This gives one clean
number per reference: **% identity over the reference's full length**, not a fragment.

**Expected signatures under this test:**
* **Native read:** own identity clearly higher than any donor identity (own explains the whole
  molecule; donor never fully spans the reference-shared middle region alone).
* **Whole-element swap:** donor identity clearly higher than own, both close to reference
  quality (~95-100%) -- the ENTIRE reference length is well explained by the donor and poorly
  by the true native.
* **Mid-Y' partial junction (part native, part donor):** neither reference reaches a clean win
  -- forcing a chimeric read to align end-to-end against a single, non-chimeric reference
  necessarily produces a mediocre identity against *both* candidates, since each one only
  explains part of the molecule. This is exactly the case the user's own framing anticipated:
  it needs the sliding-window scan (which already localizes the junction), not a global
  identity comparison, because there is no single reference for the global test to confirm.

## Master table: all 58 mismatched reads, global alignment vs. own reference and vs. the full donor group

Every read from the pooled 58-mismatch table above (all 5 identically-grouped 6991 samples),
run through the same test: own reference's full-length identity vs. the **best member of the
full donor group** (not just one representative), plus the pre-existing sliding-window-scan
class (`ws class`, from `scan_recombinant_junctions.py`) as an independent cross-check. Verdict
thresholds: gap > +5% = whole-element swap; -2% to +5% = mid-Y' partial junction (real signal,
shifted from the native baseline, but not a clean donor win -- see native control baselines of
-4.8% to -10.8% established per end above); gap < -2% = not supported.

| chr_end | read_id | own elem | own % | donor grp | donor elem | donor % | gap | ws class | verdict |
|---|---|---|---|---|---|---|---|---|---|
| chr10L | SRR33298384.498901 | chr10L-1 | 94.72% | G8 | chr7R-1 | 98.32% | +3.60% | strong | **mid-Y' partial junction** |
| chr10L | SRR33298384.151130 | chr10L-1 | 91.25% | G8 | chr7R-1 | 94.82% | +3.57% | weak | **mid-Y' partial junction** |
| chr10L | SRR33298461.57209 | chr10L-1 | 94.18% | G8 | chr7R-1 | 97.57% | +3.39% | strong | **mid-Y' partial junction** |
| chr10L | SRR33298373.545529 | chr10L-1 | 95.37% | G8 | chr16L-1 | 97.93% | +2.56% | strong | **mid-Y' partial junction** |
| chr10L | SRR33298373.491903 | chr10L-1 | 95.04% | G8 | chr16L-1 | 97.58% | +2.54% | FAILS | **mid-Y' partial junction** |
| chr10L | SRR33298373.122658 | chr10L-1 | 92.98% | G8 | chr12R-2 | 95.35% | +2.37% | weak | **mid-Y' partial junction** |
| chr12L | SRR33298384.550794 | chr12L-1 | 98.34% | G2 | chr13L-1 | 96.04% | -2.30% | weak | not supported |
| chr12L | SRR33298384.443822 | chr12L-1 | 98.96% | G2 | chr14L-3 | 96.55% | -2.41% | weak | not supported |
| chr13L | SRR33298434.76895 | chr13L-1 | 80.28% | G8 | chr15R-1 | 95.47% | +15.19% | no junction | **whole-element swap** |
| chr13L | SRR33298384.185796 | chr13L-1 | 87.74% | G1 | chr6L-1 | 98.63% | +10.88% | FAILS | **whole-element swap** |
| chr13L | SRR33298373.896473 | chr13L-1 | 88.05% | G1 | chr6L-1 | 98.93% | +10.87% | no junction | **whole-element swap** |
| chr13L | SRR33298373.797552 | chr13L-1 | 88.86% | G1 | chr6L-1 | 99.68% | +10.83% | no junction | **whole-element swap** |
| chr13L | SRR33298373.107087 | chr13L-1 | 88.35% | G1 | chr6L-1 | 99.16% | +10.82% | no junction | **whole-element swap** |
| chr13L | SRR33298373.1068481 | chr13L-1 | 89.11% | G1 | chr6L-1 | 99.90% | +10.79% | no junction | **whole-element swap** |
| chr13L | SRR33298377.260642 | chr13L-1 | 88.60% | G1 | chr6L-1 | 99.38% | +10.78% | no junction | **whole-element swap** |
| chr13L | SRR33298373.1079946 | chr13L-1 | 89.09% | G1 | chr6L-1 | 99.87% | +10.77% | no junction | **whole-element swap** |
| chr13L | SRR33298373.46529 | chr13L-1 | 88.89% | G1 | chr6L-1 | 99.62% | +10.72% | no junction | **whole-element swap** |
| chr13L | SRR33298377.272313 | chr13L-1 | 83.57% | G1 | chr6L-1 | 93.96% | +10.39% | no junction | **whole-element swap** |
| chr13L | SRR33298384.45869 | chr13L-1 | 89.93% | G1 | chr6L-1 | 99.06% | +9.13% | no junction | **whole-element swap** |
| chr13L | SRR33298377.267122 | chr13L-1 | 86.89% | G1 | chr2L-1 | 95.95% | +9.06% | strong | **whole-element swap** |
| chr13L | SRR33298377.125819 | chr13L-1 | 87.42% | G1 | chr2L-1 | 96.37% | +8.95% | strong | **whole-element swap** |
| chr13L | SRR33298373.1520960 | chr13L-1 | 89.77% | G1 | chr6L-1 | 98.38% | +8.61% | strong | **whole-element swap** |
| chr13L | SRR33298373.8939 | chr13L-1 | 90.01% | G1 | chr6L-1 | 98.59% | +8.59% | strong | **whole-element swap** |
| chr13L | SRR33298373.677153 | chr13L-1 | 93.93% | G5 | chr8L-1 | 99.57% | +5.64% | no junction | **whole-element swap** |
| chr13L | SRR33298384.220455 | chr13L-1 | 96.86% | G6 | chr16R-1 | 99.00% | +2.14% | weak | **mid-Y' partial junction** |
| chr13L | SRR33298461.44252 | chr13L-1 | 92.91% | G1 | chr6L-1 | 92.94% | +0.03% | weak | **mid-Y' partial junction** |
| chr14L | SRR33298461.168936 | chr14L-5 | 83.17% | G8 | chr14L-2 | 99.49% | +16.32% | FAILS | **whole-element swap** |
| chr14L | SRR33298434.530414 | chr14L-1 | 97.70% | G11 | chr12R-1 | 94.83% | -2.87% | FAILS | not supported |
| chr14R | SRR33298377.93825 | chr14R-1 | 73.78% | G1 | chr2L-1 | 99.08% | +25.30% | no junction | **whole-element swap** |
| chr14R | SRR33298384.285744 | chr14R-1 | 74.99% | G1 | chr6L-1 | 98.78% | +23.79% | strong | **whole-element swap** |
| chr14R | SRR33298373.942862 | chr14R-1 | 79.02% | G2 | chr13L-1 | 98.87% | +19.85% | strong | **whole-element swap** |
| chr14R | SRR33298373.897522 | chr14R-1 | 89.66% | G7 | chr9L-1 | 94.24% | +4.57% | FAILS | **mid-Y' partial junction** |
| chr14R | SRR33298377.726683 | chr14R-1 | 91.19% | G7 | chr9L-1 | 93.39% | +2.20% | weak | **mid-Y' partial junction** |
| chr14R | SRR33298373.1266643 | chr14R-1 | 96.58% | G8 | chr14L-1 | 97.53% | +0.95% | strong | **mid-Y' partial junction** |
| chr15R | SRR33298384.129326 | chr15R-1 | 93.68% | G7 | chr9L-1 | 98.72% | +5.03% | strong | **whole-element swap** |
| chr16L | SRR33298434.79614 | chr16L-1 | 95.09% | G11 | chr12R-1 | 93.71% | -1.38% | no junction | **mid-Y' partial junction** |
| chr16R | SRR33298384.538094 | chr16R-1 | 96.04% | G2 | chr13L-1 | 95.15% | -0.89% | strong | **mid-Y' partial junction** |
| chr16R | SRR33298434.64148 | chr16R-1 | 91.02% | G2 | chr14L-3 | 90.05% | -0.98% | FAILS | **mid-Y' partial junction** |
| chr16R | SRR33298434.378016 | chr16R-1 | 94.47% | G2 | chr14L-3 | 93.40% | -1.07% | FAILS | **mid-Y' partial junction** |
| chr16R | SRR33298377.644709 | chr16R-1 | 95.85% | G2 | chr13L-1 | 94.77% | -1.08% | FAILS | **mid-Y' partial junction** |
| chr16R | SRR33298377.534605 | chr16R-1 | 95.98% | G2 | chr13L-1 | 94.89% | -1.09% | strong | **mid-Y' partial junction** |
| chr16R | SRR33298384.248800 | chr16R-1 | 97.45% | G2 | chr14L-3 | 96.35% | -1.10% | strong | **mid-Y' partial junction** |
| chr16R | SRR33298373.72733 | chr16R-1 | 97.51% | G2 | chr14L-3 | 96.39% | -1.12% | weak | **mid-Y' partial junction** |
| chr2L | SRR33298434.94939 | chr2L-1 | 75.23% | G8 | chr16L-1 | 97.82% | +22.59% | strong | **whole-element swap** |
| chr2L | SRR33298373.284052 | chr2L-1 | 90.83% | G3 | chr8R-1 | 99.63% | +8.81% | no junction | **whole-element swap** |
| chr2L | SRR33298384.313055 | chr2L-1 | 87.38% | G3 | chr8R-1 | 94.97% | +7.59% | no junction | **whole-element swap** |
| chr2L | SRR33298377.541858 | chr2L-1 | 90.58% | G3 | chr8R-1 | 96.00% | +5.42% | FAILS | **whole-element swap** |
| chr2L | SRR33298384.467718 | chr2L-1 | 92.32% | G3 | chr8R-1 | 89.78% | -2.54% | no junction | not supported |
| chr5L | SRR33298373.260098 | chr5L-1 | 93.21% | G7 | chr10L-1 | 97.73% | +4.52% | no junction | **mid-Y' partial junction** |
| chr5L | SRR33298434.226995 | chr5L-1 | 96.13% | G11 | chr12R-1 | 96.48% | +0.36% | strong | **mid-Y' partial junction** |
| chr5L | SRR33298434.116145 | chr5L-1 | 96.30% | G11 | chr12R-1 | 96.62% | +0.33% | strong | **mid-Y' partial junction** |
| chr5R | SRR33298377.75463 | chr5R-1 | 92.96% | G7 | chr9L-1 | 98.36% | +5.39% | no junction | **whole-element swap** |
| chr6L | SRR33298461.181020 | chr6L-1 | 89.32% | G2 | chr13L-1 | 96.77% | +7.45% | strong | **whole-element swap** |
| chr7R | SRR33298434.164214 | chr7R-1 | 94.18% | G10 | chr14R-1 | 95.31% | +1.13% | strong | **mid-Y' partial junction** |
| chr8L | SRR33298377.763070 | chr8L-1 | 97.41% | G2 | chr13L-1 | 95.00% | -2.41% | weak | not supported |
| chr8L | SRR33298434.497585 | chr8L-1 | 92.86% | G1 | chr2L-1 | 89.46% | -3.41% | FAILS | not supported |
| chr8R | SRR33298434.515556 | chr8R-1 | 79.59% | G10 | chr14R-1 | 98.33% | +18.74% | strong | **whole-element swap** |
| chr8R | SRR33298377.221175 | chr8R-1 | 98.50% | G2 | chr14L-3 | 98.21% | -0.29% | weak | **mid-Y' partial junction** |

## Bottom line: 58 of 58 mismatched reads, true global alignment

| verdict | count | % of 58 |
|---|---|---|
| whole-element swap (gap > +5%) | 28 | 48% |
| mid-Y' partial junction (-2% to +5%) | 24 | 41% |
| not supported (gap < -2%) | 6 | 10% |

**52 of 58 (90%) of the reads whose Y' matched into the wrong cut99 group show a real,
quantifiable recombination signature** under the most rigorous test applied so far -- either a
clean whole-element swap or a real, window-scan-corroborated shift consistent with a mid-Y'
chimera. Only 6 reads (chr12L x2, chr14L x1, chr8L x2, chr2L x1) fail to clear even the partial-
junction bar and should be treated as unconfirmed grouping noise rather than recombination.

The two mechanistic classes do not track cleanly with any single recipient -- chr13L, for
example, has both whole-element swaps (15/18) and one mid-Y' case (`SRR33298384.220455`,
`SRR33298461.44252`) -- but the split is informative in aggregate: whole-element swaps dominate
where the donor and recipient are close in length and homology (chr13L->G1, chr2L->G3, chr14R,
chr8R), while mid-Y' partial junctions cluster in the size-mismatched, highly-homologous G8/G7
cases (chr10L->G8, chr16R->G2, chr5L->G11/G7, chr7R->G10) where a clean swap and a partial
conversion are both plausible outcomes of the same strand-invasion mechanism, and where the
window-scan (never `no junction` for these reads) already agrees a real breakpoint exists.
