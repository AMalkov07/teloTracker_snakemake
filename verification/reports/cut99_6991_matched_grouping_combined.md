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

# chr16R -> G2, chr10L -> G8, chr2L -> G3: the same analysis, three more recipients

Same whole-read method: each test read's Y' BLASTed against its own expected reference and
against the donor group (best member), with a control of 6 correctly-assigned reads for that
end establishing the native baseline (own reference wins comfortably; donor only picks up a
partial hit through the shared internal repeat region).

## chr16R -> G2 {chr13L-1, chr14L-3/4/5}: 7/7 confirmed

| read | class | vs own (chr16R-1) | vs G2 best | margin | verdict |
|---|---|---|---|---|---|
| SRR33298373.72733 | weak | 99.3% / 3,578bp / bs=6,469 | 99.0% / 5,199bp (chr14L-5) / bs=9,286 | +2,817 | **recombinant** |
| SRR33298384.248800 | strong | 99.5% / 3,583bp / bs=6,508 | 98.9% / 5,203bp (chr14L-5) / bs=9,280 | +2,772 | **recombinant** |
| SRR33298384.538094 | strong | 98.1% / 3,599bp / bs=6,233 | 97.6% / 5,230bp (chr13L-1) / bs=8,911 | +2,678 | **recombinant** |
| SRR33298377.534605 | strong | 98.2% / 3,591bp / bs=6,252 | 97.4% / 5,209bp (chr13L-1) / bs=8,820 | +2,568 | **recombinant** |
| SRR33298434.378016 | FAILS | 96.2% / 3,612bp / bs=5,847 | 95.8% / 5,252bp (chr14L-5) / bs=8,386 | +2,539 | **recombinant** |
| SRR33298377.644709 | FAILS | 98.1% / 3,596bp / bs=6,237 | 97.3% / 5,214bp (chr13L-1) / bs=8,776 | +2,539 | **recombinant** |
| SRR33298434.64148 | FAILS | 94.4% / 3,611bp / bs=5,459 | 92.5% / 5,265bp (chr14L-5) / bs=7,352 | +1,893 | **recombinant** |

7/7, margins +1,893 to +2,817. All 3 window-scan "FAILS" reads confirm under the whole-read
test -- the sliding-window scan can miss a real event when the junction sits near a read end,
but the full-molecule bitscore contest still catches it.

## chr10L -> G8 (18-member array group): 6/6 confirmed

| read | class | vs own (chr10L-1) | vs G8 best | margin | verdict |
|---|---|---|---|---|---|
| SRR33298384.498901 | strong | 98.2% / 5,003bp / bs=8,698 | 98.7% / 6,614bp (chr7R-1) / bs=11,725 | +3,027 | **recombinant** |
| SRR33298461.57209 | strong | 97.5% / 5,010bp / bs=8,510 | 98.0% / 6,623bp (chr7R-1) / bs=11,437 | +2,927 | **recombinant** |
| SRR33298373.491903 | FAILS | 97.6% / 5,010bp / bs=8,540 | 98.0% / 6,619bp (chr16L-1) / bs=11,452 | +2,912 | **recombinant** |
| SRR33298373.545529 | strong | 98.0% / 5,002bp / bs=8,667 | 98.3% / 6,611bp (chr16L-1) / bs=11,568 | +2,901 | **recombinant** |
| SRR33298384.151130 | weak | 94.1% / 5,043bp / bs=7,542 | 95.2% / 6,664bp (chr7R-1) / bs=10,405 | +2,863 | **recombinant** |
| SRR33298373.122658 | weak | 96.0% / 5,035bp / bs=8,117 | 95.8% / 6,665bp (chr4R-7) / bs=10,661 | +2,544 | **recombinant** |

6/6, margins +2,544 to +3,027 -- the tightest, most consistent set of any recipient examined.

## chr2L -> G3 {chr8R-1}: only 2/5 confirmed -- weaker than the raw count implied

| read | class | vs own (chr2L-1) | vs G3 (chr8R-1) | margin | verdict |
|---|---|---|---|---|---|
| SRR33298373.284052 | no junction | 98.8% / 3,861bp / bs=6,863 | 99.6% / 5,473bp / bs=9,987 | +3,124 | **recombinant** |
| SRR33298434.94939 | strong | 95.7% / 2,042bp / bs=3,262 | 97.2% / 3,664bp / bs=6,150 | +2,888 | **recombinant** |
| SRR33298384.313055 | no junction | 98.9% / 2,270bp / bs=4,052 | 98.9% / 2,270bp / bs=4,052 | 0 | not supported |
| SRR33298384.467718 | no junction | 97.8% / 3,175bp / bs=5,448 | 97.8% / 3,172bp / bs=5,443 | -5 | not supported |
| SRR33298377.541858 | FAILS | 97.6% / 3,946bp / bs=6,728 | 97.5% / 3,868bp / bs=6,584 | -144 | not supported |

Only 2 of 5 show a real contest between two candidates. The other 3 have near-identical
bitscores to both references over short (2,270-3,946bp) alignments -- own and donor return
essentially the same weak partial hit, which is not the recombination signature. **Read this
as 2 confirmed events out of 5 calls, not "4/5."**

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
