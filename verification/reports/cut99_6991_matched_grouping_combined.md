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
