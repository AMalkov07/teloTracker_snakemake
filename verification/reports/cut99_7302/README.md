# 7302 day-0: Y' copies that don't match their reference group at the 99% cutoff

Grouping rebuilt with `verification/build_cut99_groups.py`, which reproduces the pipeline's
`--stop-mode threshold --identity-threshold 99` path exactly (homopolymer-condense,
all-vs-all blastn, coverage-penalised similarity, 99.9% dedup, average linkage,
fcluster at t = 100 - 99). For `7302_day0_with_selection`: **36 elements -> 12 groups**.

Scored: every read whose Y' copy count equals its anchor end's reference array length, so
copy *i* has a positional truth. **7,112 reads scored; 14 copies in 14 reads disagree with
their reference group (0.20%).** All 14 point to a different chromosome end.

Explanations come from `verification/scan_recombinant_junctions.py`, which slides a 300 bp
window along the read's Y' and scores each window against the expected element and the
observed one. A crossover -- one winning a run of windows at one end of the read, the other
at the other end -- is the recombinant signature.

| read | end | expected | observed | group | explanation |
|---|---|---|---|---|---|
| SRR33298452.573441 | chr10L | chr10L-1 | chr14R-1 | G7→G10 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.217955 | chr14R | chr14R-1 | chr14L-1 | G10→G8 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.246202 | chr14R | chr14R-1 | chr14L-1 | G10→G8 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.241882 | chr5R | chr5R-1 | chr14R-1 | G9→G10 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.485827 | chr5R | chr5R-1 | chr10L-1 | G9→G7 | recombination part-way through the Y'; direction not resolvable |
| SRR33298452.220014 | chr6L | chr6L-1 | chr8R-1 | G1→G3 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.272227 | chr6L | chr6L-1 | chr14L-5 | G1→G2 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.400002 | chr6L | chr6L-1 | chr14L-5 | G1→G2 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.400093 | chr6L | chr6L-1 | chr14L-5 | G1→G2 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.419935 | chr6L | chr6L-1 | chr14R-1 | G1→G10 | recombination part-way through the Y'; direction not resolvable |
| SRR33298452.485661 | chr6L | chr6L-1 | chr16L-1 | G1→G8 | recombination part-way through the Y'; donor on telomere side |
| SRR33298452.59861 | chr2L | chr2L-1 | chr8R-1 | G1→G3 | **whole Y' replaced** -- donor wins across the element, no junction |
| SRR33298452.555243 | chr8R | chr8R-1 | chr2L-1 | G3→G1 | **whole Y' replaced** -- donor wins across the element, no junction |
| SRR33298452.427945 | chr6L | chr6L-1 | chr14L-5 | G1→G2 | alternating signal: possible double crossover, or noise |

Summary: **11 partial (mid-Y') recombination, 2 whole-element replacement, 1 ambiguous.**
Of the 11, 9 resolve a direction and all 9 put the donor on the telomere side -- the BIR /
template-switch expectation. The other 2 have crossovers too diffuse to orient.

## The donor is a group, not an element

The "observed" element is only the best-scoring member of a set of near-identical sequences.
The real donor cannot be narrowed below its group:

| read | donor called | indistinguishable alternatives |
|---|---|---|
| SRR33298452.217955 / .246202 | chr14L-1 | 18 others in G8 (chr4R-1..7, chr12R-2..6, chr13L-2, chr13L-4, chr14L-2, chr15R-1, chr16L-1, chr7R-1) |
| SRR33298452.485661 | chr16L-1 | 18 others in G8 |
| SRR33298452.272227 / .400002 / .400093 / .427945 | chr14L-5 | chr13L-1, chr13L-3, chr14L-3, chr14L-4 |
| SRR33298452.485827 | chr10L-1 | chr9L-1 |
| SRR33298452.555243 | chr2L-1 | chr6L-1 |
| SRR33298452.573441 / .241882 / .419935 | chr14R-1 | (unique in its group) |
| SRR33298452.220014 / .59861 | chr8R-1 | (unique in its group) |

Four donors *are* uniquely identified (chr14R-1 and chr8R-1 are alone in G10 and G3).

## Caveats

* **A native anchor is established; an untouched spacer and X element are not.** Reads were
  assigned to their end by anchor match, and their Y' copy count matches the reference, but
  the spacer and X element between anchor and Y' were not tested for recombination.
* **Library chimerism is not excluded** -- a prep chimera would also join at homology. The
  telomere-side polarity argues against it but does not rule it out.
* **0.20% is a floor**, not a rate estimate: a recombinant is only visible when donor and
  recipient fall in different groups at this cutoff, and 19 of the 36 elements sit in G8.

Files: `cut99_mismatches_annotated.tsv` (this table plus window patterns),
`cut99_mismatches.tsv` (raw), `groups_7302.json` (the grouping).
