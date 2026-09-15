# Day-0 reads whose Y' misses its reference group at the 99 % cutoff — 6991, 7172, 7302

Every read whose Y' copy count equals its anchor end's reference array length (so each copy has
a positional truth) was scored against that sample's own cut99 grouping. A mismatch is a copy
whose matched element falls in a different group from the element that positionally belongs
there. Each mismatch is then explained by splitting the read's Y' at the junction and aligning
BOTH halves to BOTH references.

| strain | scored reads | mismatches | strong | weak | FAILS | no junction | ref defect |
|---|---|---|---|---|---|---|---|
| 6991 (8 day-0 samples) | 79,593 | 530 | 52 | 22 | 30 | 38 | 388 |

6991 is eight independently assembled references, not eight replicates of one; see
`cut99_summary_6991_by_sample.md` for the per-reference breakdown. Its supported rate ranges
0.040-0.213 % across the eight, and 388 of its 530 mismatches come from one defective
reference.
| 7172 (1 sample) | 4,471 | 4 | 1 | 2 | 1 | 0 | 0 |
| 7302 (1 sample) | 7,112 | 14 | 8 | 3 | 1 | 2 | 0 |

Supported mid-Y' recombination (strong + weak), as a fraction of scored reads:

| strain | events | rate |
|---|---|---|
| 6991 | 74 | 0.093 % |
| 7172 | 3 | 0.067 % |
| 7302 | 11 | 0.155 % |

## What each label means

* **strong** — the anchor-side half favours the expected element AND the telomere-side half
  favours the donor, both by >= 1.5 %. This is the mid-Y' recombination signature.
* **weak** — same direction, margins under 1.5 %. These track low read quality as much as
  biology; treat as suggestive.
* **FAILS** — no half is better explained by the donor. The copy misses its group but a
  mid-Y' switch does not explain it. Unresolved.
* **no junction** — the donor wins across the whole element, so there is nothing to split.
  Note this is NOT necessarily a distinct mechanism: where donor and recipient are ~99 %
  identical, a junction anywhere produces the same sequence and cannot be localised.
* **reference defect** — 388 of 6991's 530 are a single artefact: `chr14L-1` is mis-assembled
  in `6991_day0_with_selection` (5,720 bp vs 6,654 elsewhere), so reads at that end match
  `chr7R-1` instead. Not recombination. Excluding it, 6991's mismatch rate falls from 0.67 %
  to 0.18 %, in line with the other strains.

## Caveats that apply to all three

* **These rates are floors.** A recombinant is only visible when donor and recipient land in
  different groups at this cutoff, and the largest group holds roughly half the elements
  (19 of 36 in 7302). Any switch within a group is invisible.
* **A native anchor is established; an untouched spacer and X element are not.** Reads were
  assigned to their end by anchor match, but the spacer and X element between anchor and Y'
  were never tested for recombination.
* **Library chimerism is not excluded.** A prep chimera would also join at homology. The
  telomere-side polarity of the junctions argues against it but does not rule it out.
* **The donor is usually a group, not an element** — the named donor is the best-scoring member
  of a near-identical set. Per-read alternatives are listed in each sample's files.

## Files

Per strain: `cut99_summary_<strain>.tsv` (full table) and `.md` (rendered).
Per sample under `cut99_<sample>/`: `groups.json` (the grouping), `mismatches.tsv` (raw),
`flank_identities.tsv` (per-half identities), `pair_homology.tsv` (crossover homology per pair).

Scripts: `build_cut99_groups.py`, `cut99_mismatch_report.py`, `scan_recombinant_junctions.py`,
`flank_identities.py`, `pair_homology.py`, `merge_cut99_report.py`.
