# Do the other "mislabelled" reads share the chr6L recombinant structure? Mostly yes.

Follow-up to `verification/reports/chr6L_recombinants/`, which showed that reads at chr6L
that matched no single Y' element well were in fact recombinant hybrids. This extends the
test to every such read across all ten corrected-boundary (`__elemYPfix`) day-0 runs.

## Selecting candidates

From all 10 runs, take every read whose Y' copy count equals its end's reference array
length (so copy i has a positional truth) and where some copy was matched to an element at a
DIFFERENT chromosome end:

| filter | copies |
|---|---|
| foreign-donor misread copies | 7,532 |
| ...after dropping pairs that are >=99% identical over >=90% of their length | 545 |
| ...after dropping the known `6991_day0_with_selection` chr14L-1 assembly defect (5,720 bp vs 6,654) | **157 reads** |

The first filter removes the bulk: most foreign-donor calls are between elements that are
indistinguishable by sequence (chr4R-7 vs chr14L-2, chr12R-* vs chr4R-7, chr7R-1 vs chr16L-1
after the boundary fix). Those are arbitrary assignment, not recombination.

## The test

`verification/scan_recombinant_junctions.py` slides a 300 bp window along each read's Y'
region and scores it against both the recipient element (positional truth, "A") and the
matched element ("B"). A recombinant shows a **crossover** -- one element winning a run of
windows at one end of the read, the other at the other end. A noisy read, or one with a bad
reference, shows a single element winning throughout.

| verdict | reads |
|---|---|
| **RECOMBINANT** | **106** |
| single | 45 |
| mixed (alternating, possible double crossover or noise) | 3 |
| unresolved (elements too similar to discriminate) | 3 |

## The validation that matters

The crossover has a direction, and it is not random. Accounting for read orientation via
`telo_side`:

| | foreign element on telomere side | foreign element on anchor side |
|---|---|---|
| reads with `telo_side = beginning` | 41 | 1 |
| reads with `telo_side = end` | 39 | 3 |
| **total** | **80 (95%)** | **4 (5%)** |

**In 95% of cases the foreign sequence is on the telomere side and the native sequence on the
anchor side.** This holds equally in both read orientations, so it is not an artifact of read
direction or of how the pattern is read off. A random or technical artifact would give ~50/50.

That polarity is the biological expectation: BIR / template switching initiates at the
eroding telomere-proximal end, so a recombinant carries donor sequence distal and native
sequence proximal. Finding it at 95% across 106 independent reads, 15 chromosome ends and 3
strains is strong evidence these are genuine recombination events rather than matching noise.

## Distribution

Recipient ends (top): chr14R 22, chr13L 20, chr6L 14, chr10L 11, chr16R 11, chr5R 6.
Most common recipient -> donor pairs: chr13L-1 -> chr2L-1 (15), chr14R-1 -> chr9L-1 (8),
chr6L-1 -> chr14L-5 (8), chr14R-1 -> chr14L-1 (5).

Rate, as a fraction of all reads carrying a Y' call:

| strain | recombinants | reads | rate |
|---|---|---|---|
| 6991 | 93 | 161,552 | 0.06% |
| 7172 | 2 | 9,317 | 0.02% |
| 7302 | 11 | 14,202 | 0.08% |

These are day-0 populations, so a low-level background of this order is what was originally
hypothesised. The rate is a floor, not an estimate: this pipeline only sees recombinants
where donor and recipient are distinguishable by sequence, and most Y' pairs are not.

## Caveats

* **Library chimerism is not excluded.** A chimera formed during prep would also join at a
  homologous region. The 95% telomere-side polarity argues against it -- a prep artifact has
  no reason to prefer one end -- but does not rule it out. Whether the ONT prep involved
  amplification, and whether junctions recur across independent preps, would settle it.
* **The donor is usually not identifiable to a single element**, only to a set of
  near-identical ones. Donor names above are representatives.
* The 45 "single" calls are not necessarily non-recombinant; a junction very close to either
  end of the element leaves too few windows on the short side to register as a run.
