# 7302 day-0 alone: which Y' are mismatched, and what for

14 rows, 11 supported (8 strong + 3 weak), 2 no-junction, 1 FAILS. 7,112 reads scored across
18 ends. **At n=11 almost nothing here is statistically robust** — one finding survives that
bar, the rest are observations.

## The one solid finding: chr6L

| recipient | class | supported | scored reads | per 10k |
|---|---|---|---|---|
| **chr6L** | Short | **7** | 575 | **121.7** |
| chr14R | Long | 2 | 535 | 37.4 |
| chr10L | Long | 1 | 359 | 27.9 |
| chr5R | Long | 1 | 580 | 17.2 |
| chr16R, chr5L, chr8L, chr12L, … | | 0 | ~500-600 each | 0 |

chr6L holds **8.1 % of the reads but 7 of 11 events**. Against a null where events follow read
depth, P(≥7 of 11 at chr6L) = **5.6 × 10⁻⁶**. Depth is near-uniform across these ends
(~350-600 reads each), so this is not a coverage artefact.

**chr6L draws from four different donors**, so this is not one clonal haplotype inflating a
single end:

| donor | events |
|---|---|
| chr14L-5 | 4 |
| chr16L-1 | 1 |
| chr14R-1 | 1 |
| chr8R-1 | 1 |

## Observations that do NOT reach significance at n=11

**chr14L is the most frequent donor end** (6 of 14 rows: chr14L-5 ×4, chr14L-1 ×2). It is the
only 7302 end supplying donors in both size classes — Long at copy 1, Short at copy 5 — which
would explain versatility, but 6 events cannot establish it.

**Size class is respected in 9 of 11** (Long→Long 4, Short→Short 5, Short→Long 2). Consistent
with the 83 % seen across all strains, but on its own 9/11 is weak.

**No reciprocity.** chr6L is a recipient 7 times and never a donor; chr14L is a donor 6 times
and never a recipient. Suggestive of direction, but with most pairs appearing once or twice the
absence of a reverse event is expected anyway.

**The two no-junction rows are a reciprocal pair** — chr2L→chr8R and chr8R→chr2L, one read
each. These are the only reciprocal pair in the strain, and chr2L-1/chr8R-1 is also the only
7302 pair at 99.14 % identity, where a junction cannot be localised. Their being reciprocal is
likely a consequence of that near-identity rather than of two separate events.

## Full table

| read | recipient | expected | donor | groups | evidence |
|---|---|---|---|---|---|
| SRR33298452.246202 | chr14R | chr14R-1 | chr14L-1 | G10→G8 | strong |
| SRR33298452.217955 | chr14R | chr14R-1 | chr14L-1 | G10→G8 | strong |
| SRR33298452.241882 | chr5R | chr5R-1 | chr14R-1 | G9→G10 | strong |
| SRR33298452.485661 | chr6L | chr6L-1 | chr16L-1 | G1→G8 | strong |
| SRR33298452.419935 | chr6L | chr6L-1 | chr14R-1 | G1→G10 | strong |
| SRR33298452.400002 | chr6L | chr6L-1 | chr14L-5 | G1→G2 | strong |
| SRR33298452.272227 | chr6L | chr6L-1 | chr14L-5 | G1→G2 | strong |
| SRR33298452.220014 | chr6L | chr6L-1 | chr8R-1 | G1→G3 | strong |
| SRR33298452.573441 | chr10L | chr10L-1 | chr14R-1 | G7→G10 | weak |
| SRR33298452.427945 | chr6L | chr6L-1 | chr14L-5 | G1→G2 | weak |
| SRR33298452.400093 | chr6L | chr6L-1 | chr14L-5 | G1→G2 | weak |
| SRR33298452.485827 | chr5R | chr5R-1 | chr10L-1 | G9→G7 | FAILS |
| SRR33298452.59861 | chr2L | chr2L-1 | chr8R-1 | G1→G3 | no junction |
| SRR33298452.555243 | chr8R | chr8R-1 | chr2L-1 | G3→G1 | no junction |

## Caveat carried from the cross-strain analysis

Four ends cannot register a mismatch at this cutoff because their elements share a group with
their likely partners — chr4R, chr12R and chr14L (all in G8) and chr9L (with chr10L in G7).
Their absence from the recipient column is structural. chr14L appearing only as a donor and
never as a recipient is an instance of this, not evidence of directionality.
