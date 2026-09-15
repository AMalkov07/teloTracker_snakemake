# 7172: every cut99 Y' group mismatch (per-half identity view)

| sample | read | end | expected | donor | anchor half exp/don | anchor margin | telo half exp/don | telo margin | pair homology | evidence |
|---|---|---|---|---|---|---|---|---|---|---|
| 7172 | SRR33298432.177252 | chr16R | chr16R-1 | chr14L-3 | 99.95%/98.34% | +1.61 | 96.91%/99.88% | +2.97 | 3588bp @98.55% | strong |
| 7172 | SRR33298432.315476 | chr16R | chr16R-1 | chr14L-3 | 99.97%/98.60% | +1.37 | 98.97%/100.00% | +1.03 | 3588bp @98.55% | weak |
| 7172 | SRR33298432.328531 | chr5R | chr5R-1 | chr13L-1 | 96.25%/89.24% | +7.01 | 97.34%/98.09% | +0.75 | 3656bp @99.62% | weak |
| 7172 | SRR33298432.77931 | chr6L | chr6L-1 | chr14L-3 | 99.61%/97.93% | +1.68 | 96.66%/96.60% | -0.06 | 3884bp @97.35% | FAILS |

Per-half identity view (no split-point columns). "anchor half exp/don" = % identity of the anchor-side half to the expected element / to the donor element; "telo half exp/don" likewise for the telomere-side half. Positive anchor margin = anchor half favours the expected element; positive telo margin = telo half favours the donor.
