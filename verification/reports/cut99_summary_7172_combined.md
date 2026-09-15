# 7172: every cut99 Y' group mismatch (full detail)

| sample | read | end | expected | donor | anchor half exp/don | telo half exp/don | split anchor/telo bp | cut % from anchor | anchor margin | telo margin | evidence |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 7172 | SRR33298432.177252 | chr16R | chr16R-1 | chr14L-3 | 99.95%/98.34% | 96.91%/99.88% | 2219/3375 | 39.7% | +1.61 | +2.97 | strong |
| 7172 | SRR33298432.315476 | chr16R | chr16R-1 | chr14L-3 | 99.97%/98.60% | 98.97%/100.00% | 5100/401 | 92.7% | +1.37 | +1.03 | weak |
| 7172 | SRR33298432.328531 | chr5R | chr5R-1 | chr13L-1 | 96.25%/89.24% | 97.34%/98.09% | 1154/4650 | 19.9% | +7.01 | +0.75 | weak |
| 7172 | SRR33298432.77931 | chr6L | chr6L-1 | chr14L-3 | 99.61%/97.93% | 96.66%/96.60% | 3675/1889 | 66.0% | +1.68 | -0.06 | FAILS |

Combines the per-half identity percentages with the split-point location.
"anchor/telo half exp/don" = % identity of that half to the expected element / to
the donor element. "split anchor/telo bp" and "cut % from anchor" locate the cut:
0% = right at the anchor, 100% = right at the telomere; quantised to the 150 bp
window-scan step, so approximate. Blank split/cut columns = no located junction.
