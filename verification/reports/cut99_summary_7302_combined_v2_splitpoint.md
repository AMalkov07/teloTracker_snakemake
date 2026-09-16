# 7302: every cut99 Y' group mismatch (split-point view)

| sample | read | end | expected | donor | split (anchor/telo bp) | cut at % from anchor | anchor margin | telo margin | evidence |
|---|---|---|---|---|---|---|---|---|---|
| 7302 | SRR33298452.246202 | chr14R | chr14R-1 | chr14L-1 | 3983/3000 | 57.0% | +2.35 | +3.79 | strong |
| 7302 | SRR33298452.217955 | chr14R | chr14R-1 | chr14L-1 | 2775/4194 | 39.8% | +2.94 | +3.54 | strong |
| 7302 | SRR33298452.241882 | chr5R | chr5R-1 | chr14R-1 | 2250/4629 | 32.7% | +4.03 | +3.10 | strong |
| 7302 | SRR33298452.485661 | chr6L | chr6L-1 | chr16L-1 | 1770/5250 | 25.2% | +2.43 | +2.60 | strong |
| 7302 | SRR33298452.419935 | chr6L | chr6L-1 | chr14R-1 | 1425/5332 | 21.1% | +4.27 | +3.49 | strong |
| 7302 | SRR33298452.400002 | chr6L | chr6L-1 | chr14L-5 | 2394/3450 | 41.0% | +2.24 | +2.83 | strong |
| 7302 | SRR33298452.272227 | chr6L | chr6L-1 | chr14L-5 | 2250/3396 | 39.9% | +2.25 | +2.76 | strong |
| 7302 | SRR33298452.220014 | chr6L | chr6L-1 | chr8R-1 | 2789/3000 | 48.2% | +2.09 | +2.37 | strong |
| 7302 | SRR33298452.573441 | chr10L | chr10L-1 | chr14R-1 | 3975/2928 | 57.6% | +1.22 | +1.99 | weak |
| 7302 | SRR33298452.427945 | chr6L | chr6L-1 | chr14L-5 | 4710/1125 | 80.7% | +0.64 | +0.99 | weak |
| 7302 | SRR33298452.400093 | chr6L | chr6L-1 | chr14L-5 | 2250/3394 | 39.9% | +2.09 | +1.34 | weak |
| 7302 | SRR33298452.485827 | chr5R | chr5R-1 | chr10L-1 | 4875/2600 | 65.2% | +2.01 | -3.94 | FAILS |
| 7302 | SRR33298452.59861 | chr2L | chr2L-1 | chr8R-1 | — | — | — | — | no junction |
| 7302 | SRR33298452.555243 | chr8R | chr8R-1 | chr2L-1 | — | — | — | — | no junction |

`cut at % from anchor`: where the anchor-side half ends and the telomere-side
half begins, as a percentage of the read's total Y' span, measured from the
anchor. 0% = the cut sits right at the anchor; 100% = right at the telomere.
Located from the sliding-window scan, quantised to its 150 bp step -- approximate,
not base-pair precise. Blank for rows with no located junction.
