# 7302: every cut99 Y' group mismatch (full detail)

| sample | read | end | expected | donor | anchor half exp/don | telo half exp/don | split anchor/telo bp | cut % from anchor | anchor margin | telo margin | evidence |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 7302 | SRR33298452.246202 | chr14R | chr14R-1 | chr14L-1 | 99.71%/97.36% | 95.48%/99.27% | 3983/3000 | 57.0% | +2.35 | +3.79 | strong |
| 7302 | SRR33298452.217955 | chr14R | chr14R-1 | chr14L-1 | 99.34%/96.40% | 96.07%/99.61% | 2775/4194 | 39.8% | +2.94 | +3.54 | strong |
| 7302 | SRR33298452.241882 | chr5R | chr5R-1 | chr14R-1 | 99.35%/95.32% | 96.45%/99.55% | 2250/4629 | 32.7% | +4.03 | +3.10 | strong |
| 7302 | SRR33298452.485661 | chr6L | chr6L-1 | chr16L-1 | 99.48%/97.05% | 97.22%/99.82% | 1770/5250 | 25.2% | +2.43 | +2.60 | strong |
| 7302 | SRR33298452.419935 | chr6L | chr6L-1 | chr14R-1 | 99.77%/95.50% | 95.78%/99.27% | 1425/5332 | 21.1% | +4.27 | +3.49 | strong |
| 7302 | SRR33298452.400002 | chr6L | chr6L-1 | chr14L-5 | 99.86%/97.62% | 96.72%/99.55% | 2394/3450 | 41.0% | +2.24 | +2.83 | strong |
| 7302 | SRR33298452.272227 | chr6L | chr6L-1 | chr14L-5 | 99.04%/96.79% | 95.39%/98.15% | 2250/3396 | 39.9% | +2.25 | +2.76 | strong |
| 7302 | SRR33298452.220014 | chr6L | chr6L-1 | chr8R-1 | 96.92%/94.83% | 95.79%/98.16% | 2789/3000 | 48.2% | +2.09 | +2.37 | strong |
| 7302 | SRR33298452.573441 | chr10L | chr10L-1 | chr14R-1 | 99.55%/98.33% | 97.68%/99.67% | 3975/2928 | 57.6% | +1.22 | +1.99 | weak |
| 7302 | SRR33298452.427945 | chr6L | chr6L-1 | chr14L-5 | 97.58%/96.94% | 98.42%/99.41% | 4710/1125 | 80.7% | +0.64 | +0.99 | weak |
| 7302 | SRR33298452.400093 | chr6L | chr6L-1 | chr14L-5 | 95.08%/92.99% | 93.96%/95.30% | 2250/3394 | 39.9% | +2.09 | +1.34 | weak |
| 7302 | SRR33298452.485827 | chr5R | chr5R-1 | chr10L-1 | 98.86%/96.85% | 96.28%/92.34% | 4875/2600 | 65.2% | +2.01 | -3.94 | FAILS |
| 7302 | SRR33298452.59861 | chr2L | chr2L-1 | chr8R-1 | — | — | — | — | — | — | no junction |
| 7302 | SRR33298452.555243 | chr8R | chr8R-1 | chr2L-1 | — | — | — | — | — | — | no junction |

Combines the per-half identity percentages with the split-point location.
"anchor/telo half exp/don" = % identity of that half to the expected element / to
the donor element. "split anchor/telo bp" and "cut % from anchor" locate the cut:
0% = right at the anchor, 100% = right at the telomere; quantised to the 150 bp
window-scan step, so approximate. Blank split/cut columns = no located junction.
