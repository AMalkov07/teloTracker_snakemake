# 7302 day-0: Y' copies that miss their reference group at the 99 % cutoff

Samples: 7302_day0_with_selection

| evidence | reads | meaning |
|---|---|---|
| **strong** | 8 | mid-Y' recombination, both halves clearly favour the right reference |
| **weak** | 3 | mid-Y' recombination, correct direction but margins < 1.5 % |
| **FAILS** | 1 | no half is better explained by the donor -- unexplained |
| **no junction** | 2 | donor wins across the element, or signal too weak to split |

| sample | read | end | expected | donor | anchor half exp/don | Δ | telo half exp/don | Δ | pair homology | evidence |
|---|---|---|---|---|---|---|---|---|---|---|
| 7302_day0_with_selection | SRR33298452.246202 | chr14R | chr14R-1 | chr14L-1 | 99.71%/97.36% | +2.35 | 95.48%/99.27% | +3.79 | 4973bp @96.70% | strong |
| 7302_day0_with_selection | SRR33298452.217955 | chr14R | chr14R-1 | chr14L-1 | 99.34%/96.40% | +2.94 | 96.07%/99.61% | +3.54 | 4973bp @96.70% | strong |
| 7302_day0_with_selection | SRR33298452.241882 | chr5R | chr5R-1 | chr14R-1 | 99.35%/95.32% | +4.03 | 96.45%/99.55% | +3.10 | 4962bp @96.41% | strong |
| 7302_day0_with_selection | SRR33298452.485661 | chr6L | chr6L-1 | chr16L-1 | 99.48%/97.05% | +2.43 | 97.22%/99.82% | +2.60 | 2047bp @97.36% | strong |
| 7302_day0_with_selection | SRR33298452.419935 | chr6L | chr6L-1 | chr14R-1 | 99.77%/95.50% | +4.27 | 95.78%/99.27% | +3.49 | 1963bp @96.23% | strong |
| 7302_day0_with_selection | SRR33298452.400002 | chr6L | chr6L-1 | chr14L-5 | 99.86%/97.62% | +2.24 | 96.72%/99.55% | +2.83 | 3884bp @97.35% | strong |
| 7302_day0_with_selection | SRR33298452.272227 | chr6L | chr6L-1 | chr14L-5 | 99.04%/96.79% | +2.25 | 95.39%/98.15% | +2.76 | 3884bp @97.35% | strong |
| 7302_day0_with_selection | SRR33298452.220014 | chr6L | chr6L-1 | chr8R-1 | 96.92%/94.83% | +2.09 | 95.79%/98.16% | +2.37 | 3864bp @97.67% | strong |
| 7302_day0_with_selection | SRR33298452.573441 | chr10L | chr10L-1 | chr14R-1 | 99.55%/98.33% | +1.22 | 97.68%/99.67% | +1.99 | 4902bp @97.63% | weak |
| 7302_day0_with_selection | SRR33298452.427945 | chr6L | chr6L-1 | chr14L-5 | 97.58%/96.94% | +0.64 | 98.42%/99.41% | +0.99 | 3884bp @97.35% | weak |
| 7302_day0_with_selection | SRR33298452.400093 | chr6L | chr6L-1 | chr14L-5 | 95.08%/92.99% | +2.09 | 93.96%/95.30% | +1.34 | 3884bp @97.35% | weak |
| 7302_day0_with_selection | SRR33298452.485827 | chr5R | chr5R-1 | chr10L-1 | 98.86%/96.85% | +2.01 | 96.28%/92.34% | -3.94 | 4999bp @96.86% | FAILS |
| 7302_day0_with_selection | SRR33298452.59861 | chr2L | chr2L-1 | chr8R-1 | — | — | — | — | 3855bp @99.14% | no junction |
| 7302_day0_with_selection | SRR33298452.555243 | chr8R | chr8R-1 | chr2L-1 | — | — | — | — | 3855bp @99.14% | no junction |
