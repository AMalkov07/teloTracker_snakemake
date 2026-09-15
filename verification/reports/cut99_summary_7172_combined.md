# 7172: every cut99 Y' group mismatch

| sample | read | end | expected | donor | split (anchor/telo bp) | cut at % from anchor | anchor margin | telo margin | evidence |
|---|---|---|---|---|---|---|---|---|---|
| 7172 | SRR33298432.177252 | chr16R | chr16R-1 | chr14L-3 | 2219/3375 | 39.7% | +1.61 | +2.97 | strong |
| 7172 | SRR33298432.315476 | chr16R | chr16R-1 | chr14L-3 | 5100/401 | 92.7% | +1.37 | +1.03 | weak |
| 7172 | SRR33298432.328531 | chr5R | chr5R-1 | chr13L-1 | 1154/4650 | 19.9% | +7.01 | +0.75 | weak |
| 7172 | SRR33298432.77931 | chr6L | chr6L-1 | chr14L-3 | 3675/1889 | 66.0% | +1.68 | -0.06 | FAILS |

`cut at % from anchor`: where the anchor-side half ends and the telomere-side
half begins, as a percentage of the read's total Y' span, measured from the
anchor. 0% = the cut sits right at the anchor (an almost entirely foreign Y');
100% = right at the telomere (an almost entirely native Y'). Located from the
sliding-window scan (`scan_recombinant_junctions.py`), quantised to its 150 bp
step, so treat it as approximate rather than base-pair precise. Blank for rows
with no located junction (`no junction`, `reference defect`).
