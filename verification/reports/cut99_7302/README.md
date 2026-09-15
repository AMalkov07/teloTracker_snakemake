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

## Per-read evidence: how each half of the Y' matches

`verification/flank_identities.py` locates the junction from the window scan, splits the Y' there,
and aligns each half to BOTH references. A genuine mid-Y' switch must show the anchor-side half
favouring the expected element AND the telomere-side half favouring the donor. Halves are named by
biological side using `telo_side`, so reads sequenced in either direction are comparable.
Delta is (expected - donor) on the anchor half and (donor - expected) on the telomere half; both
must be positive. "strong" = both margins >= 1.5 %.

| read | end | expected | donor | anchor half: exp / donor | Δ | telo half: exp / donor | Δ | evidence |
|---|---|---|---|---|---|---|---|---|
| SRR33298452.217955 | chr14R | chr14R-1 | chr14L-1 | 99.34% / 96.40% | +2.94 | 96.07% / 99.61% | +3.54 | **strong** |
| SRR33298452.246202 | chr14R | chr14R-1 | chr14L-1 | 99.71% / 97.36% | +2.35 | 95.48% / 99.27% | +3.79 | **strong** |
| SRR33298452.241882 | chr5R | chr5R-1 | chr14R-1 | 99.35% / 95.32% | +4.03 | 96.45% / 99.55% | +3.10 | **strong** |
| SRR33298452.220014 | chr6L | chr6L-1 | chr8R-1 | 96.92% / 94.83% | +2.09 | 95.79% / 98.16% | +2.37 | **strong** |
| SRR33298452.272227 | chr6L | chr6L-1 | chr14L-5 | 99.04% / 96.79% | +2.25 | 95.39% / 98.15% | +2.76 | **strong** |
| SRR33298452.400002 | chr6L | chr6L-1 | chr14L-5 | 99.86% / 97.62% | +2.24 | 96.72% / 99.55% | +2.83 | **strong** |
| SRR33298452.419935 | chr6L | chr6L-1 | chr14R-1 | 99.77% / 95.50% | +4.27 | 95.78% / 99.27% | +3.49 | **strong** |
| SRR33298452.485661 | chr6L | chr6L-1 | chr16L-1 | 99.48% / 97.05% | +2.43 | 97.22% / 99.82% | +2.60 | **strong** |
| SRR33298452.573441 | chr10L | chr10L-1 | chr14R-1 | 99.55% / 98.33% | +1.22 | 97.68% / 99.67% | +1.99 | weak |
| SRR33298452.400093 | chr6L | chr6L-1 | chr14L-5 | 95.08% / 92.99% | +2.09 | 93.96% / 95.30% | +1.34 | weak |
| SRR33298452.427945 | chr6L | chr6L-1 | chr14L-5 | 97.58% / 96.94% | +0.64 | 98.42% / 99.41% | +0.99 | weak |
| SRR33298452.485827 | chr5R | chr5R-1 | chr10L-1 | 98.86% / 96.85% | +2.01 | 96.28% / 92.34% | -3.94 | **FAILS** |
| SRR33298452.59861 | chr2L | chr2L-1 | chr8R-1 | — | — | — | — | no junction |
| SRR33298452.555243 | chr8R | chr8R-1 | chr2L-1 | — | — | — | — | no junction |

| evidence | reads |
|---|---|
| **strong** -- both halves clearly favour the right reference | **8** |
| weak -- correct direction, margins under 1.5 % | 3 |
| **FAILS** -- telomere half still favours the expected element | 1 |
| no junction -- whole Y' replaced, donor wins throughout | 2 |

### What the flank test changed

The window scan alone called 11 of these recombinant. Measuring both halves against both
references revises that:

* **SRR33298452.485827** (chr5R, donor chr10L-1) **fails**. Its anchor half behaves correctly
  (+2.01 for the expected element) but its telomere half *also* favours the expected element
  (-3.94). There is no half that the donor explains better, so a mid-Y' switch to chr10L-1 is not
  supported. It was one of the 11; it should not be counted.
* **SRR33298452.427945**, previously "mixed/ambiguous", now shows the correct direction on both
  halves, though weakly (+0.64 / +0.99).
* The 2 whole-element replacements have no junction to split at, as expected.

So the defensible count for 7302 day-0 is **8 strong + 3 weak mid-Y' recombination events**, 2
whole-element replacements, and 1 unexplained -- out of 7,112 reads scored.

Note the weak calls cluster with low overall identity: SRR33298452.400093 sits at 93-95 % against
both references, so its +2.09 / +1.34 margins rest on a noisy read rather than a clean signal.

## Is there homology for a crossover in every pair? Yes -- but at 96-98 %, not 99 %+

`verification/pair_homology.py` aligns each recipient/donor element pair and reports the
high-identity blocks available for strand invasion.

| recipient | donor | longest block | identity | blocks | total homologous bp | evidence class |
|---|---|---|---|---|---|---|
| chr14R-1 | chr14L-1 | 4,973 bp | 96.70% | 2 | 6,920 | strong |
| chr5R-1 | chr14R-1 | 4,962 bp | 96.41% | 2 | 6,909 | strong |
| chr6L-1 | chr8R-1 | 3,864 bp | 97.67% | 2 | 5,898 | strong |
| chr6L-1 | chr14L-5 | 3,884 bp | 97.35% | 2 | 5,923 | strong |
| chr6L-1 | chr14R-1 | 1,963 bp | 96.23% | 1 | 1,963 | strong |
| chr6L-1 | chr16L-1 | 2,047 bp | 97.36% | 3 | 5,741 | strong |
| chr10L-1 | chr14R-1 | 4,902 bp | 97.63% | 2 | 6,839 | weak |
| chr5R-1 | chr10L-1 | 4,999 bp | 96.86% | 2 | 7,038 | FAILS |
| chr2L-1 | chr8R-1 | 3,855 bp | **99.14%** | 2 | 5,889 | no junction |
| chr8R-1 | chr2L-1 | 3,855 bp | **99.14%** | 2 | 5,889 | no junction |

**No pair lacks a homologous region.** Every one has at least ~2 kb of contiguous homology,
most have 4-5 kb, and total homologous sequence runs 1.9-7.0 kb. There is always somewhere a
crossover could occur.

But the identity is **96-98 %, not 99 %+**. At a >= 98 % threshold, eight of the ten pairs
have no qualifying block at all. These are homeologous rather than identical partners:
~3 % divergence, which is well within range for homologous recombination in yeast but is not
the near-perfect identity the word "homologous area" might suggest.

### The one pair that is different, and what it may explain

`chr2L-1` / `chr8R-1` is the only pair reaching **99.14 %** -- and it is exactly the pair
behind both "whole Y' replaced, no junction" reads (SRR33298452.59861 and .555243, which are
reciprocal: chr2L->chr8R and chr8R->chr2L).

That is consistent, though not proven here: when two elements are 99.1 % identical across
3.8 kb, a crossover anywhere in that block produces almost the same sequence, so there is no
localisable junction to find. "Whole element replaced" and "junction present but unlocalisable"
are indistinguishable at this identity. These two reads should not be reported as a
mechanistically distinct class on this evidence alone.

### What could not be established

Mapping each junction back to a coordinate in the recipient element, to test whether it falls
*inside* the homology block, did not produce reliable positions -- several junctions mapped to
position ~1, which is an alignment-endpoint artifact rather than a real junction. That
containment test needs a proper pairwise alignment of read to both references (not BLAST HSP
endpoints) and is not reported here.
