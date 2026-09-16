# chr16R->G2, chr10L->G8, chr2L->G3: recombination confirmation and homology ranking

Same method as the chr13L->G1 analysis: each test read's whole Y' BLASTed against its own
expected reference and against the donor group (best member), with a control of 6
correctly-assigned reads for that end run through the identical test to establish the native
baseline (own reference wins comfortably, donor gets only a partial hit).

## Part 1: is it real recombination?

### chr16R -> G2 {chr13L-1, chr14L-3/4/5}: 7/7 confirmed

| read | window-scan class | vs own (chr16R-1) | vs G2 best | margin | verdict |
|---|---|---|---|---|---|
| SRR33298373.72733 | weak | 99.3% / 3,578bp / bs=6,469 | 99.0% / 5,199bp (chr14L-5) / bs=9,286 | +2,817 | **recombinant** |
| SRR33298384.248800 | strong | 99.5% / 3,583bp / bs=6,508 | 98.9% / 5,203bp (chr14L-5) / bs=9,280 | +2,772 | **recombinant** |
| SRR33298384.538094 | strong | 98.1% / 3,599bp / bs=6,233 | 97.6% / 5,230bp (chr13L-1) / bs=8,911 | +2,678 | **recombinant** |
| SRR33298377.534605 | strong | 98.2% / 3,591bp / bs=6,252 | 97.4% / 5,209bp (chr13L-1) / bs=8,820 | +2,568 | **recombinant** |
| SRR33298434.378016 | FAILS | 96.2% / 3,612bp / bs=5,847 | 95.8% / 5,252bp (chr14L-5) / bs=8,386 | +2,539 | **recombinant** |
| SRR33298377.644709 | FAILS | 98.1% / 3,596bp / bs=6,237 | 97.3% / 5,214bp (chr13L-1) / bs=8,776 | +2,539 | **recombinant** |
| SRR33298434.64148 | FAILS | 94.4% / 3,611bp / bs=5,459 | 92.5% / 5,265bp (chr14L-5) / bs=7,352 | +1,893 | **recombinant** |

**7 of 7**, margins +1,893 to +2,817. Notably the 3 reads the sliding-window scan called
"FAILS" all confirm here -- the same pattern already seen for one chr13L read
(SRR33298384.185796): the window scan can miss a real event when the junction is close to
one end of the read, but the whole-read bitscore test still catches it.

### chr10L -> G8 (18-member array group): 6/6 confirmed

| read | window-scan class | vs own (chr10L-1) | vs G8 best | margin | verdict |
|---|---|---|---|---|---|
| SRR33298384.498901 | strong | 98.2% / 5,003bp / bs=8,698 | 98.7% / 6,614bp (chr7R-1) / bs=11,725 | +3,027 | **recombinant** |
| SRR33298461.57209 | strong | 97.5% / 5,010bp / bs=8,510 | 98.0% / 6,623bp (chr7R-1) / bs=11,437 | +2,927 | **recombinant** |
| SRR33298373.491903 | FAILS | 97.6% / 5,010bp / bs=8,540 | 98.0% / 6,619bp (chr16L-1) / bs=11,452 | +2,912 | **recombinant** |
| SRR33298373.545529 | strong | 98.0% / 5,002bp / bs=8,667 | 98.3% / 6,611bp (chr16L-1) / bs=11,568 | +2,901 | **recombinant** |
| SRR33298384.151130 | weak | 94.1% / 5,043bp / bs=7,542 | 95.2% / 6,664bp (chr7R-1) / bs=10,405 | +2,863 | **recombinant** |
| SRR33298373.122658 | weak | 96.0% / 5,035bp / bs=8,117 | 95.8% / 6,665bp (chr4R-7) / bs=10,661 | +2,544 | **recombinant** |

**6 of 6**, margins +2,544 to +3,027 -- the tightest, most consistent set of the four ends
examined so far. Again the one "FAILS" read (491903) confirms under the whole-read test.

### chr2L -> G3 {chr8R-1}: only 2/5 confirmed -- the weakest case of the four

| read | window-scan class | vs own (chr2L-1) | vs G3 (chr8R-1) | margin | verdict |
|---|---|---|---|---|---|
| SRR33298373.284052 | no junction | 98.8% / 3,861bp / bs=6,863 | 99.6% / 5,473bp / bs=9,987 | +3,124 | **recombinant** |
| SRR33298434.94939 | strong | 95.7% / 2,042bp / bs=3,262 | 97.2% / 3,664bp / bs=6,150 | +2,888 | **recombinant** |
| SRR33298384.313055 | no junction | 98.9% / 2,270bp / bs=4,052 | 98.9% / 2,270bp / bs=4,052 | **0** | not supported |
| SRR33298384.467718 | no junction | 97.8% / 3,175bp / bs=5,448 | 97.8% / 3,172bp / bs=5,443 | -5 | not supported |
| SRR33298377.541858 | FAILS | 97.6% / 3,946bp / bs=6,728 | 97.5% / 3,868bp / bs=6,584 | -144 | not supported |

Only 2 of 5 show a real contest. The other 3 have **near-identical bitscores against both
references** (margins of 0, -5, -144, against short alignments of only 2,270-3,946 bp) --
own and donor are not competing for different parts of the read, they are returning
essentially the same weak, partial hit. That is not the recombination signature; it looks
like a poor-quality or low-information read at that locus, not a confirmed event either way.

**chr2L -> G3 should be treated as 2 confirmed events out of 5 calls, not 4/5** as the raw
mismatch count implied.

## Part 2: why these donor groups? Homology ranking

Each recipient's own reference BLASTed against every other element in the library, ranked by
bitscore, with elements in the recipient's own cut99 group (invisible to detection, as
established for chr13L) marked.

### chr16R-1: G2 is the #1-4 most homologous group

| rank | element | identity | bitscore |
|---|---|---|---|
| 1-3 | **chr14L-3/4/5 (G2, observed donor)** | 98.55% | 6,324 |
| (tied) | **chr13L-1 (G2, observed donor)** | 98.61% | 6,333 |
| 5 | chr6L-1 | 97.59% | 6,307 |
| 6 | chr12L-1 | 98.41% | 6,298 |
| 7 | chr2L-1 | 97.82% | 6,176 |
| 8 | chr8R-1 | 98.27% | 6,074 |

chr16R-1 has no partner in its own cut99 group (it is a singleton, G6), so nothing is hidden
from this ranking. **G2 occupies the top 4 positions outright.** This is the opposite of the
chr13L case: here the dominant donor is genuinely the most homologous available option.

### chr10L-1: G8 is the #1 detectable group by a wide margin

| rank | element | bitscore | in G8? |
|---|---|---|---|
| 1 | chr9L-1 (own group, G7 -- invisible) | 12,597 | no (own group) |
| 2-18 | **chr7R-1, chr14L-1, chr16L-1, chr4R-1..7, chr15R-1, chr14L-2, chr12R-2..7 (all G8)** | 8,399-8,425 | **yes, all 17** |
| 19 | chr14R-1 (not in G8 -- closest outside miss) | 8,388 | no |
| 20 | chr5R-1 | 8,301 | no |
| -- | (next tier, chr5L-1 and below) | <=5,051 | -- |

After excluding chr10L-1's own invisible group-mate, **every single one of G8's 17
detectable members occupies the top of the ranking**, edging out the nearest non-member
(chr14R-1) by only 11-37 bitscore points before the ranking drops sharply (8,301 -> 5,051, a
3,000+ point cliff). G8 is not just the most homologous choice -- it is essentially the only
plausible one; nothing else comes close.

### chr2L-1: chr8R-1 (G3) is the #1 detectable candidate

| rank | element | bitscore | note |
|---|---|---|---|
| 1 | chr6L-1 (own group, G1 -- invisible) | 10,717 | own group |
| 2 | **chr8R-1 (G3, observed donor)** | **6,937** | |
| 3-5 | chr14L-3/4/5 | 6,541 | |
| 6 | chr13L-1 | 6,514 | |

Excluding the invisible own-group partner, chr8R-1 is the single most homologous detectable
element to chr2L-1. So even though only 2 of the 5 chr2L->G3 calls hold up as confirmed
recombination, the ones that do are consistent with homology-driven pairing -- unlike
chr13L, where the dominant donor was one of the *least* homologous options.

## The contrast with chr13L

| recipient | donor rank among detectable candidates | recombination confirmed |
|---|---|---|
| chr16R -> G2 | **#1 (tied top 4)** | 7/7 |
| chr10L -> G8 | **#1 (by a wide margin)** | 6/6 |
| chr2L -> G3 | **#1** | 2/5 |
| chr13L -> G1 | **#8-9 of ~12** (one of the least homologous) | 13/15 |

Three of four recipients recombine with their single most sequence-homologous available
partner -- exactly what a homology-driven strand-invasion model predicts, and a clean,
sufficient answer for why those donor groups dominate. **chr13L is the outlier**, and the
earlier finding stands: its dominant donor (G1) is not explained by homology at all, and the
best current explanation for chr13L specifically is standing variation from a small number of
historical events rather than an ongoing homology-driven preference.

This also sharpens what "why does recombination pick the same group" means: for most
recipients examined so far, the answer is straightforwardly "because that group is the best
available sequence match." chr13L needed a different explanation because that answer is
false for it specifically -- it is not true of every recipient in this dataset.
