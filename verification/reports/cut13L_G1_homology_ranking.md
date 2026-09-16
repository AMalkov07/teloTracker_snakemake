# Is G1 (chr2L-1/chr6L-1) the most homologous group to chr13L-1? No -- it isn't even close.

## Full ranking, chr13L-1 vs every other element in the 6991 library

| rank | element | identity | aligned length | class |
|---|---|---|---|---|
| own | chr13L-1 | 100.0% | 5,483 (full) | -- |
| 1-3 | **chr14L-3/4/5** (chr13L-1's own G2 group) | **99.60%** | 5,482 (full) | Short |
| 4 | chr8R-1 | 98.56% | 5,495 (full) | Short |
| 5 | chr12L-1 | 99.17% | 5,188 (95% of length) | Short |
| 6 | chr8L-1 | 98.14% | 5,059 (92% of length) | Short |
| 7 | chr16R-1 | 98.61% | 3,588 (65% of length) | Short |
| **8** | **chr6L-1 (G1)** | **97.27%** | 3,884 (71% of length) | Short |
| **9** | **chr2L-1 (G1)** | **97.12%** | 3,883 (71% of length) | Short |
| lower | chr5L-1, chr10L-1, chr9L-1, chr14R-1 | 95.6-96.4% | ~2,000 (37% of length) | Long |

**G1's two members rank 8th and 9th out of the ~12 groups in the genome** -- lower identity
and lower coverage than chr8R-1, chr12L-1, chr8L-1, and chr16R-1, every one of which is a
*different*, independently detectable group. If recombination simply favoured whichever
sequence was most alike, chr13L should be dominated by chr14L-3/4/5, then chr8R-1/chr12L-1/
chr8L-1, with chr2L-1/chr6L-1 barely registering. The opposite happened: G1 accounts for 15 of
18 events; chr8R-1, chr12L-1 and chr8L-1 together account for 1.

## The catch: the single most homologous group is invisible to this whole method

chr13L-1's actual closest relatives, chr14L-3/4/5 at 99.60% (near-perfect, full-length), sit
**inside chr13L-1's own cut99 group (G2)**. A read whose Y' converted to chr14L-3/4/5-type
sequence would still match group G2 -- indistinguishable from a read carrying native chr13L-1 --
so it is never counted as a mismatch at all. **This detection method structurally cannot see
recombination with a read's own closest relatives**, only with sequences dissimilar enough to
land in a different group. That is a real limitation of the whole cut99 mismatch-counting
approach, not specific to chr13L: whatever recombines most easily with each recipient's
single nearest neighbour is the recombination this method is least able to detect.

So the honest comparison is chr2L-1/chr6L-1 against the other *detectable* (different-group)
candidates -- chr8R-1, chr12L-1, chr8L-1, chr16R-1 -- and G1 still loses that comparison on
identity while winning overwhelmingly on event count.

## What this means

Sequence homology of the Y' element itself does not explain why chr13L recombines
preferentially with chr2L/chr6L. Something else is driving that specific pairing --
candidates, none confirmed here:

* **Flanking sequence, not the Y' itself.** The X element, spacer, or anchor region adjacent
  to the Y' could be what actually pairs during strand invasion/synapsis, with the Y' just
  carried along; those regions were not compared here.
* **A historical event now fixed as standing variation.** The earlier investigation showed
  most of these events (9 of 15, the whole-element-replacement signature) recur identically
  across three independent sequencing runs -- consistent with one or a few ancestral exchange
  events between the chr13L and chr2L/chr6L loci that became fixed in the population, rather
  than an ongoing preference being exercised fresh each generation. A fixed historical event
  does not need to have used the most-homologous available partner at the time it occurred;
  it only needs to have happened once.
* **Spatial/nuclear proximity** independent of sequence -- some subtelomeres cluster more
  than others in the nucleus, which would favour a specific pairing regardless of Y' identity.

This was not resolved here. What is resolved is that "G1 is simply the closest sequence
match" is not the explanation -- it is measurably one of the worse matches among the
detectable candidates.
