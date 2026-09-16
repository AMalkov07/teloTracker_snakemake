# Which ends get mistaken, and for what — normalised

From the cut99 bundle: 548 rows, 160 excluding the reference defect, **88 supported**
(strong + weak). Raw counts are misleading here because ends differ enormously in both read
depth (647 to 8,286 scored reads) and element count (1 to 7), so both sides are normalised.

## Recipients — events per 10,000 scored reads at that end

| end | supported | scored reads | per 10k |
|---|---|---|---|
| chr14R | 18 | 6,582 | **27.3** |
| chr6L | 15 | 5,556 | **27.0** |
| chr10L | 10 | 4,118 | **24.3** |
| chr13L | 13 | 6,469 | **20.1** |
| chr16R | 9 | 8,100 | 11.1 |
| chr8R | 2 | 1,928 | 10.4 |
| chr5R | 4 | 7,654 | 5.2 |
| chr12L | 4 | 7,938 | 5.0 |
| chr16L | 2 | 4,090 | 4.9 |
| chr5L | 4 | 8,223 | 4.9 |
| chr7R | 2 | 5,225 | 3.8 |
| chr2L | 2 | 6,190 | 3.2 |
| chr8L | 2 | 8,286 | 2.4 |
| chr15R | 1 | 4,205 | 2.4 |
| **chr9L** | **0** | 3,744 | **0** |
| **chr14L** | **0** | 1,264 | **0** |
| **chr4R** | **0** | 957 | **0** |
| **chr12R** | **0** | 647 | **0** |

A clear top tier — chr14R, chr6L, chr10L, chr13L at 20-27 per 10k — roughly 5-10x the
bottom of the distribution. chr13L's raw count (41 mismatches) overstates it: normalised it
sits below chr14R and chr6L.

**The four zeros are structural, not biological.** chr4R, chr12R and chr14L are the
multi-copy tandem arrays whose elements all sit in the large G8 group; a read at one of them
matching another G8 element is not a mismatch at all. chr9L is paired with chr10L in G7 for
the same reason. These ends are largely *unable* to register a mismatch at this cutoff, so
their zeros carry no information about recombination.

## Donors — events per element the end contributes

| end | as donor | elements (mean) | per element |
|---|---|---|---|
| chr2L | 12 | 1.0 | **12.0** |
| chr13L | 11 | 1.3 | **8.5** |
| chr7R | 8 | 1.0 | **8.0** |
| chr9L | 6 | 1.0 | 6.0 |
| chr14L | 25 | 4.8 | 5.2 |
| chr14R | 5 | 1.0 | 5.0 |
| chr16L | 6 | 1.2 | 5.0 |
| chr10L | 3 | 1.0 | 3.0 |
| **chr12R** | **3** | **6.9** | **0.4** |
| **chr4R** | **0** | **6.3** | **0** |

chr14L is the top donor by raw count (25) but only mid-table once its 4.8 elements are
accounted for. The striking result is the opposite one: **chr4R and chr12R contribute the most
elements of any end (6-7 each) and are essentially never donors.**

## Two patterns that explain most of this

**1. Size class is respected — 83 %.**

| | |
|---|---|
| Short → Short | 42 |
| Long → Long | 31 |
| Long → Short | 10 |
| Short → Long | 5 |

73 of 88 events stay within a size class. This is why chr4R and chr12R (all Long, all in G8)
rarely donate: the high-rate recipients chr6L, chr13L and chr16R are Short, and their
best-matching partners are the Short elements at chr14L (chr14L-3/4/5), chr2L and chr13L.
chr14L's prominence as a donor comes from carrying *both* classes — Long at copies 1-2, Short
at copies 3-5 — so it can serve either.

**2. The direction is almost entirely one-way.**

| pair | A→B | B→A |
|---|---|---|
| chr13L → chr2L | 11 | 0 |
| chr6L → chr14L | 10 | 0 |
| chr16R → chr14L | 6 | 0 |
| chr14R → chr14L | 5 | 0 |
| chr14R → chr9L | 5 | 0 |
| chr10L → chr7R | 4 | 0 |
| chr14R → chr7R | 4 | 1 |
| chr5L → chr12R | 3 | 0 |

Nine of the top twelve pairs are strictly one-directional. If these were symmetric sequence
confusions between similar elements, reciprocal counts would be expected; they are not seen.

## What this does and does not establish

The one-way direction and the size-class constraint are consistent with directed
recombination — a recipient end acquiring sequence from a donor — rather than symmetric
mis-assignment noise.

But two cautions apply. First, **detectability is not uniform across pairs**: a pair is only
visible when the two elements fall in different cut99 groups, and group membership is not
symmetric in its consequences, so some apparent one-wayness may be structural. Second, 59 % of
the supported reads belong to haplotypes recurring across independent preps, so much of this
pattern describes **standing variation already in the cultures**, not events occurring at
day 0 — the donor/recipient asymmetry may be a historical record rather than an active process.
