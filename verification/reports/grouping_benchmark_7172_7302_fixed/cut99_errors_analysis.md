# Reads the `cut99` scheme mislabels — 7172 and 7302 day-0 (post-fix)

18 misread reads across 11,583 scored non-recombinant reads (0.16%). Files:
`cut99_errors_reads.tsv` (all 18 reads, ONT UUIDs included), `cut99_errors_confusion.tsv`
(enriched with element lengths, size class and pairwise similarity).

## Every confusion, ranked by similarity to the true element

| strain | true | matched | copies | true len | matched len | class | similarity | verdict |
|---|---|---|---|---|---|---|---|---|
| 7302 | chr6L-1 | chr14L-5 | 4 | 5975 | 5485 | Short/Short | 63.03% | same class |
| 7172 | chr16R-1 | chr14L-3 | 2 | 5302 | 5485 | Short/Short | 64.46% | same class |
| 7302 | chr14R-1 | chr14L-1 | 2 | 6530 | 6654 | Long/Long | 72.38% | same class |
| 7172 | chr5R-1 | chr13L-1 | 1 | 6693 | 5483 | Long/Short | 54.34% | CROSS class |
| 7172 | chr6L-1 | chr14L-3 | 1 | 5975 | 5485 | Short/Short | 63.03% | same class |
| 7302 | chr10L-1 | chr14R-1 | 1 | 6869 | 6530 | Long/Long | 69.67% | same class |
| 7302 | chr2L-1 | chr8R-1 | 1 | 5975 | 5469 | Short/Short | 63.98% | same class |
| 7302 | chr5R-1 | chr14R-1 | 1 | 6693 | 6530 | Long/Long | 71.14% | same class |
| 7302 | chr5R-1 | chr10L-1 | 1 | 6693 | 6869 | Long/Long | 70.49% | same class |
| 7302 | chr6L-1 | chr8R-1 | 1 | 5975 | 5469 | Short/Short | 63.10% | same class |
| 7302 | chr6L-1 | chr14R-1 | 1 | 5975 | 6530 | Short/Long | 28.97% | CROSS class |
| 7302 | chr6L-1 | chr16L-1 | 1 | 5975 | 6654 | Short/Long | 30.04% | CROSS class |
| 7302 | chr8R-1 | chr2L-1 | 1 | 5469 | 5975 | Short/Short | 63.98% | same class |

**None of the 18 involve a near-identical pair** (best is 72%; most are 55-64%). Compare
chr16L's own defect, where the confusion was between elements 1-3 bp apart at ≥99.94%
similarity. That distinction is the basis for the answer to "could these be real day-0
recombination at low frequency":

**Argument against real recombination.** A genuine template-switch copies the donor's
actual sequence, so a real recombinant Y' should match its donor at high identity — the
same logic the pipeline itself uses to call recombination on later time points. A "match"
at 55-72% identity over 5-6 kb is not evidence the read's copy came from that donor; it is
the least-bad option among references that are all a poor fit. If any of the 18 came from
a real day-0 template switch, the true donor would be expected to show up in this table at
>=95% identity somewhere, and none do.

## Robustness across schemes (per-read)

| scheme | misread reads | shared with cut99 |
|---|---|---|
| cut99 | 18 | — |
| cut97 | 18 | 100% |
| silhouette | 18 | 100% |
| condensed | 29 | 100% (cut99's 18 are a subset) |
| curated_family | 18 | 89% |
| **wrong under all five** | | **16 of 18** |

Nearly every misread read is wrong under every scheme simultaneously — this is a per-read
matching failure common to all resolutions, not something a particular grouping choice
produces or could fix. Read length is not a factor either: misread reads have a median of
23,658 bp against 15,932 bp for the whole scored population — if anything longer, not
shorter/lower-quality.

## chr6L is the one real trend: 8 of 18 (44%)

chr6L accounts for far more than its share. Its only strong reference partner is
`chr2L-1` (99.0% identical, already grouped together in every scheme) — beyond that its
next-best matches are all ~60-64%, a cliff rather than a gradient. And 6 of its 8
misreads have a striking pattern: the actual RepeatMasker hit on the read spans only
~5,360-5,490 bp, systematically ~500-600 bp shorter than `chr6L-1`'s own reference length
(5,975 bp) and close instead to the Short-class chr14L reference length (5,485 bp) it gets
mismatched to:

| read | matched as | hit span | chr6L-1 ref | shortfall |
|---|---|---|---|---|
| SRR33298452.220014 | chr8R-1 | 5,478 bp | 5,975 bp | 497 bp |
| SRR33298452.272227 | chr14L-5 | 5,455 bp | 5,975 bp | 520 bp |
| SRR33298452.400002 | chr14L-5 | 5,468 bp | 5,975 bp | 507 bp |
| SRR33298452.400093 | chr14L-5 | 5,428 bp | 5,975 bp | 547 bp |
| SRR33298452.427945 | chr14L-5 | 5,492 bp | 5,975 bp | 483 bp |
| SRR33298432.77931 | chr14L-3 | 5,361 bp | 5,975 bp | 614 bp |

This is worth a closer look but is NOT yet diagnosed here: it could mean these particular
reads carry a genuinely shorter Y' copy at chr6L (~500 bp partial/truncated), a
sequencing/basecalling gap breaking the alignment mid-copy, or (less likely, given
chr6L-1/chr2L-1 already agree with each other at 99%) a boundary issue in the chr6L-1
reference itself analogous to chr16L's. Any of these would produce exactly the low-identity,
short-span "mismatch" seen here without it being recombination. The two chr6L misreads NOT
in this pattern (spans 6,514 and 6,633 bp, matched to Long-class elements) don't fit it and
remain unexplained.

## Bottom line

These are very unlikely to be real day-0 recombination events. The identity argument is the
strongest one: recombination in this pipeline is called by high-identity donor matches, and
none of these 18 mismatches have one. They look like a residual per-read Y' matching floor —
the same kind of thing found in 6991's cut97 analysis (also zero near-identical confusions,
also scheme-independent) — with chr6L as a specific, reproducible hotspot worth investigating
further if finer resolution at that locus matters.
