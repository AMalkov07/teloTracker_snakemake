# Defending (and correcting) the "2 ancestral molecules" claim

Direct answers first, then the full evidence and a correction to the earlier claim.

## Are the ancestral haplotypes present in all 5 matched samples?

**No.** Of the 5 samples with an identical cut99 grouping, chr13L->G1 events appear in only
**4**: `6991_day0`, `6991_day0_reference_promethion`, `6991_day0_with_selection_repeat`,
`6991_day0_with_selection_repeat2`. `6991_day0_TeloTag_with_selection` has zero -- its one
chr13L mismatch that sample did produce went to a different group (G8) entirely.

Read counts per sample: `6991_day0` 1, `reference_promethion` 4, `with_selection_repeat` 8,
`with_selection_repeat2` 2.

## Are "descendants" found within the same sample as other copies of the same haplotype?

Yes, and that within-sample recurrence is expected, not surprising -- `with_selection_repeat`
alone contributes 8 of the 15 reads. The real question is whether within-sample reads are
*more* similar to each other than cross-sample reads are (which would suggest a library
artefact -- PCR duplication or barcode crosstalk within one run -- rather than true standing
variation in the cell population). Checked directly:

| | n pairs | mean identity | range |
|---|---|---|---|
| within-sample pairs | 34 | 97.7% | 90.4-99.8% |
| cross-sample pairs | 44 | 96.2% | 92.5-99.3% |

The two distributions overlap almost completely. If these were PCR/library duplicates,
within-sample pairs would cluster near 100% and sit clearly above cross-sample pairs; they
don't. This is what standing population variation should look like: any two reads drawn from
the same underlying pool of cells are about equally similar whether they came from the same
sequencing run or a different one.

## Correction to the earlier claim: it is not cleanly 2 haplotypes

The earlier writeup reported "9 whole-element reads, one shared haplotype" and "5
partial-junction reads, one shared haplotype" based on each group's *best hit* to a single
reference read. Redone as the full 15x15 all-pairs matrix (every read against every other,
identity and query coverage both reported), two things change the picture:

**Two reads have a coverage artefact and should be set aside.** `SRR33298461.44252`
(6991_day0) and `SRR33298384.45869` (with_selection_repeat2) align to every other read at only
~70-74% query coverage, not the ~100% every other pair shows. That is a boundary-calling
issue with how their Y' span was extracted (their coordinates likely under- or over-shot),
not a biological signal -- comparing across a partial, misaligned region produces an
unreliable identity number. Both are removed from the analysis below pending a fix to their
coordinates.

**Among the 13 reads with full coverage, identity is not bimodal -- it is a spread from 90.4%
to 99.8%, with a dominant tight cluster and a few more divergent members:**

|  | tight cluster (7 reads) | close to the cluster (2 reads) | more divergent (3 reads) |
|---|---|---|---|
| reads | 260642(W), 1068481(W), 107087(W), 1079946(W), 46529(W), 797552(W), 896473(W) | 1520960(P), 8939(P) | 125819(P), 267122(P), 272313(W) |
| mutual identity | 97.5-99.8% | 97.3-99.0% to the cluster; 99.0% to each other | 90.4-96.5% to everything |

The "W" (no-junction) and "P" (strong/weak, partial-junction) labels from the window scan do
**not** map cleanly onto this structure. Two P reads (1520960, 8939) sit almost as close to
the W cluster as W reads sit to each other. One W read (272313) sits out with the divergent
group. So the evidence-category split from the window scan and the sequence-similarity split
visible here are not the same partition.

## What this actually supports, stated carefully

* **There is one clearly real, recurring lineage**: the 7-read tight cluster plus its 2 close
  neighbours (9 of 13 full-coverage reads), spanning all 3 samples that have multiple events
  (`reference_promethion`, `with_selection_repeat`, `with_selection_repeat2`) at 97.3-99.8%
  mutual identity. That range is consistent with ONT read-to-read noise (1-3% per read, so
  ~2-6% between two reads of the same true sequence) on a single underlying molecule. This
  part of the earlier claim holds.
* **The 3 more divergent reads (90.4-96.5%) are a separate question the earlier writeup
  answered too confidently.** That spread is wider than ONT noise alone comfortably explains,
  so it is more likely a second, related lineage than pure error -- but with only 3 reads
  (2 from one sample, 1 from another) this is not established with the same confidence as the
  dominant cluster.
* **"Two ancestral molecules" should be read as "at least one clearly recurring lineage,
  plausibly a second, smaller one" -- not two cleanly resolved haplotypes.** The core
  argument this was built to support -- that chr13L->G1 is dominated by recurring standing
  variation rather than 15 independent fresh recombination events -- still holds, because
  9 of 13 usable reads collapse into one lineage regardless of exactly how the remainder
  resolves. But the precise count of distinct ancestral haplotypes is not as settled as
  "exactly 2" implied.

## What would tighten this further

* Re-extract `SRR33298461.44252` and `SRR33298384.45869` with corrected coordinates and
  re-run the comparison -- their current ~70% coverage numbers are uninformative.
* Build a proper multiple alignment (not pairwise best-hit) across all 13 full-coverage reads
  and call variant sites, rather than inferring cluster structure from a pairwise identity
  matrix -- that would settle whether 125819/267122/272313 are a real second lineage or the
  tail of ONT error on the dominant one.
