# Recipient -> donor relationships, per dataset

For each sample: which chromosome end was mismatched (recipient), and which donor GROUP
its Y' was matched into instead. Donor is given as a group because at this cutoff the
pipeline genuinely cannot separate a group's members -- naming one element would overstate
precision. Group labels are local to each sample.

The defect rows in 6991_day0_with_selection (chr14L-1, 388 copies) are shown but marked,
since they are an assembly artefact rather than a real recipient->donor relationship.

## 6991_day0  (4 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr10L | 1 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 100% |
| chr6L | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |
| chr13L | 1 | G1 (chr2L-1, chr6L-1) | 1 | 100% |
| chr14L | 1 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 100% |

## 6991_day0_TeloTag  (6 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr14R | 3 | G9 (chr14L-1, chr16L-1, chr7R-1) | 2 | 67% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 33% |
| chr7R | 1 | G10 (chr5R-1) | 1 | 100% |
| chr13L | 1 | G3 (chr8R-1) | 1 | 100% |
| chr16R | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |

## 6991_day0_TeloTag_with_selection  (11 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr5L | 2 | G11 (chr12R-1) | 2 | 100% |
| chr16R | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 2 | 100% |
| chr2L | 1 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 100% |
| chr7R | 1 | G10 (chr14R-1) | 1 | 100% |
| chr8R | 1 | G10 (chr14R-1) | 1 | 100% |
| chr14L | 1 | G11 (chr12R-1) | 1 | 100% |
| chr8L | 1 | G1 (chr2L-1, chr6L-1) | 1 | 100% |
| chr13L | 1 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 100% |
| chr16L | 1 | G11 (chr12R-1) | 1 | 100% |

## 6991_day0_reference  (27 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr14R | 5 | G7 (chr10L-1, chr9L-1) | 3 | 60% |
|  |  | G12 (chr5L-1) | 1 | 20% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr14L-1, chr14L-...) | 1 | 20% |
| chr13L | 4 | G1 (chr2L-1, chr6L-1) | 2 | 50% |
|  |  | G7 (chr10L-1, chr9L-1) | 2 | 50% |
| chr6L | 4 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 4 | 100% |
| chr12L | 3 | G6 (chr16R-1) | 1 | 33% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 33% |
|  |  | G1 (chr2L-1, chr6L-1) | 1 | 33% |
| chr2L | 3 | G3 (chr8R-1) | 2 | 67% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 33% |
| chr16L | 2 | G12 (chr5L-1) | 1 | 50% |
|  |  | G7 (chr10L-1, chr9L-1) | 1 | 50% |
| chr5L | 2 | G9 (chr5R-1) | 1 | 50% |
|  |  | G11 (chr12R-1) | 1 | 50% |
| chr16R | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |
| chr8L | 1 | G1 (chr2L-1, chr6L-1) | 1 | 100% |
| chr15R | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |
| chr8R | 1 | G1 (chr2L-1, chr6L-1) | 1 | 100% |

## 6991_day0_reference_promethion  (12 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr13L | 4 | G1 (chr2L-1, chr6L-1) | 4 | 100% |
| chr16R | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 2 | 100% |
| chr14R | 2 | G7 (chr10L-1, chr9L-1) | 1 | 50% |
|  |  | G1 (chr2L-1, chr6L-1) | 1 | 50% |
| chr8L | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |
| chr8R | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |
| chr2L | 1 | G3 (chr8R-1) | 1 | 100% |
| chr5R | 1 | G7 (chr10L-1, chr9L-1) | 1 | 100% |

## 6991_day0_with_selection  (439 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr14L **[defect]** | 388 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 388 | 100% |
| chr13L | 18 | G1 (chr2L-1, chr6L-1) | 12 | 67% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 4 | 22% |
|  |  | G3 (chr8R-1) | 1 | 6% |
|  |  | G9 (chr5R-1) | 1 | 6% |
| chr14R | 8 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 4 | 50% |
|  |  | G7 (chr10L-1, chr9L-1) | 4 | 50% |
| chr10L | 5 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 3 | 60% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 20% |
|  |  | G9 (chr5R-1) | 1 | 20% |
| chr5R | 3 | G7 (chr10L-1, chr9L-1) | 3 | 100% |
| chr6L | 3 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 3 | 100% |
| chr2L | 3 | G3 (chr8R-1) | 3 | 100% |
| chr16L | 2 | G7 (chr10L-1, chr9L-1) | 2 | 100% |
| chr16R | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 2 | 100% |
| chr8L | 2 | G1 (chr2L-1, chr6L-1) | 1 | 50% |
|  |  | G7 (chr10L-1, chr9L-1) | 1 | 50% |
| chr8R | 2 | G1 (chr2L-1, chr6L-1) | 2 | 100% |
| chr12L | 2 | G1 (chr2L-1, chr6L-1) | 1 | 50% |
|  |  | G5 (chr8L-1) | 1 | 50% |
| chr7R | 1 | G3 (chr8R-1) | 1 | 100% |

## 6991_day0_with_selection_repeat  (18 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr13L | 9 | G1 (chr2L-1, chr6L-1) | 8 | 89% |
|  |  | G5 (chr8L-1) | 1 | 11% |
| chr10L | 3 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 3 | 100% |
| chr14R | 3 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 33% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 33% |
|  |  | G7 (chr10L-1, chr9L-1) | 1 | 33% |
| chr16R | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |
| chr2L | 1 | G3 (chr8R-1) | 1 | 100% |
| chr5L | 1 | G7 (chr10L-1, chr9L-1) | 1 | 100% |

## 6991_day0_with_selection_repeat2  (13 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr13L | 3 | G1 (chr2L-1, chr6L-1) | 2 | 67% |
|  |  | G6 (chr16R-1) | 1 | 33% |
| chr10L | 2 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 2 | 100% |
| chr16R | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 2 | 100% |
| chr12L | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 2 | 100% |
| chr2L | 2 | G3 (chr8R-1) | 2 | 100% |
| chr14R | 1 | G1 (chr2L-1, chr6L-1) | 1 | 100% |
| chr15R | 1 | G7 (chr10L-1, chr9L-1) | 1 | 100% |

## 7172_day0_with_selection  (4 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr16R | 2 | G2 (chr13L-1, chr14L-3) | 2 | 100% |
| chr5R | 1 | G2 (chr13L-1, chr14L-3) | 1 | 100% |
| chr6L | 1 | G2 (chr13L-1, chr14L-3) | 1 | 100% |

## 7302_day0_with_selection  (14 mismatched copies)

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr6L | 7 | G2 (chr13L-1, chr13L-3, chr14L-3, chr14L-4, chr14L-5) | 4 | 57% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr13L-2, chr13L-...) | 1 | 14% |
|  |  | G10 (chr14R-1) | 1 | 14% |
|  |  | G3 (chr8R-1) | 1 | 14% |
| chr14R | 2 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr13L-2, chr13L-...) | 2 | 100% |
| chr5R | 2 | G10 (chr14R-1) | 1 | 50% |
|  |  | G7 (chr10L-1, chr9L-1) | 1 | 50% |
| chr10L | 1 | G10 (chr14R-1) | 1 | 100% |
| chr2L | 1 | G3 (chr8R-1) | 1 | 100% |
| chr8R | 1 | G1 (chr2L-1, chr6L-1) | 1 | 100% |

---

# Combined: the 5 samples with an identical cut99 grouping  (58 copies)

Samples: 6991_day0, 6991_day0_TeloTag_with_selection, 6991_day0_reference_promethion, 6991_day0_with_selection_repeat, 6991_day0_with_selection_repeat2

Group labels are identical across all five, verified by direct membership comparison,
so these counts pool without any label-translation.

| recipient | n | donor group (members) | count | % of recipient |
|---|---|---|---|---|
| chr13L | 18 | G1 (chr2L-1, chr6L-1) | 15 | 83% |
|  |  | G6 (chr16R-1) | 1 | 6% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 6% |
|  |  | G5 (chr8L-1) | 1 | 6% |
| chr16R | 7 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 7 | 100% |
| chr10L | 6 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 6 | 100% |
| chr14R | 6 | G1 (chr2L-1, chr6L-1) | 2 | 33% |
|  |  | G7 (chr10L-1, chr9L-1) | 2 | 33% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 17% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 17% |
| chr2L | 5 | G3 (chr8R-1) | 4 | 80% |
|  |  | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 20% |
| chr5L | 3 | G11 (chr12R-1) | 2 | 67% |
|  |  | G7 (chr10L-1, chr9L-1) | 1 | 33% |
| chr8R | 2 | G10 (chr14R-1) | 1 | 50% |
|  |  | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 50% |
| chr8L | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 50% |
|  |  | G1 (chr2L-1, chr6L-1) | 1 | 50% |
| chr12L | 2 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 2 | 100% |
| chr14L | 2 | G8 (chr12R-2, chr12R-3, chr12R-4, chr12R-5, chr12R-6, chr12R-7, chr14L-...) | 1 | 50% |
|  |  | G11 (chr12R-1) | 1 | 50% |
| chr6L | 1 | G2 (chr13L-1, chr14L-3, chr14L-4, chr14L-5) | 1 | 100% |
| chr7R | 1 | G10 (chr14R-1) | 1 | 100% |
| chr15R | 1 | G7 (chr10L-1, chr9L-1) | 1 | 100% |
| chr16L | 1 | G11 (chr12R-1) | 1 | 100% |
| chr5R | 1 | G7 (chr10L-1, chr9L-1) | 1 | 100% |

---

# Is the relationship reproducible across datasets?

For each recipient, the dominant donor group in each dataset it appears in, expressed as that
group's set of chromosome ends (labels differ per sample, end-sets do not):

| recipient | datasets | dominant donor end-set, per dataset | verdict |
|---|---|---|---|
| **chr16R** | 8 | `{chr13L, chr14L}` x8 | **identical in all 8** |
| **chr6L** | 5 | `{chr13L, chr14L}` x5 | **identical in all 5** |
| **chr2L** | 7 | `{chr8R}` x6, `{chr12R,chr14L,chr15R,chr16L,chr4R,chr7R}` x1 | 6 of 7 |
| **chr13L** | 8 | `{chr2L, chr6L}` x6, `{chr8R}` x1, `{chr12R,...}` x1 | 6 of 8 |
| chr10L | 5 | `{chr12R,chr14L,chr15R,chr16L,chr4R,chr7R}` x4, `{chr14R}` x1 | 4 of 5 |
| chr8L | 4 | `{chr2L, chr6L}` x3, `{chr13L, chr14L}` x1 | 3 of 4 |
| chr8R | 5 | `{chr2L, chr6L}` x3, `{chr14R}` x1, `{chr13L, chr14L}` x1 | 3 of 5 |
| chr5R | 4 | `{chr10L,chr9L}` x2, `{chr13L,chr14L}` x1, `{chr14R}` x1 | 2 of 4 |
| chr14R | 7 | five different end-sets | **no consistent donor** |
| chr14L, chr7R, chr5L, chr16L, chr12L, chr15R | 2-3 | all different | too few events |

**chr16R -> {chr13L, chr14L} and chr6L -> {chr13L, chr14L} hold in every single dataset they
appear in** -- 8 and 5 independent samples respectively, spanning all three strains. These are
not sample-specific quirks.

`chr2L -> {chr8R}` and `chr13L -> {chr2L, chr6L}` hold in 6 of 7 and 6 of 8 datasets. Note
chr13L and chr6L point at each other's groups, and chr2L and chr8R point at each other's
groups -- two reciprocal pairs.

chr14R is the clear negative: 24 events across 7 datasets with five different dominant donor
end-sets and no group exceeding 38 %. Whatever drives the other recipients does not apply to it.

## Caveats

* Most per-dataset cells rest on 1-3 events. The "identical in all N" claims for chr16R and
  chr6L are strong because they never deviate, but several other rows are single-event
  observations that would not survive more data.
* `6991_day0_with_selection`'s chr14L row is the assembly defect (chr14L-1 truncated to
  5,720 bp), not a real relationship; it is marked in the per-dataset table above.
* Donor is a group, not an element -- `{chr13L, chr14L}` means the pipeline cannot separate
  chr13L-1 from chr14L-3/4/5 at this cutoff, so "chr16R was mistaken for one of those four".
