# Column dictionary — `ALL_<strain>_combined.tsv` and `<sample>_annotated_summary.tsv`

One row = **one Y' copy** that landed in a different group from the one that positionally
belongs there. A read with several bad copies contributes several rows, so `read_id` is **not**
a unique key — `(read_id, copy)` is.

| # | column | meaning |
|---|---|---|
| 1 | `strain` | 6991, 7172 or 7302. |
| 2 | `sample` | Which day-0 preparation. Present in the 6991 combined file because it pools eight independently assembled references; each has its own library and its own grouping, so rows are only comparable within a sample. |
| 3 | `read_id` | SRA-style read identifier. Not unique — see the note above. |
| 4 | `chr_end` | The chromosome end the read is **anchored** to, i.e. where it physically comes from. This is the recipient. |
| 5 | `copy` | Which Y' copy in the array, as `n/N`. `3/5` = the third of five. **Use this to index `program_assigned_yprime`** — the *n*-th comma-separated entry is this row's copy. |
| 6 | `expected_element` | What *should* be at that position: the reference element at array index *n* of `chr_end`. The positional truth. |
| 7 | `expected_group` | Its cut99 group label (`G1`…`G13`). Labels are per-sample; `G8` in one sample is unrelated to `G8` in another. |
| 8 | `donor_element` | The element the program actually matched **for this copy**. The mismatch is that `donor_group != expected_group`. |
| 9 | `donor_group` | Its cut99 group. |
| 10 | `program_assigned_yprime` | The pipeline's **whole** assigned array for the read, all copies, comma-separated. Equals `donor_element` only on single-copy reads. |
| 11 | `program_status` | The pipeline's own recombination verdict. Nearly always `1st Y' Change` (159 of 160 non-defect rows), so it does **not** discriminate between our evidence classes. |
| 12 | `program_source` | The chromosome end the pipeline attributed the Y' to. |
| 13 | `program_mechanism` | Its mechanism label: `donor_transfer`, `subtelomere_switch`, `unmatched_array`. |
| 14 | `program_path` | Its reconstructed origin path, e.g. `chr14L[1]:E-chr14L-1`. |
| 15 | `anchorHalf_vs_expected` | % identity of the **anchor-side half** of the read's Y' to the expected element. Should be high if the anchor side is still native. |
| 16 | `anchorHalf_vs_donor` | Same half against the donor. Should be lower. |
| 17 | `anchor_margin` | `anchorHalf_vs_expected − anchorHalf_vs_donor`. **Positive** = the anchor side favours the read's own reference. |
| 18 | `teloHalf_vs_expected` | % identity of the **telomere-side half** to the expected element. Should be lower. |
| 19 | `teloHalf_vs_donor` | Same half against the donor. Should be higher. |
| 20 | `telo_margin` | `teloHalf_vs_donor − teloHalf_vs_expected`. **Positive** = the telomere side favours the donor. |
| 21 | `pair_longest_homology_bp` | Longest contiguous ≥95 %-identity block between the expected and donor **elements** — how much homology a crossover had to work with. A property of the reference pair, not of the read. |
| 22 | `pair_homology_identity` | Identity of that block. Typically 96–98 %. |
| 23 | `evidence` | The verdict. See below. |
| 24 | `explanation` | One-sentence reading of `evidence`. |
| 25 | `telo_side` | `beginning` or `end` — which end of the read the telomere is on. Needed to convert read coordinates into anchor-side/telomere-side; already applied in columns 15–20. |

## `evidence` values

| value | meaning |
|---|---|
| `strong` | Both margins ≥ 1.5 %: anchor half favours the expected element **and** telomere half favours the donor. The mid-Y' recombination signature. **The only class I would treat as established.** |
| `weak` | Correct direction, margins < 1.5 %. Tracks low read quality as much as biology — suggestive only. |
| `FAILS` | No half is better explained by the donor. The copy misses its group but a mid-Y' switch does not explain it. Of 31 such rows, 13 have whole-read BLAST clearly favouring the *expected* element, i.e. a RepeatMasker mis-assignment rather than biology. |
| `no junction` | The donor wins across the whole element, nothing to split. **Not necessarily a distinct mechanism** — where donor and recipient are ~99 % identical a junction anywhere gives the same sequence and cannot be located. |
| `reference defect` | Known mis-assembly, not recombination. All 388 are `chr14L-1` in `6991_day0_with_selection` (5,720 bp vs 6,654 elsewhere). |

Columns 15–22 are blank for `no junction` and `reference defect` rows — there is no junction to
split at, so the per-half comparison is undefined.

## Reading a row

```
chr_end=chr14R  copy=1/1  expected=chr14R-1 (G10)  donor=chr14L-1 (G8)
anchor half: 99.71 % vs chr14R-1 / 97.36 % vs chr14L-1   → +2.35, own reference wins
telo   half: 95.48 % vs chr14R-1 / 99.27 % vs chr14L-1   → +3.79, donor wins
evidence=strong
```

A read anchored at chr14R whose Y' is native chr14R-1 on the anchor side and chr14L-1-like on
the telomere side — a junction part-way through the element. The program reported the whole Y'
as chr14L (`program_source=chr14L`), having no way to express a partial switch.

## Two caveats worth carrying

* **The donor is usually a group, not an element.** `donor_element` is the best-scoring member
  of a set of near-identical sequences; several others would score equally. Check the group
  membership in `<sample>_groups.json` before naming a specific donor.
* **Group labels are per-sample.** They come from each reference's own clustering, so never
  compare a `G` label across samples.
