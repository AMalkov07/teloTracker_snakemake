# Comparison with Supplementary Data 6 (template switching) — 7302 day 4

Source: `41467_2026_72032_MOESM8_ESM.xlsx`, sheet "mph1 template switching"
(Tsai et al. 2026; 2245 rows = one row per Y' copy per read, strains 7172 + 7302, PD 28 and 33).

## PD <-> our sample (established by read identity, not assumption)

Matching the 189 7302 ONT read UUIDs against every sample's `<sample>_read_id_map.tsv`:

| paper PD | our sample | reads matched |
|---|---|---|
| 28 | `7302_day3_with_selection` | 103 / 103 |
| 33 | `7302_day4_with_selection` | 86 / 86 |

No read matched any other 7302 sample, so **PD 33 = day 4** (and PD 28 = day 3).

## Read-level agreement (PD 33 = day 4)

* 84 / 86 paper reads are present in our day-4 features; **chr_end agrees for 84/84**.
* Y' copy count identical for 80 / 84 (we see 1–4 fewer copies in 4 long reads).
* The 2 absent reads: `528c9417-…` (chr14R) is anchored in our run but dropped by the
  telomere filter (no qualifying telomere repeat); `6d347b8d-…` (chr2L) never passed the
  anchor BLAST filter.

## Their "template switching" flag vs our path prediction

The flag is **per read** (one value on the first row of each read), 7 of 86 at day 4.
Their criterion is a *change of donor template within the read*, not merely several Y' IDs:
`ID7,ID8,ID7,ID8,ID7` at chr10R is **N** (one template, chr13L, copied repeatedly), while
`ID8,ID7,ID3` is **Y** (chr13L then chr14L). That is exactly what `yprime_path.py` encodes,
so the comparable statistic is "the path has >= 2 donor blocks" (consecutive segments whose
candidate-donor sets intersect are merged).

| input to our path parser | paper Y caught | paper N agreed | agreement |
|---|---|---|---|
| the paper's own per-copy Y' calls, oriented by our `telo_side` | **7 / 7** | 76 / 79 | **97 %** |
| our pipeline's output, curated library | 5 / 7 | 72 / 79 | 93 % |
| our pipeline's output, silhouette library | 1 / 7 | 74 / 79 | 90 % |

The three reads we call a switch and they do not look like genuine switches:
`chr2R ID3,ID2,ID2,ID2,ID2` (chr14L then chr12R/chr4R), `chr2L ID3,ID2,ID2`,
`chr12R ID1,ID2x6,ID3,ID3` (own array then two chr14L short copies).

## Why our own runs catch fewer: the library labels, not the matching

The paper labels the library entry `Y_Prime_chr13L4` as its own group **ID8_Brown**; our copy of
`repeatmasker_7302_all_y_primes.fasta` labels the *same entry* `ID2_Red-Light`, together with the
chr4R/chr12R long copies. With the paper's labels chr13L reads `ID7,ID8,ID7,ID8` — both IDs unique
to chr13L — which is why their chr13L-derived switches are unambiguous.

**RepeatMasker already makes this distinction in our runs**: both libraries carry chr13L's long
copies as separate entries (`Y_Prime_chr13L2` / `Y_Prime_chr13L4` in the curated library,
`Y_Prime_chr13L2,4` in the silhouette library). We discard it when we collapse the matched entry
to its group ID. Keeping the entry name in `y_prime_positions` (or an entry-level
`--y-prime-id-level`) would recover the paper's resolution without changing any clustering;
it needs one re-run of step 11 because the entry name is not stored in the current TSVs.

**Caveat on that resolution.** chr13L2 and chr13L4 are 100 % identical to each other, 99.955 % to
`chr12R6` (3 differences in 6652 bp) and 99.85 % to the chr4R entries, and the SW-score
distributions of their ID2 and ID8 copies overlap almost completely. In the paper's own data, at
the chr13L end 21 of 26 long copies are assigned to chr13L4 but 5 to chr4R entries — a ~19 %
per-copy error rate, so a two-copy chr13L fingerprint is right ~65 % of the time on entry
assignment alone. The ITS lengths are independent of this and much more reliable (chr13L 152/163 bp
vs chr4R 10 bp), which is why the (ID, ITS) path agrees with the paper at 97 % when fed the same
per-copy calls.

Files: `day4_switch_comparison.tsv` (per read: paper array, our blocks, both flags),
`day4_perread_join.tsv` (full join of the paper's rows with our day-4 features).
