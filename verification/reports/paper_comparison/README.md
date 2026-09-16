# Comparison with the paper's template-switching / "onion skin" analysis

Two supplementary tables, one row per Y' copy per read:

| file | sheet | strain | populations |
|---|---|---|---|
| `41467_2026_72032_MOESM7_ESM.xlsx` (Supp Data 5) | `Sheet 1` | 6991 (WT) | PD 28 (240 reads), PD 33 (244) |
| `41467_2026_72032_MOESM8_ESM.xlsx` (Supp Data 6) | `mph1 template switching` | 7172, 7302 (mph1Δ) | PD 28, PD 33 |

The `Template Switching (Y/N)` flag is **per read** (one value, on the read's first row) and means
*the donor template changed within the read* -- not "several Y' IDs are present".
`ID7,ID8,ID7,ID8,ID7` at one end is **N** (chr13L copied repeatedly); `ID8,ID7,ID3` is **Y**
(chr13L then chr14L). Our equivalent is the number of donor blocks in `y_prime_path`:
consecutive path segments whose candidate-donor sets intersect are merged, and >= 2 blocks
means the template changed.

## Which of our samples the paper's reads are (by ONT read UUID, not assumption)

| paper | our sample(s) | matched |
|---|---|---|
| 7302 PD 28 | `7302_day3_with_selection` | 103 / 103 |
| 7302 PD 33 | `7302_day4_with_selection` | 86 / 86 |
| 7172 PD 28 | `7172_day3_with_selection` | 60 / 60 |
| 7172 PD 33 | `7172_day4_with_selection` **+** `..._repeat` pooled | 41 + 84 = 125 / 125 |
| 6991 PD 28, PD 33 | **none of our samples** | 0 / 484 |

7172 PD 33 pools both of our day-4 libraries, i.e. the paper treats them as one population --
worth remembering, because chr11L->chr11R is fixed in one and near-absent in the other.

The 484 wild-type reads are not in any of the 33 samples we have. Our 6991 data is only the
eight day-0 references; the paper's WT populations are later time points that were never
downloaded. From the BioProject run table (PRJNA1254968, 133 runs, 56 of them 6991) the
candidates are `SRR33298425` / `SRR33298459` / `SRR33298414` (day 3, ~8 GB) and `SRR33298442` /
`SRR33298440` / `SRR33298441` (day 4, ~7 GB); by the mapping above PD 28 = day 3 and PD 33 = day 4
in both mph1 strains. Downloading and processing those would make a read-level WT comparison
possible.

## Two ways to compare

**A. Read level** (`compare_paper_switching.py`) -- our pipeline's own calls for the same reads.
Needs the reads, so 7302 and 7172 only. Run against the curated-library results.

| | reads found | chr_end agrees | Y' count agrees | switches caught | extra | agreement |
|---|---|---|---|---|---|---|
| 7302 PD 28 | 103/103 | 103 | 102 | 7 / 9 | 7 | 91 % |
| 7302 PD 33 | 83/86 | 83 | 80 | 6 / 7 | 11 | 86 % |
| 7172 PD 28 | 59/60 | 59 | 53 | 9 / 9 | 15 | 75 % |
| 7172 PD 33 | 119/125 | 119 | 114 | 22 / 24 | 17 | 86 % |

**B. Call level** (`compare_paper_paths.py`) -- our path parser applied to the paper's own
per-copy Y' calls. Needs only the Excel plus a day-0 BED and the curated library, so it works
for the wild type too.

| strain | PD | reads | paper switches | caught | missed | extra | agreement |
|---|---|---|---|---|---|---|---|
| 6991 (WT) | 28 | 240 | 42 | 41 | 1 | 25 | 89.2 % |
| 6991 (WT) | 33 | 244 | 59 | 59 | 0 | 17 | 93.0 % |
| 7172 | 28 | 60 | 9 | 9 | 0 | 7 | 88.3 % |
| 7172 | 33 | 125 | 24 | 24 | 0 | 15 | 88.0 % |
| 7302 | 28 | 103 | 9 | 8 | 1 | 4 | 95.1 % |
| 7302 | 33 | 86 | 7 | 6 | 1 | 5 | 93.0 % |

**147 of the paper's 150 switching reads are recovered**, overall agreement 92.2 %.

Two details make this work, both validated:

* *Read orientation.* The sheet lists copies in read coordinates and does not say which end the
  telomere is on. When the telomere is at the read start the array begins almost immediately
  (first copy within 265 bp in every read we can check); otherwise anchor + spacer + X come
  first and the array starts at least 4.5 kb in. A 2 kb threshold separates the two perfectly
  (366/366 where our own `telo_side` is known, and 187/187 + 179/179 in the runs above).
* *Vocabulary.* The reference is relabelled with the entry -> group map taken from the Excel
  itself. For 7302 this changes 3 positions -- the paper calls chr13L's long copies
  `ID8_Brown`, our copy of the curated file calls them `ID2_Red-Light` -- and without it those
  copies are unexplainable and fabricate switches (7302 agreement 61 % -> 95 %). 7172 and 6991
  need no relabelling: our curated files match the paper's vocabulary exactly.

## Variant level is the right granularity

Collapsing the colour shade (ID2_Red-Light + ID2_Red-Dark -> ID2) removes some false calls but
loses far more real ones, across all six strain/PD groups:

| rule | caught (of 150) | missed | extra | overall agreement |
|---|---|---|---|---|
| variant (ID + shade) | **147** | 3 | 73 | **92.2 %** |
| family (ID only) | 98 | 52 | 48 | 89.7 % |
| hybrid (variant and family must agree) | 98 | 52 | 43 | 90.2 % |

So the curated runs' `--y-prime-id-level variant` is the correct setting.

## The one systematic artefact

Most of our extra calls are a single pattern: a chr14L read whose array is its own
(`ID5, ID2, ID3, ID3, ...`) but whose first copy carries the Blue-**Dark** shade (chr7R/chr16L)
instead of Blue-Light (chr14L's own), so its native array reads as foreign --
23 of 42 in the wild type, and the same signature in both mph1 strains. Note the shade there is
assigned in the paper's own per-copy calls, and the paper still scores those reads N, which is
why its switch call must tolerate within-family shade differences. `ID5_Blue-Light` and
`ID5_Blue-Dark` differ by a handful of SNPs over 6.6 kb, at the edge of what a single ONT read
can support.

Files: `paper_vs_ours_per_read.tsv` / `paper_vs_ours_summary.tsv` (A),
`<strain>_paper_paths[.md|.tsv]` and `..._family.*` (B).
