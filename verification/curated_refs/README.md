# Curated references (copied from TeloTracker/references/<strain>_features/ on Argon)

One local correction:

* `7172_features/repeatmasker_7172_all_y_primes.fasta`: header `>Y_Prime_chr5R11` -> `>Y_Prime_chr5R1`.
  chr5R carries exactly one Y' in both the curated 7172 BED and our 7172 day-0 reference, and the
  same entry is spelled `chr5R1` in the 6991 and 7302 curated libraries. Left as `chr5R11` the
  library declares a Y' at chr5R position 11 and none at position 1, which makes the strict-library
  check abort (a BED Y' that cannot be resolved silently turns every read at that end into a
  "1st Y' Change"). The file on Argon is unchanged.

## Known defect: `Y_Prime_chr16L1` is 82 bp over-extended (6991)

`repeatmasker_6991_all_y_primes.fasta` gives `Y_Prime_chr16L1` as **6737 bp** while its
near-identical partners are 6656 (`chr7R1`) and 6657 (`chr14L1`) — the three are 99.94–99.97 %
identical over their shared 6655 bp, i.e. 1–3 bp apart. The excess is at chr16L1's **5'
(anchor-proximal) end**: BLASTing it against every partner ≥99 % identical puts the shared start at
query position 83 in each case, and the leading 82 bp instead align into what chr7R has annotated
as `x_variable_element`. That region measures **0 % telomeric repeat**, whereas every genuine
Y'-to-Y' ITS in the 6991 reference measures 0.89–1.00.

Consequences, measured on 52,182 day-0 reads (`verification/reports/grouping_benchmark/`):

* the labelling step reproduces the boundary in every assembly (6731 ± 1 bp in all eight
  independently built 6991 day-0 references);
* the clustering similarity is coverage-penalised including terminal overhangs, so 77 bp of end
  offset drops a 99.98 %-identical pair to 98.8 % and splits chr16L into a group of its own at the
  99.9 % dedup and any 99 % cut;
* because the entry carries 77 bp of real subtelomeric sequence that reads also carry,
  RepeatMasker prefers it — **806/806 chr7R reads and 192/204 chr14L anchor-proximal copies match
  `chr16L1` instead of their own element**, giving a 6.9–7.0 % false recombination-call rate at the
  `cut99` / `condensed` resolutions (0.16 % at the shipped silhouette default, where all three
  already share `ID1_Gray`).

**The original files here are unchanged and remain the reference for the Supplementary Data 5/6
comparisons**, which were all run against them. A corrected copy with the 82 bp removed is at
`_pipeline/references/repeatmasker_6991_all_y_primes.chr16Lfix.fasta` (same 23 entries, chr16L1
6737 → 6655 bp); no pipeline default points at it. The equivalent correction to the derived
reference BEDs is made by `verification/fix_yprime_boundary.py`, which measures the offset per
assembly rather than assuming it (77 bp on all eight 6991 day-0 references, 10–17 partners agreeing).

Not fixable by a boundary correction: the curated scheme calls `chr7R1`/`chr16L1` `ID5_Blue-Dark`
and `chr14L1` `ID5_Blue-Light`, a distinction resting on 1–2 bp across 6.6 kb — below ONT read
accuracy, so no grouping scheme can read it reliably from a read. This is the source of most of the
spurious wild-type switch calls in the Supp Data 5 comparison.
