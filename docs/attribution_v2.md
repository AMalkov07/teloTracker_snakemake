# Recombination source attribution, version 2

`analyze_features.py` (step 11) decides, for every telomere-anchored read, whether its
subtelomere differs from the day-0 reference and, if so, **which other chromosome end the
new material came from** (`recombination_source`) and by what mechanism. This document
describes the v2 attribution introduced in September 2026, why it was needed, and what changed
in the output. The previous behaviour is still available with `--attribution-mode legacy`
(Snakemake: `attribution_mode: legacy` in `config.yaml`).

## Why

Verification against the 7302 time course showed that the legacy attribution was blind to
the most interesting class of event, **Y'-only gains**:

* only the spacer, the x-element and the supplementary alignment voted for the source; the
  Y' array never did, and only the *first* Y' ID of a read was ever compared with the
  reference ends;
* the day-0 chr13L end of strain 7302 carries a unique alternating short/long array
  (`ID2,ID1,ID2,ID1`). Reads at other ends that had gained that exact array were attributed
  to chr13L in **1.4 %** of cases (1 of 72); the rest went to whichever chromosome the
  supplementary alignment happened to hit (`chr2`, `chr5`, ...), or to `ambiguous`;
* the spacer switch confidence was numerically meaningless (−0.02 … 0.05 for every switch in
  every sample) because it compared mean identities of two *different* chunk sets;
* the spacer chunk walk ran over the **whole read**, and the "spacer" library
  (`pairings_for_spacers/`, built with `--fixed-50kb`: 50 kb inward from each telomere) is
  mostly Y' sequence for tandem-array ends, so reads with a large Y' array were called
  "spacer switch to chr12R / chr4R / chr14L" from chunks inside their Y' array (in 7302 day 5,
  9 of 12 spacer switch positions lay inside the Y' region). Once the confidence was repaired
  this artefact started out-voting the Y' fingerprint, which is how it was noticed. The 50-kb
  window also did not reach the spacer at all when the array exceeded 50 kb (7172 chr12R,
  ~60 kb), giving 51 spurious chr12R→chr3R switches at day 0;
* several smaller bugs: arm-less supplementary contigs (`chr4`) could not be compared with
  arm-resolved ends (`chr4R`, flagged as "complex"); `'chr1' in 'chr10L'` substring tests;
  ties resolved by dictionary order; a spacer with a real switch was labelled with the
  plurality end instead of the post-breakpoint end; a Y' library / BED mismatch silently
  turned every read into a "1st Y' Change"; `most_common_source` in the summary counted the
  literal `ambiguous`.

## What v2 does

1. **Full-array fingerprint** (`find_fingerprint`). The gained segment of the read's Y'
   array (everything from the divergence point to the telomere) is matched against every
   reference end's array: `contiguous` (verbatim piece of that end), `rotation` (a cyclic
   rotation — a circle excised from that end and re-inserted in another phase), or
   `periodic` (repeats of a unit that is a piece of that end — rolling-circle
   amplification). Rotations rank equal to contiguous matches; ties are broken by how much
   of the window the end carries verbatim, so `ABABAB` goes to the end that holds `ABAB`
   rather than to one with a lone `AB`. Ragged reads are handled by trying shorter windows
   (structured windows before homopolymer runs). A single Y' is a fingerprint only when its
   variant exists at exactly one end (weight 0.5) — e.g. `ID7` at chr13L in the curated 7302
   library. If the read's **own** end explains the gained segment, the event is
   `tandem_amplification_same_end` and the source is the end itself.
2. **Weighted votes** (`reconcile_features_v2`): spacer 1.0, unique Y' fingerprint 1.0
   (≥ 3 Y') / 0.6 (2 Y') / 0.5 (1 unique Y'), 2–3 Y' candidate ends 0.3 shared, x-element
   0.5, supplementary 0.4 (0.25 when arm-less). Ties are broken deterministically (spacer >
   Y' > x-element > supplementary) and recorded (`source_tie`).
3. **Proximal donor first.** A *confident* spacer / x-element switch (confidence ≥ 0.2)
   names the source, because the breakpoint is anchor-side of the Y' array. The Y'
   fingerprint is then reported separately as `y_prime_donor`; if it points elsewhere the
   event is `is_complex_event` (the donor end's own array had changed before the transfer —
   e.g. 7172 chr11L→chr11R reads whose Y' came from chr14L).
4. **Arm-less evidence never wins** over an arm-resolved end; it is folded into the same
   chromosome's candidate if one exists, otherwise kept with `source_resolution=chromosome`.
   A supplementary hit to the read's *own* chromosome is not evidence of a donor.
5. **Spacer walk restricted to the spacer interval** of the read (anchor → x element / Y'
   start; `spacer_interval`, `chunks_in_interval`), so Y' chunks can no longer produce spacer
   switches. **Spacer switch confidence** = 0.9 × min(1, mean per-chunk identity advantage of
   the new end over the expected end / 10). The spacer source is the post-breakpoint segment
   (`spacer_plurality_source` keeps the old value).
6. **Y' Loss confirmation.** With `--telo-tsv/--probe-tsv` (passed by the Snakefile) every
   read gets `telomere_end_confirmed` (adapter after the telomere and ≥ 30 bp of repeat —
   the `read_summary.tsv` "qualifying" definition), `telomere_repeat_length`,
   `y_prime_probe_count` and `y_prime_tail_bp`. A Loss on an unconfirmed end keeps its
   status but its confidence is halved and `qc_flags` gains `loss_unconfirmed_end`.
7. **Mechanism tag** `recombination_mechanism`: `donor_transfer`,
   `donor_transfer_candidates:<n>`, `tandem_amplification_same_end`, `subtelomere_switch`,
   `array_contraction`, `array_contraction_unconfirmed`, `single_y_prime_change`,
   `unmatched_array`.
8. **(ID, ITS) path** (`yprime_path.py`). Every Y' array is written as (variant ID, length of
   the ITS that follows) tokens; ITS lengths measured from the RepeatMasker hit gaps reproduce
   the reference ITS to ±2 bp, so they separate ends whose Y' IDs are identical (a run of long
   Y' with 10-bp ITS is chr4R, with ~172-bp ITS chr12R). The gained part of the read is parsed
   into the fewest donor segments that explain it (dynamic programme over all splits, not a
   greedy longest-first walk); each segment is a piece of one end's array, either linear or as
   a circle copied in tandem (any phase; the junction ITS must be the same at every repeat and,
   when it equals the ITS flanking the piece in the donor, it counts as verified — real chr13L
   circles carry chr13L's 152/163-bp ITS at their junctions). Linear explanations are
   preferred over circles at equal evidence. Each piece is reported with its copy positions in
   the donor (`chr13L[1-2]`), and each circle with a support level: `strong` = more copies
   than the donor array holds (no single copy can produce it), `weak` = one verbatim copy of
   the donor explains it equally well (the linear reading is given as `alt`), `moderate` =
   otherwise (a two-piece linear reading, i.e. two events or a whole-array circle, is shown as
   `alt` when it exists). A piece of a single Y' copy whose ID exists at more than one end is
   only `tentative` (`chr13L?[1]`) and never names the donor. Output: `y_prime_path` (e.g.
   `chr13L[1-2]:ID2,ID1,ID2,ID1,ID2(circ x2.5 strong) > chr4R[1-2]:ID1,ID1`),
   `y_prime_path_primary_donor`, `y_prime_path_circles` (`chr13L[1-2]:ID2,ID1x2.5:strong`),
   ITS verified/checked counts. ITS are compared by length only (±8 bp), not by sequence. The
   path's primary donor is used as the Y' vote when the ID-only fingerprint is not unique (ITS
   as tie-breaker): on 7302 day 5 it resolves ~40 additional multi-Y' gains on top of the 100
   resolved by IDs alone. `verification/path_report.py` tabulates donors and circles;
   `plot_yprime_copies.py` prints the path in each row's label.
9. **Library guard.** Reference Y' features that cannot be resolved in the Y' library abort
   the run (`--strict-lib`, default in v2; `--no-strict-lib` to override). The committed
   `config.yaml` no longer carries a `y_prime_lib_override`.
10. `--reprocess-tsv` re-derives the Y' comparison, the path and the reconciliation from an
   existing `*_features.tsv` without BLAST / RepeatMasker (everything except the spacer
   confidence, spacer interval and post-breakpoint spacer source, which need the chunk hits).

### Y' ID granularity

`--y-prime-id-level family|variant` (config `y_prime_id_level`). The pipeline's own libraries
have one level (colour is a function of the ID number). Curated libraries carry a second level,
the colour shade (`ID2_Red-Light` vs `ID2_Red-Dark`, a few SNPs apart over 6.6 kb); `variant`
keeps that level as the ID used in arrays, fingerprints and paths. On the 7302 curated pilot
(family level) the chr13L truth set reached 96.9 % (31/32) vs 92.3 % with the pipeline library,
because chr13L's short copy is its own variant (`ID7`) there.

### Spacer library

`label_regions.sh` now builds the spacer library from the actual `space_between_anchor` regions
(the script's default). The earlier `--fixed-50kb` window (50 kb inward from each telomere)
was mostly Y' sequence for tandem-array ends and did not even reach the spacer when the array
exceeded 50 kb (7172 chr12R, ~60 kb), which produced spurious "spacer switches" — see the
`v2b`/`v2c` rows below. Spacer-only libraries for the 7302 and 7172 day-0 references are in
`verification/spacer_lib_variable/` (1108 chunks each, no Y' sequence).

### New columns

`*_features.tsv`: `spacer_plurality_source`, `y_prime_gained_segment`,
`y_prime_fingerprint_window`, `y_prime_fingerprint_len`, `y_prime_fingerprint_kind`,
`y_prime_array_matches`, `y_prime_self_match`, `y_prime_fingerprint_source`,
`y_prime_fingerprint_specificity`, `source_votes`, `source_tie`, `source_resolution`,
`recombination_mechanism`, `y_prime_donor`, `telomere_end_confirmed`,
`telomere_repeat_length`, `y_prime_probe_count`, `y_prime_tail_bp`, `y_prime_path`,
`y_prime_path_n_segments`, `y_prime_path_primary_donor`, `y_prime_path_primary_len`,
`y_prime_path_circles`, `y_prime_path_its_verified`, `y_prime_path_its_checked`.
No existing column was renamed or removed.

`*_recombination_summary.tsv`: `n_y_prime_gain`, `n_y_prime_loss`,
`n_y_prime_loss_confirmed_end`, `n_first_y_prime_change`, `n_y_prime_recombination`,
`n_ambiguous`, `most_common_mechanism`, `most_common_source_legacy`.
`most_common_source` now excludes `ambiguous`.

## Run history on Argon (results/<sample>/_pipeline/)

| dir | code | note |
|---|---|---|
| `recombination_v1/` | legacy | the original September 2026 runs |
| `recombination_v2a/` | v2 without the spacer-interval restriction | positive control fell to 85 % and the truth set to 74 % because the newly confident (but spurious) spacer switches out-voted the X element / Y' fingerprint |
| `recombination_v2b/` | v2 + spacer walk restricted to the spacer interval (`recomb_v2b`) | exposed the 50-kb library gap: 51 spurious chr12R→chr3R spacer switches at 7172 day 0 |
| `recombination/` | v2 + spacer-only library + (ID, ITS) path columns (`recomb_v2c`) | current |
| `results/<sample>__curatedYP/` | v2 + curated strain Y' library (`y_prime_lib_override`) + path | pilot on 7302 day 0 / 4 / 5 (family-level IDs) |

## Effect on results

| check | legacy (v1) | v2c, pipeline library | v2 + curated 7302 library (pilot) |
|---|---|---|---|
| 7302 chr13L alternating-array truth set, strict (a full period): reads → chr13L | 1/39 (2.6 %) | 36/39 (92.3 %) | 31/32 (96.9 %) |
| same, all windows ≥ 3 | 1.5 % (68 reads) | 85.3 % (68 reads) | 95.9 % (49 reads) |
| 7172 chr11L → chr11R positive control (day 4 / 6 / 9) | 97 / 91 / 91 % | 97 / 92.5 / 96.2 % | — |
| 7302 day-0 self-run, reads called recombinant with a non-Loss Y' status | 0.57 % | 0.61 % (86 / 14202) | 0.96 % (137 / 14202; 53 of the extra are chr15R reads whose `ID2_Red-Dark` copy RepeatMasker labels `ID5_Blue`, 99.2 % identical) |
| 7172 day-0 spacer switches | 16 | 38 (v2b with the 50-kb library: 83, of which 51 chr12R→chr3R) | — |
| arm-less `most_common_source` values (`chr13`, `chr14`, ...) | common | none | none |
| 7302 day-5 multi-Y' gains with a unique donor (IDs / IDs+ITS / unresolved / self) | — | 100 / 40 / 39 / 23 of 202 | 117 / 26 / 37 / 23 of 203 |

The curated library resolves more donors (chr13L's short copy is its own variant, `ID7`) but
costs per-read noise where curated variants are nearly identical: on single ONT reads
RepeatMasker sometimes prefers a sister variant of another family, which the array
comparison then reports as a "1st Y' Change" (chr15R above). A library clustered at a high
identity threshold (~99.5 %, keeping `ID7`-like distinct variants but merging 99.9 %-identical
sisters) is the likely sweet spot; not done yet.

Detection (`recombination_detected`, `y_prime_recombination_status`) is identical between
legacy and v2 for every read except reads whose only "evidence" was a spacer switch (called
from Y' chunks in v1, or from a spacer chunk in v2c); only the *source*, its confidence, and
the new columns change otherwise. All numbers above are from the Argon runs, reproduced
locally by `--reprocess-tsv` (Part B reports in `verification/reports/`).

Things to state when reporting results produced with v2:

* `ambiguous` counts and `most_common_source` differ from earlier tables (a supplementary
  hit to the read's own chromosome is no longer a donor; arm-less donors are resolved or
  folded).
* the day-0 "null" is not zero at tandem-array ends: 4–11 % of day-0 reads at chr4R, chr12R,
  chr13L and chr14L already carry fewer Y' copies with a confirmed telomere end
  (sub-clonal copy-number heterogeneity), and ≈ 0.6 % of reads genome-wide show an extra or
  changed Y'. Time points are compared against that per-end baseline.
* day-6 chr13L in 7302 is a real contraction (`ID2,ID1,ID2,ID1` → `ID2,ID1`; 85 of 95 Loss
  reads have a confirmed telomere end), not a read-length artifact.
* `7172_day4_with_selection` and `7172_day4_with_selection_repeat` are different
  populations (chr11L→chr11R fixed in one, chr6L→chr5R fixed in the other).
* the pipeline's silhouette clustering of Y' variants is coarser than the curated scheme
  (8 clusters vs 9–12 curated variants; the merged variants are 98.3–99.97 % identical);
  the curated library resolves more donors (see the pilot).

Unit tests: `python _pipeline/tests/test_attribution.py` and
`python _pipeline/tests/test_yprime_path.py` (or `pytest`).
Verification scripts and reports: `verification/`.
