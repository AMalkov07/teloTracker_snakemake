# Verification of the day-0 references and the recombination step

Everything here works on a **snapshot** of pipeline outputs copied from Argon
(`/nfsscratch/amalkova/telo_sra_runs`, mounted at `~/argon_scratch/telo_sra_runs`).
`verification/snapshot*/` are git-ignored; the scripts, reports and curated references are
committed.

| path | what |
|---|---|
| `snapshot/<ref>/pretelomeric_labels/` | day-0 labels: `*_simp.bed`, `extracted_yprimes_*.fasta`, probe/anchor/Y' BLAST tables, clustering dir |
| `snapshot/<ref>/run_config.yaml` | the run dir's `_pipeline/config.yaml` (library provenance) |
| `snapshot/<sample>/recombination/` | `*_features.tsv`, `*_recombination_summary.tsv` of the ORIGINAL runs (v1, `recombination_v1/` on Argon) |
| `snapshot/<sample>/<sample>_post_telo_trimming.tsv`, `_post_y_prime_probe.tsv`, `_read_summary.tsv` | per-read telomere / probe tables |
| `snapshot_argon_v2/`, `snapshot_argon_v2b/`, `snapshot_argon_v2c/` | the v2a / v2b / v2c Argon re-runs (see `docs/attribution_v2.md`, run history) |
| `snapshot_v2c_path/` | v2c with the (ID, ITS) path columns re-derived locally (`reprocess_snapshot.sh`) |
| `snapshot_curated*/` | the curated-library pilot (7302 day 0/4/5) |
| `curated_refs/<strain>_features/` | curated beds + Y' libraries (copied from `TeloTracker/references/`) |
| `spacer_lib_variable/<day0>/` | spacer-only chunk libraries (variable-size method) for 7302 / 7172 |
| `reports/partA/` | Part A: day-0 reference verification (one report per reference + `partA_summary.tsv`) |
| `reports/v1/`, `reports/argon_v2*/`, `reports/v2c_path/` | Part B: recombination checks per run version (`partB_report.md` + TSVs) |
| `reports/diff_*/` | per-read differences between run versions (`diff_runs.py`) |
| `reports/*/path/`, `reports/*/circles/`, `reports/*/circle_plots/` | path report, circle report, Y'-copy plots |
| `truth_sets/7302_circular_gains*.tsv` | reads at non-chr13L ends carrying the chr13L alternating array |
| `argon/` | the Argon re-run scripts (array scripts, Snakefile patch, manifests) |

## Part A — day-0 references vs the curated references

`run_partA.sh` runs `_pipeline/scripts/verify_day0_reference.py` for the 8 × 6991, 7172 and
7302 day-0 references against `curated_refs/<strain>_features/` (curated bed + curated Y'
library).

* **A1 per-end Y' counts**: exact for 7172 and 7302 and for 6 of 8 6991 references; the two
  known exceptions are `6991_day0_TeloTag` (chr4R 0 of 7) and `6991_day0_reference`
  (chr12R 6 of 7).
* **A2 feature coordinates** (anchor-relative): no major difference (> 500 bp) in 8 of 10
  references; only three recurring minor ones (a 25-bp chr5L ITS we do not label; the
  chr12R ITS 1-2 boundary 33 bp shorter; the chr13R anchor 6 bp longer).
  `6991_day0_with_selection` has a real defect at chr14L (Y'_1 936 bp short, the ITS 1-2
  absorbing ~1.8 kb) that the two repeat references do not show; `6991_day0_TeloTag` has
  the chr4R array missing.
* **A3 Y' grouping**: the pipeline's silhouette clustering is **coarser** than the curated
  scheme. Adjusted Rand Index vs the curated *variants* 0.54 (7302), 0.73 (7172), 0.54
  (6991); 8 clusters vs 9–12 curated variants. The curated variants that we merge are
  98.3–99.97 % identical to each other (e.g. `ID2_Red-Light` vs `ID2_Red-Dark` 99.96 %,
  `ID5_Blue` 99.97 %), i.e. a handful of SNPs over 6.6 kb — below what the clustering
  threshold (97 %) and the silhouette stop resolve. The `CONFLICT` elements are the curated
  `ID1_Gray` members, a label that spans both Long and Short elements (it behaves as an
  "unassigned singleton" bucket rather than a sequence group). The pre-clustering
  best-hit IDs (labeling BLASTs every Y' against the curated 6991 library) reproduce the
  curated assignment for 83–100 % of elements, so the reference-level variant calls are
  fine; the coarsening happens in the clustering step. Size class (Long/Short) agrees for
  every element.
* **A4 provenance**: every run used its own extracted per-strain library, no override, and
  BED ↔ library locations are a bijection.

## Part B — recombination

`run_partB.sh <snapshot> <report_dir>` runs `_pipeline/scripts/verify_recombination.py`:

* **B1 null** (day-0 self-runs): ~0.6–0.7 % of reads called recombinant with a non-Loss Y'
  status; per-end this exceeds 1 % at a handful of ends (max ~2.5 %). Y' Loss at the tandem
  ends (chr4R 11 %, chr12R 7 %, chr13L 4 %, chr14L 10 % in 7302) is mostly on confirmed
  telomere ends = sub-clonal copy-number heterogeneity, recorded as the per-end baseline.
* **B2 replicates**: the two 6991 day-0 replicate pairs agree within 5 points at 31 of 32
  ends. `7172_day4` vs `7172_day4_repeat` differ at ~10 ends because they are different
  populations (chr11L→chr11R 100 % vs 3 %; chr6L→chr5R 7 % vs 100 %).
* **B3 positive control** 7172 chr11L→chr11R: 97 / 93 / 96 % of reads at day 4 / 6 / 9 (v2).
* **B4 Y' Loss**: day-6 7302 chr13L Loss = 90 % confirmed-end real contraction; no end where
  the probe count contradicts RepeatMasker in more than 20 % of Loss reads.
* **B5 truth set** (chr13L alternating array at other ends): strict set 36 / 39 → chr13L
  with v2 (1 / 39 with v1); 31 / 32 with the curated library.

`diff_runs.py` compares two snapshots read-by-read; `reprocess_snapshot.sh` re-derives the
v2 columns (incl. the path) from any snapshot without cluster time.

## Argon re-runs (step 11+ only)

Scripts in `argon/`. Each array copies the staged code into every sample's run dir, patches
the run-dir Snakefile, keeps the previous outputs under `recombination_<tag>/`, deletes only
the step-11+ outputs and runs `snakemake all --rerun-triggers mtime` (essential: otherwise
steps 0–6 would be repeated). Stage with

```
cp _pipeline/scripts/{analyze_features,yprime_path,aggregate_recombination,recombination_utils}.py \
   verification/argon/patch_snakefile_v2.py  ~/argon_scratch/telo_sra_runs/v2_stage/path/
cp verification/argon/array_recomb_*.sh verification/argon/*_manifest.txt ~/argon_scratch/telo_sra_runs/v2_stage/
# on Argon (LANG=C for the JSV):
cd /nfsscratch/amalkova/telo_sra_runs && export LANG=C && qsub -v PREV_TAG=v2c v2_stage/array_recomb_rerun_v2c.sh
```

Run history: `recombination_v1/` (legacy), `recombination_v2a/` (v2 before the spacer-interval
fix), `recombination_v2b/` (spacer-interval fix, 50-kb library), `recombination/` (v2c: spacer-only
library + path). The curated-library pilot writes to `results/<sample>__curatedYP/`.
