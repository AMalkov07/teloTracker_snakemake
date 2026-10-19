# Day-0 reference genomes

Curated day-0 genomes for the three strains analysed with TeloTracker, with the chromosome-end labels made from
them. They are the references behind the cut99.5 recombination analyses (code `8be3cfa`, 2026-09-28) and the
published results built on them.

| Folder | Genome | Labels |
|---|---|---|
| `6991_curated/` | `6991_curated.fasta` (16 chromosomes) | `pretelomeric_regions_6991_curated_simp.bed` |
| `7172_curated/` | `7172_curated.fasta` (16 chromosomes) | `pretelomeric_regions_7172_curated_simp.bed` |
| `7302_curated/` | `7302_curated.fasta` (16 chromosomes) | `pretelomeric_regions_7302_curated_simp.bed` |

## How the labels were made

`label_regions` on each genome, run through `run_pipeline.py` with:

```yaml
base_name: "<strain>_curated"
strain: "<strain>_curated"
reference_fasta: <this genome>
anchor_set: "telomerase_shutoff_anchors"
yprime_stop_mode: "threshold"      # Y' groups at a fixed 99.5 % identity cutoff
yprime_identity_threshold: 99.5
yprime_linkage: "average"
```

Y' group IDs (ID1, ID2, ...) depend on this grouping; the same genome grouped at 97 % or by silhouette gives
different IDs.

## Coordinates

The `_simp.bed` files are written **1-based** (start + 1), not standard 0-based BED. A planned fix moves the
label output to standard BED; these files will be regenerated and re-committed then.

## Using one as the day-0 reference

```yaml
base_name: "6991_curated"
strain: "6991_curated"
reference_fasta: "_pipeline/references/day0/6991_curated/6991_curated.fasta"
yprime_grouping: 99.5              # same as yprime_stop_mode: threshold + yprime_identity_threshold: 99.5
```

Run `--steps label_regions recombination` (see the main README). `label_regions` rebuilds the Y' / spacer / X
libraries recombination needs under `results/<base_name>/_pipeline/`.

## Provenance

The genomes are identical (MD5) to `/Shared/malkova_lab/Andrey/teloTracker_runs/cut995/reference/<strain>.fasta`
and to the `reference_fasta` of each cut99.5 run; the BEDs are identical to those runs' label output.

| File | MD5 |
|---|---|
| `6991_curated.fasta` | `a6a31115b59e709119ef7a80e16924a1` |
| `7172_curated.fasta` | `57c6b36dddcaa7c826875f821d949506` |
| `7302_curated.fasta` | `13999da08608bc7bf6da4b18b4743d38` |
| `pretelomeric_regions_6991_curated_simp.bed` | `279286a9d747c29a598df310f9b8da2b` |
| `pretelomeric_regions_7172_curated_simp.bed` | `e93bbedbee8bcc663cf6799fbb4b78ba` |
| `pretelomeric_regions_7302_curated_simp.bed` | `6e92130bbfa37e2a1277e8253fe3a647` |

## Other references in this repository

- `results/dorado_7302_day0_PromethION_no_tag_yes_rejection/` holds an older **create_ref build** of 7302 (added
  March 2026, `637cc9b`). It predates ends-only polishing and is **not** the curated 7302 genome above. It is
  in the old folder layout (`results/<sample>/assembly_7302/`, no `_pipeline/`), so the current pipeline does
  not find it on its own; it stays because the README and config template use that sample name as their example.
- Files in `_pipeline/references/` that the pipeline reads by default: the anchor sets
  (`telomerase_shutoff_anchors.fasta`, `test_anchors.fasta`, with BLAST databases), the create_ref base reference
  `6991_only_to_anchors.fasta`, the Y' probe `y_prime_probe.fasta`, the adapter sequences, and the labelling
  libraries `repeatmasker_6991_all_y_primes.fasta` and `6991_xprimes.fasta`.
- `pairings_for_spacers/`, `pairings_for_x_element_ends/`, `extracted_yprimes_7302.fasta` and
  `7302_features/` are libraries from an earlier layout. Current runs build these per day-0 sample under
  `results/<day0>/_pipeline/pretelomeric_labels/`; the pipeline no longer reads them.
