# TeloTracker Snakemake Pipeline

A Snakemake-based pipeline for analyzing telomeric regions in yeast Oxford Nanopore sequencing data. It identifies chromosome anchors, trims reads to telomeric regions, detects Y′ elements, X elements, and spacers, and tracks recombination events across time points.

The pipeline is driven by a single Python wrapper, `run_pipeline.py`, so you don't need to hand-edit shell scripts or `config.yaml` every time a parameter changes.

---

## Pipeline Overview

The pipeline has three steps that map onto the three `--steps` choices of `run_pipeline.py`:

| Step | What it does | Script(s) involved |
|---|---|---|
| `create_ref` | Runs the Snakemake pipeline through `through_y_prime_analysis`, then assembles and polishes a strain-specific telomeric reference using Flye and Dorado. | `create_ref.sh` (which internally runs `snakemake`) |
| `label_regions` | Labels pre-telomeric features (Y′ elements, X elements, spacers) on the new reference and performs automatic Y′ clustering. | `label_regions.sh` |
| `recombination` | For future time points, runs `snakemake recombination_summary` against the day-0 reference to detect recombination events. | `Snakefile` (direct Snakemake invocation, one call per time point) |

Typical workflow for a new strain: run `create_ref` + `label_regions` **once** on the day-0 sample, then run `recombination` for each later time point.

---

## Setup

### 1. Conda environment

```bash
conda env create -f consensus.yaml
conda activate consensus
```

Installs everything needed (BLAST, minimap2, samtools, Flye, Snakemake, RepeatMasker, PyYAML, Python deps).

### 2. Dorado Singularity image (not in repo)

`create_ref.sh` polishes with Dorado. You need a Singularity `.sif` image on disk:

```bash
singularity pull dorado.sif docker://nanoporetech/dorado:latest
```

Point at it with `dorado_image` in the config.

### 3. Input data

- Place your BAM (or FASTQ / FASTQ.GZ) in `samples_dorado_basecalled/` — the filename's stem must match `base_name` in the config.
- The `references/` directory is pre-populated with anchor FASTAs, the Y′ probe, adapter file, and strain-specific pairing libraries.

### 4. (Cluster) SGE job headers

`create_ref.sh` and `label_regions.sh` have SGE headers at the top (`#$ -q UI-GPU`, `#$ -pe smp 56`, etc.) so they can be submitted with `qsub`. See [Running on the cluster](#running-on-the-cluster) below.

---

## Quick Start

```bash
# 1. Generate a sample config to edit:
python run_pipeline.py --init-config

# 2. Edit pipeline_config.yaml (base_name, strain, dorado_image, samples, etc.)

# 3. Build the day-0 reference (default: create_ref + label_regions):
python run_pipeline.py --config pipeline_config.yaml

# 4. Later, run recombination analysis on future time points:
python run_pipeline.py --config pipeline_config.yaml --steps recombination

# Or do everything in one shot:
python run_pipeline.py --config pipeline_config.yaml --steps all

# Preview what would happen without executing:
python run_pipeline.py --config pipeline_config.yaml --dry-run
```

---

## `run_pipeline.py` Commands

### Core flags

| Flag | Purpose |
|---|---|
| `--config FILE` | Load settings from a YAML config file. |
| `--init-config [FILE]` | Write a fresh `pipeline_config.yaml` template and exit. Default filename: `pipeline_config.yaml`. |
| `--steps STEP ...` | Which step(s) to run. Choices: `create_ref`, `label_regions`, `recombination`, `day0`, `all`. Default: `day0` (= `create_ref` + `label_regions`). `all` = `day0` + `recombination`. |
| `--dry-run` | Print what would be executed — patched scripts, generated `config.yaml`, Snakemake command — without running anything. Use this first. |
| `--no-organize` | Skip the post-step pass that creates sibling symlinks (`reference.fasta`, `labels.bed`, …) at the top of `results/<base_name>/` pointing into `_pipeline/`. See [Output layout](#output-layout). |

### Per-parameter overrides (win over the config file)

| Flag | What it sets |
|---|---|
| `--base-name NAME` | Day-0 sample base name (no file extension, no `.bam`). |
| `--strain ID` | Yeast strain number (e.g. `7302`). |
| `--bam-dir DIR` | Directory holding input BAM/FASTQ. Default `samples_dorado_basecalled`. |
| `--threads N` | Core count. Default 56. |
| `--anchor-set NAME` | Anchor FASTA set (default `telomerase_shutoff_anchors`). |
| `--dorado-mode {docker,local}` | Docker (Singularity) or local PATH dorado. |
| `--dorado-image PATH` | `.sif` path (only needed when mode is `docker`). |
| `--dorado-model MODEL` | Basecaller model (e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`). |
| `--reference-fasta PATH` | Custom reference for `label_regions` (defaults to `create_ref`'s output). |
| `--samples NAME ...` | Time-point samples for `recombination`. Overrides `timepoint_samples` in the config. |

---

## Configuration Reference

`pipeline_config.yaml` — one file, all settings.

```yaml
# Required
base_name: "dorado_7302_day0_PromethION_no_tag_yes_rejection"
strain: "7302"
bam_dir: "samples_dorado_basecalled"

# Dorado (required for create_ref)
dorado_mode: "docker"                # or "local"
dorado_image: "/path/to/dorado.sif"  # only needed for docker mode
dorado_model: "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"

# Optional — shown with their defaults
threads: 56
anchor_set: "telomerase_shutoff_anchors"
min_raw_gapped_score: 5000
yprime_identity_threshold: 97.0

# Optional — override the auto-extracted Y′ library with a curated FASTA
# y_prime_lib_override: "references/7302_features/repeatmasker_7302_all_y_primes.fasta"

# Optional — custom reference for label_regions
# reference_fasta: "path/to/custom_reference.fasta"

# Recombination — day-0 reference and time-point samples
# day0_base_name defaults to base_name if omitted
# day0_base_name: "dorado_7302_day0_PromethION_no_tag_yes_rejection"

timepoint_samples:
  - "dorado_7302_day3_PromethION_no_tag_yes_rejection"
  - "dorado_7302_day6_PromethION_no_tag_yes_rejection"
```

---

## Output layout

Every pipeline run writes into `results/<base_name>/`. To keep the top of that directory readable, the pipeline puts all its real output inside a `_pipeline/` subfolder and creates sibling **symlinks** at the top for the handful of files users actually care about:

```
results/<base_name>/
├── _pipeline/                       ← every file produced by Snakemake + the shell scripts
│   ├── assembly_<strain>/           (create_ref.sh output)
│   ├── pretelomeric_labels/         (label_regions.sh output)
│   ├── recombination/               (snakemake recombination_summary output)
│   ├── blast/, porechop/, graphs/   (intermediates — feel free to ignore)
│   └── ...                          (everything else the pipeline produces)
├── reference.fasta                  → _pipeline/assembly_<strain>/assembly_<strain>_dorado_reference.fasta
├── labels.bed                       → _pipeline/pretelomeric_labels/pretelomeric_regions_<strain>_simp.bed
├── labels.gff3                      → _pipeline/pretelomeric_labels/pretelomeric_regions_<strain>.gff3
├── quality_report.txt               → _pipeline/pretelomeric_labels/pretelomeric_regions_<strain>_quality_report.txt
├── recombination_summary.tsv        → _pipeline/recombination/<base_name>_recombination_summary.tsv
└── MANIFEST.md                      (short note about what's where)
```

Important:

- **Deleting a symlink does not delete the underlying file.** The symlinks are just pointers; the real files live inside `_pipeline/`.
- **A symlink only appears after its target has been produced.** If `create_ref.sh` has finished but `label_regions.sh` hasn't, you'll see `reference.fasta` but not `labels.bed`.
- **`scripts/organize_outputs.py`** creates and refreshes these symlinks. It's run automatically after every step by `run_pipeline.py`; disable with `--no-organize` if you'd rather work with the raw `_pipeline/` layout.
- **Migration from older runs:** if you have existing `results/<base>/` directories produced before this change, move them by hand:
  ```bash
  cd results/<base_name>
  mkdir _pipeline
  # shopt -s dotglob  # uncomment if you also want hidden files moved
  mv * _pipeline/ 2>/dev/null
  python ../../scripts/organize_outputs.py <base_name>
  ```

---

## Running on the cluster

The user's workflow is to run `run_pipeline.py` **on the cluster** (not over SSH from a laptop). Two options:

### Option 1: Interactive node

```bash
ssh cluster
conda activate consensus
cd /path/to/teloTracker_snakemake
python run_pipeline.py --config pipeline_config.yaml
```

Fine for small/quick jobs. The Snakemake pipeline still parallelizes across `--threads` cores, but the whole thing runs as your interactive session.

### Option 2: Submit the wrapper itself via `qsub` (recommended for long jobs)

Wrap the wrapper in a tiny SGE script:

```bash
#!/bin/bash
#$ -q UI-GPU
#$ -pe smp 56
#$ -l gpu=true
#$ -j y
#$ -cwd

# qsub spawns a non-interactive shell that does NOT source .bashrc, so
# `conda activate` is not defined until we source conda.sh explicitly.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus

python run_pipeline.py --config pipeline_config.yaml
```

Then `qsub run_wrapper.sh`. The wrapper's `subprocess.run("bash create_ref.sh")` calls inherit the SGE allocation.

Note: `run_pipeline.py` itself does **not** call `qsub` — it runs the shell scripts with `bash` directly via `subprocess`. The SGE headers on `create_ref.sh` / `label_regions.sh` only matter if you submit **those files** (not the wrapper) as jobs.

---

## Under the hood — what `run_pipeline.py` actually changes

### `config.yaml` (written, not patched)

Each run, the wrapper overwrites `config.yaml` at the repo root with a freshly generated file based on your `pipeline_config.yaml` values. The keys it emits are the ones the Snakefile reads at module load (`base_name`, `day0_base_name`, `anchor_set`, `bam_dir`, `strain`, `min_raw_gapped_score`, `yprime_identity_threshold`, `references.*`).

For the `recombination` step with multiple `timepoint_samples`, `config.yaml` is **rewritten between samples** so each `snakemake recombination_summary` call targets the right time point while keeping `day0_base_name` pinned to the reference.

### `create_ref.sh` and `label_regions.sh` (patched in memory)

The originals on disk are never modified. For each run:

1. The wrapper reads the script file into memory.
2. It does regex substitutions on bash variable assignments:
   - `BASE_NAME=...` → `BASE_NAME="<your value>"`
   - `STRAIN_ID=...`, `THREADS=...`, `DORADO_MODE=...`, `DORADO_IMAGE=...`, `DORADO_MODEL=...` (for `create_ref`)
   - `REFERENCE_DIR=...`, `REFERENCE_FASTA=...` (for `label_regions`, pointing at `create_ref`'s output or your `--reference-fasta`)
3. The patched content is written to a temporary file in the repo root (`telo_create_ref_*.sh` or `telo_label_regions_*.sh`).
4. `bash <tempfile>` is executed with the repo root as `cwd`.
5. The temp file is deleted on successful exit; on failure the path is printed so you can inspect.

The regex (in `substitute_var`) handles quoted values, unquoted values, and trailing comments. It matches line-start, so commented-out assignments (`#VAR=...`) are left alone.

### `snakemake recombination_summary` (direct call)

For the `recombination` step the wrapper doesn't touch a shell script — it just calls `snakemake recombination_summary -c <threads>` directly, once per time-point sample, rewriting `config.yaml` between calls.

### Files that are NEVER modified

- `create_ref.sh` (on disk) — patches happen in memory only.
- `label_regions.sh` (on disk) — same.
- `Snakefile` — never patched; all Snakemake customization goes through `config.yaml`.
- Anything in `references/`, `scripts/`, or `samples_dorado_basecalled/`.

---

## Troubleshooting

- **`CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'`** — qsub (and subprocess bash) spawn non-interactive shells that don't source `.bashrc`, so the `conda` shell function isn't loaded. Fix by sourcing `conda.sh` in your qsub wrapper:
  ```bash
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate consensus
  ```
  `run_pipeline.py` already injects the same `source` at the top of patched `create_ref.sh`/`label_regions.sh`, so you only need this in the outer qsub wrapper.
- **"PyYAML is required"** — activate the `consensus` conda env first.
- **"could not find 'VAR=' in ... to substitute"** — a variable name in the wrapper no longer matches the shell script (e.g. you renamed a bash var). The wrapper will warn but continue; the script will run with its default value.
- **Snakemake fails at module load with `KeyError: ...`** — your `pipeline_config.yaml` is missing a key the Snakefile expects. Run with `--dry-run`, inspect the generated `config.yaml`, and verify it has `references.y_prime_lib`, `spacer_lib_dir`, `x_element_lib_dir`, etc.
- **Broken `config.yaml` after a `recombination` run** — the wrapper overwrites `config.yaml` with the **last** sample's settings. Re-run `--steps create_ref` or `--steps label_regions` (both call the config writer) to reset it for the day-0 sample.
- **Want to see what would happen without risk** — `--dry-run` prints the patched scripts and the Snakemake commands without executing.

---

## Post-processing — per-read structure diagrams

After `recombination` has produced `*_features.tsv`, you can render a schematic of any individual read's pretelomeric structure with `_pipeline/scripts/plot_read_structure.py`. This is the easiest way to sanity-check a recombination call by eye.

### Running it

The script has two modes: **list** the reads available for plotting, then **plot** one of them. Both require either `--features <one_features.tsv>` or `--recombination-dir <dir> --base-name <base>` so it knows where to look.

```bash
# 1. List reads that look interesting (e.g. multi-Y' recombination events)
python _pipeline/scripts/plot_read_structure.py \
    --recombination-dir results/<base_name>/_pipeline/recombination \
    --base-name <base_name> \
    --list-reads --min-yprimes 4 --recombination-only

# 2. Plot one
python _pipeline/scripts/plot_read_structure.py \
    --recombination-dir results/<base_name>/_pipeline/recombination \
    --base-name <base_name> \
    --read-id <read_id_from_list> \
    --output diagram.png
```

Optional filters for `--list-reads`: `--chr-end chr10R` (single chromosome end), `--recombination-only` (skip `recombination_detected = False`), `--min-yprimes N`.

#### Y′ colors

The script auto-detects the Y′ library at `results/*/_pipeline/pretelomeric_labels/extracted_yprimes_*.fasta` and reads the per-ID color names from the FASTA headers. If you have several runs with different libraries, pin the right one explicitly:

```bash
--yprime-lib results/<day0_base_name>/_pipeline/pretelomeric_labels/extracted_yprimes_<strain>.fasta
```

### Reading the diagram

Reads are always drawn **telomere on the outer end → chromosome on the chr-label end**. For an L-arm read (`chr_end` ending in `L`), the telomere is on the left and the chromosome label is on the right; R-arms are mirrored. The horizontal rail represents the read backbone, drawn from the telomere outer edge to the chromosome label.

| Feature | Shape | Color | Meaning |
|---|---|---|---|
| **Telomere** | Square box at the outer end | Black (same as inter-Y′ spacer) | TG/AC repeat tract. Tract length in bp printed below the box (sourced from `<base>_post_telo_trimming.tsv`). |
| **Y′ element** | Pentagon arrow, pointing toward the telomere | Per-ID color, from the Y′ library | A Y′ element. Label above is positional (`Y-1` is innermost / chromosome-side, `Y-N` is outermost / telomere-side). Label below is the Y′ ID (`ID1`, `ID2`, …) from clustering. |
| **inter-Y′ spacer** | Small square between Y′ arrows | Black | Gap between two consecutive Y′ elements on the read. The number below is the gap length in bp. |
| **X element** | Dark grey square | Dark grey | X element. Label above shows the reference chromosome end it matched (e.g. `X chr12R`); `bp` underneath is its size. |
| **Spacer** | Blue-grey box | `#7399AB` | Pretelomeric spacer between the X element and the anchor. Label above shows the reference chromosome end it matched. |
| **Anchor** | Dark `anc` box | Black | Subtelomeric anchor sequence. Chromosome assignment comes from which anchor matched (printed as the chromosome label). |
| **★ + red outline** | On a Y′ arrow | — | Divergence point. The first Y′ position (innermost → outermost) where the read disagrees with the day-0 reference. The title bar reports `expected_at_divergence` and `found_at_divergence`. |

### Reading the `Y' recombination status`

The title bar always shows three things: the read ID, then `<chr_end>  ·  <N> Y' elements  ·  <status>`, and optionally `source: <chr_end>` (the reference chromosome end that the X element and/or spacer best match — telling you which donor end the recombined material likely came from).

Possible `<status>` values, produced by [`analyze_features.py`](_pipeline/scripts/analyze_features.py):

| Status | Meaning | What the diagram looks like |
|---|---|---|
| `No Change` | The Y′ array on the read matches the day-0 reference 1:1, position by position. | No ★ marker; Y′ colors should be the same as the reference for this chromosome end. |
| `1st Y' Change` | The **innermost** Y′ (`Y-1`, adjacent to the X element) differs from the reference. This is the most common ALT-recombination signal — the read still has the right chromosome end, but its first Y′ is new. | ★ on `Y-1` (and a red outline). `expected_at_divergence` = reference's `Y-1` ID; `found_at_divergence` = the new ID actually on the read. |
| `Y' Recombination` | Same as above but the divergence is at `Y-2` or further out — the innermost Y′ still matches reference but a downstream one differs. Rarer. | ★ on whichever `Y-k` is the divergent one. |
| `Y' Gain` | The read has **more Y′ elements** than the reference for this chromosome end (extension events). | ★ on the first "extra" Y′ beyond the reference's count. `expected_at_divergence` = `None`. |
| `Y' Loss` | The read has **fewer Y′ elements** than the reference (contraction events, or reference has Y′ but read has none). | If any Y′ are present, ★ is on the position where the array runs out. `found_at_divergence` = `None`. |

Two additional `*_features.tsv` columns help interpret the status — they're not directly drawn but are visible in `--list-reads` output:

- `recombination_source` — the reference chromosome end that the X-element + spacer best match. If `source: chrXY` differs from `chr_end: chrZW`, the read's chromosome end has likely been replaced by material from a different donor end. `ambiguous` means the X/spacer didn't pick a single donor confidently.
- `y_prime_compatible_ends` — chromosome ends whose `Y-1` ID matches the read's `Y-1` ID. Useful for working out where a `1st Y' Change` event's new Y′ came from.

### Worked example

A read whose title reads:

```
SRR33298393.133768
chr12R  ·  10 Y' elements  ·  Y' Gain
source: chr4
```

…is an R-arm read assigned to chr12R, with 10 Y′ elements on it. The day-0 chr12R reference has fewer than 10 Y′, so the call is `Y' Gain`. The X-element + spacer match best to chr4, suggesting the read's pretelomere is the recombined product of a chr12R anchor stitched to chr4-derived pretelomeric material that has been elongated. The ★ on the diagram will sit on the first "extra" Y′ position beyond the original chr12R array.

---

## Repository layout

```
.
├── run_pipeline.py            # Single-entry wrapper
├── pipeline_config.yaml       # User-facing config
├── config.yaml                # Auto-generated by run_pipeline.py (do not hand-edit)
├── Snakefile                  # Full snakemake pipeline
├── create_ref.sh              # Reference assembly + polishing
├── label_regions.sh           # Pre-telomeric feature labeling
├── consensus.yaml             # Conda environment spec
├── references/                # Anchors, probe, adapters, strain pairing libraries
├── scripts/                   # Python scripts called by Snakefile / label_regions
│   └── organize_outputs.py    # Creates sibling symlinks under results/<base>/
├── samples_dorado_basecalled/ # Your input BAM/FASTQ files
└── results/                   # Per-sample outputs — see "Output layout" above
    └── <base_name>/
        ├── _pipeline/         # Everything produced by the pipeline
        ├── reference.fasta    # Symlinks to the user-facing finals
        ├── labels.bed
        ├── labels.gff3
        ├── quality_report.txt
        ├── recombination_summary.tsv
        └── MANIFEST.md
```
