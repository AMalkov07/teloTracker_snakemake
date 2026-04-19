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

- **"PyYAML is required"** — activate the `consensus` conda env first.
- **"could not find 'VAR=' in ... to substitute"** — a variable name in the wrapper no longer matches the shell script (e.g. you renamed a bash var). The wrapper will warn but continue; the script will run with its default value.
- **Snakemake fails at module load with `KeyError: ...`** — your `pipeline_config.yaml` is missing a key the Snakefile expects. Run with `--dry-run`, inspect the generated `config.yaml`, and verify it has `references.y_prime_lib`, `spacer_lib_dir`, `x_element_lib_dir`, etc.
- **Broken `config.yaml` after a `recombination` run** — the wrapper overwrites `config.yaml` with the **last** sample's settings. Re-run `--steps create_ref` or `--steps label_regions` (both call the config writer) to reset it for the day-0 sample.
- **Want to see what would happen without risk** — `--dry-run` prints the patched scripts and the Snakemake commands without executing.

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
└── samples_dorado_basecalled/ # Your input BAM/FASTQ files
```
