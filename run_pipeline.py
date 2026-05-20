#!/usr/bin/env python3
"""
TeloTracker Pipeline Wrapper
============================
Eliminates the need to manually edit config.yaml, create_ref.sh, and
label_regions.sh separately.

Configure once via CLI args or a single YAML config file, then run
any step (or all steps) in sequence.

Usage examples
--------------
  # Generate a sample config file to fill in:
  python run_pipeline.py --init-config

  # Default: build the day-0 reference (create_ref -> label_regions):
  python run_pipeline.py --config pipeline_config.yaml

  # Run recombination analysis on later time points:
  python run_pipeline.py --config pipeline_config.yaml --steps recombination

  # Full pipeline in one shot (day0 + recombination):
  python run_pipeline.py --config pipeline_config.yaml --steps all

  # Run only a single step:
  python run_pipeline.py --config pipeline_config.yaml --steps create_ref

  # Override any value on the command line:
  python run_pipeline.py --config pipeline_config.yaml --threads 32 --steps create_ref

  # Fully CLI-driven (no config file needed for basic use):
  python run_pipeline.py \\
      --base-name dorado_7302_day0_PromethION_no_tag_yes_rejection \\
      --strain 7302 \\
      --dorado-mode local

  # Dry run to preview substitutions without executing:
  python run_pipeline.py --config pipeline_config.yaml --dry-run
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import yaml
    _YAML_AVAILABLE = True
except ImportError:
    _YAML_AVAILABLE = False

# Root directory of the pipeline (same dir as this script)
PIPELINE_DIR = Path(__file__).resolve().parent

# ──────────────────────────────────────────────────────────────────────────────
# Config helpers
# ──────────────────────────────────────────────────────────────────────────────

def load_yaml_config(config_path: str) -> dict:
    if not _YAML_AVAILABLE:
        print("ERROR: PyYAML is required to use --config.")
        print("       Install with: pip install pyyaml")
        sys.exit(1)
    with open(config_path) as f:
        return yaml.safe_load(f) or {}


def merge_configs(file_cfg: dict, args) -> dict:
    """Merge file config with CLI args; CLI args take priority."""
    cfg = dict(file_cfg)

    cli_map = {
        "base_name":      args.base_name,
        "strain":         args.strain,
        "bam_dir":        args.bam_dir,
        "threads":        args.threads,
        "anchor_set":     args.anchor_set,
        "dorado_mode":    args.dorado_mode,
        "dorado_image":   args.dorado_image,
        "dorado_model":   args.dorado_model,
        "reference_fasta": args.reference_fasta,
    }
    for key, val in cli_map.items():
        if val is not None:
            cfg[key] = val

    # Samples for the analyze step
    if args.samples:
        cfg["timepoint_samples"] = args.samples

    return cfg


def validate_cfg(cfg: dict, steps: list):
    errors = []

    needs_base = any(s in steps for s in ("create_ref", "label_regions"))
    if needs_base:
        if not cfg.get("base_name"):
            errors.append("--base-name (or base_name in config) is required for create_ref / label_regions")
        if not cfg.get("strain"):
            errors.append("--strain (or strain in config) is required for create_ref / label_regions")

    if "create_ref" in steps:
        mode = cfg.get("dorado_mode", "docker")
        if mode == "docker" and not cfg.get("dorado_image"):
            print("WARNING: dorado_mode is 'docker' but no dorado_image specified. "
                  "create_ref.sh will use its own hardcoded DORADO_IMAGE value.")

    if "recombination" in steps:
        if not cfg.get("timepoint_samples") and not cfg.get("base_name"):
            errors.append("At least one sample is required for the recombination step. "
                          "Use --samples or set timepoint_samples in the config file.")

    if errors:
        for e in errors:
            print(f"ERROR: {e}")
        sys.exit(1)


# ──────────────────────────────────────────────────────────────────────────────
# Shell script patching
# ──────────────────────────────────────────────────────────────────────────────

def substitute_var(content: str, var_name: str, new_value: str, script_name: str = "") -> str:
    """
    Replace a bash variable assignment line like:
        VAR_NAME="old"           → VAR_NAME="new_value"
        VAR_NAME=old             → VAR_NAME="new_value"
        VAR_NAME="old"  # comment  → VAR_NAME="new_value"  # comment

    Handles quoted values, unquoted values, and preserves trailing comments.
    Anchored to line-start so commented-out lines (#VAR_NAME=...) are ignored.
    """
    # Group 1: "VAR="
    # Group 2: the current value — either a double-quoted string or a bare word
    # Group 3: optional trailing whitespace + comment
    pattern = rf'^({re.escape(var_name)}=)("(?:[^"\\]|\\.)*"|\S*)(\s*#.*)?$'
    replacement = rf'\g<1>"{new_value}"\g<3>'
    new_content, n = re.subn(pattern, replacement, content, flags=re.MULTILINE)
    if n == 0:
        print(f"  WARNING: could not find '{var_name}=' in {script_name} to substitute.")
    return new_content


def patch_script(script_path: Path, subs: dict) -> str:
    """Apply multiple variable substitutions to a script file."""
    content = script_path.read_text()
    for var, value in subs.items():
        content = substitute_var(content, var, str(value), script_path.name)
    return content


# ──────────────────────────────────────────────────────────────────────────────
# Snakemake config.yaml writer
# ──────────────────────────────────────────────────────────────────────────────

def write_snakemake_config(cfg: dict, dest: Path = None):
    """Write the config.yaml consumed by the Snakemake pipeline.

    Emits the key structure the current Snakefile reads at module load time.
    """
    if dest is None:
        dest = PIPELINE_DIR / "_pipeline" / "config.yaml"

    base_name      = cfg["base_name"]
    strain         = cfg["strain"]
    bam_dir        = cfg.get("bam_dir", "samples_dorado_basecalled")
    anchor_set     = cfg.get("anchor_set", "telomerase_shutoff_anchors")
    day0_base_name = cfg.get("day0_base_name", base_name)
    min_gapped     = cfg.get("min_raw_gapped_score", 5000)
    yp_threshold   = cfg.get("yprime_identity_threshold", 97.0)
    yp_override    = cfg.get("y_prime_lib_override")

    lines = [
        "# config.yaml — auto-generated by run_pipeline.py",
        f'base_name: "{base_name}"',
        f'day0_base_name: "{day0_base_name}"',
        f'anchor_set: "{anchor_set}"',
        f'bam_dir: "{bam_dir}"',
        f'strain: "{strain}"',
        "",
        "# BLAST parameter",
        f"min_raw_gapped_score: {min_gapped}",
        "",
        "# Y prime clustering threshold (percent identity for grouping)",
        f"yprime_identity_threshold: {yp_threshold}",
        "",
        "# Paths to references",
        "references:",
        '  anchors: "_pipeline/references/telomerase_shutoff_anchors.fasta"',
        '  adapters: "_pipeline/references/nanopore_sqk-slk114_adapter_sequence_truncated.txt"',
        '  probe: "_pipeline/references/y_prime_probe.fasta"',
        "",
        f'  y_prime_lib: "results/{day0_base_name}/_pipeline/pretelomeric_labels/extracted_yprimes_{{strain}}.fasta"',
    ]
    if yp_override:
        lines.append(f'  y_prime_lib_override: "{yp_override}"')
    lines.extend([
        f'  spacer_lib_dir: "results/{day0_base_name}/_pipeline/pretelomeric_labels/pairings_for_spacers/{{strain}}_pairings"',
        f'  x_element_lib: "results/{day0_base_name}/_pipeline/pretelomeric_labels/clustered_x_elements_{{strain}}.fasta"',
        "",
    ])
    dest.write_text("\n".join(lines))
    print(f"  Written: {dest}")


# ──────────────────────────────────────────────────────────────────────────────
# Script runner
# ──────────────────────────────────────────────────────────────────────────────

# Prepended to every patched shell script before execution so that
# `conda activate` works even though subprocess bash is non-interactive
# (shell functions aren't inherited from the parent — env vars are,
# so `CONDA_SHLVL` alone is NOT a reliable "conda is ready" signal).
# Check for the `conda` shell function itself. Runs above any `set -u`
# / `set -e` so a missing conda.sh falls through silently rather than
# aborting before we see the real error.
CONDA_INIT_PROLOGUE = (
    "# Auto-injected by run_pipeline.py — initialize conda for non-interactive shells\n"
    "if ! declare -F conda > /dev/null 2>&1; then\n"
    '    _conda_base="$(conda info --base 2>/dev/null)"\n'
    '    if [ -n "$_conda_base" ] && [ -f "$_conda_base/etc/profile.d/conda.sh" ]; then\n'
    "        # shellcheck disable=SC1091\n"
    '        . "$_conda_base/etc/profile.d/conda.sh"\n'
    "    fi\n"
    "fi\n"
    "\n"
)


def run_script(patched_content: str, label: str, extra_args: list = None, dry_run: bool = False):
    """Write patched content to a temp file and execute it with bash."""
    extra_args = extra_args or []
    cmd_display = f"bash {label} {' '.join(extra_args)}".strip()

    patched_content = CONDA_INIT_PROLOGUE + patched_content

    if dry_run:
        print(f"\n  [DRY RUN] Would execute: {cmd_display}")
        print("  --- patched script (first 60 lines) ---")
        for i, line in enumerate(patched_content.splitlines()[:60], 1):
            print(f"  {i:3d}  {line}")
        if patched_content.count("\n") > 60:
            print("  ... (truncated)")
        print("  --- end ---")
        return

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sh", delete=False,
        prefix=f"telo_{label.replace('.sh','')}_", dir=PIPELINE_DIR
    ) as tmp:
        tmp.write(patched_content)
        tmp_path = tmp.name

    os.chmod(tmp_path, 0o755)
    print(f"\n  Executing: {cmd_display}")
    print(f"  (temp script: {tmp_path})")

    try:
        subprocess.run(["bash", tmp_path] + extra_args, check=True, cwd=PIPELINE_DIR)
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: {label} failed (exit code {e.returncode})")
        print(f"  Temp script retained at: {tmp_path}")
        sys.exit(e.returncode)
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


# ──────────────────────────────────────────────────────────────────────────────
# Post-step organization
# ──────────────────────────────────────────────────────────────────────────────

def organize_outputs(base_name: str, strain: str, enabled: bool = True, dry_run: bool = False):
    """Run _pipeline/scripts/organize_outputs.py to (re)create sibling symlinks + MANIFEST."""
    if not enabled or dry_run:
        return
    script = PIPELINE_DIR / "_pipeline" / "scripts" / "organize_outputs.py"
    cmd = [sys.executable, str(script), base_name, "--strain", strain, "--quiet"]
    try:
        subprocess.run(cmd, check=True, cwd=PIPELINE_DIR)
    except subprocess.CalledProcessError as e:
        print(f"  WARNING: organize_outputs exited with code {e.returncode} — symlinks may be incomplete")


# ──────────────────────────────────────────────────────────────────────────────
# Pipeline steps
# ──────────────────────────────────────────────────────────────────────────────

def step_create_ref(cfg: dict, dry_run: bool = False):
    print("\n" + "=" * 70)
    print("STEP 1: Create Reference  (create_ref.sh)")
    print("=" * 70)

    base_name   = cfg["base_name"]
    strain      = cfg["strain"]
    threads     = cfg.get("threads", 56)
    dorado_mode = cfg.get("dorado_mode", "docker")
    dorado_image = cfg.get("dorado_image", "")
    dorado_model = cfg.get("dorado_model", "dna_r10.4.1_e8.2_400bps_sup@v5.2.0")

    # create_ref.sh expects BASE_NAME without extension
    if base_name.endswith(".bam"):
        base_name = base_name[:-4]

    print(f"  base_name:    {base_name}")
    print(f"  strain:       {strain}")
    print(f"  threads:      {threads}")
    print(f"  dorado_mode:  {dorado_mode}")
    if dorado_image:
        print(f"  dorado_image: {dorado_image}")
    print(f"  dorado_model: {dorado_model}")

    # 1. Write Snakemake config.yaml
    print("\n  Writing config.yaml for Snakemake...")
    if not dry_run:
        write_snakemake_config(cfg)

    # 2. Patch create_ref.sh
    subs = {
        "BASE_NAME":    base_name,
        "STRAIN_ID":    strain,
        "THREADS":      str(threads),
        "DORADO_MODE":  dorado_mode,
        "DORADO_MODEL": dorado_model,
    }
    if dorado_image:
        subs["DORADO_IMAGE"] = dorado_image

    patched = patch_script(PIPELINE_DIR / "_pipeline" / "create_ref.sh", subs)
    run_script(patched, "create_ref.sh", dry_run=dry_run)


def step_label_regions(cfg: dict, dry_run: bool = False):
    print("\n" + "=" * 70)
    print("STEP 2: Label Pre-Telomeric Regions  (label_regions.sh)")
    print("=" * 70)

    base_name = cfg["base_name"]
    strain    = cfg["strain"]
    threads   = cfg.get("threads", 56)

    print(f"  base_name: {base_name}")
    print(f"  strain:    {strain}")
    print(f"  threads:   {threads}")

    # Determine reference FASTA path
    default_ref = PIPELINE_DIR / f"results/{base_name}/_pipeline/assembly_{strain}/assembly_{strain}_dorado_reference.fasta"
    custom_ref  = cfg.get("reference_fasta", "")

    if custom_ref:
        ref_fasta = custom_ref
        ref_dir   = str(Path(custom_ref).parent)
        print(f"  reference: {ref_fasta}  (custom)")
    elif default_ref.exists():
        ref_fasta = str(default_ref)
        ref_dir   = str(default_ref.parent)
        print(f"  reference: {ref_fasta}")
    else:
        print(f"  WARNING: default reference not found at:\n"
              f"             {default_ref}")
        print("  The script will use whatever REFERENCE_FASTA is set in label_regions.sh.")
        ref_fasta = None
        ref_dir   = None

    # Patch label_regions.sh
    subs = {
        "BASE_NAME":        base_name,
        "STRAIN_ID":        strain,
        "THREADS":          str(threads),
        "YPRIME_LINKAGE":   cfg.get("yprime_linkage", "average"),
        "YPRIME_STOP_MODE": cfg.get("yprime_stop_mode", "silhouette"),
    }
    patched = patch_script(PIPELINE_DIR / "_pipeline" / "label_regions.sh", subs)

    # Also fix REFERENCE_DIR and REFERENCE_FASTA (script has hardcoded test values)
    if ref_dir:
        patched = substitute_var(patched, "REFERENCE_DIR",  ref_dir,   "label_regions.sh")
    if ref_fasta:
        patched = substitute_var(patched, "REFERENCE_FASTA", ref_fasta, "label_regions.sh")

    run_script(patched, "label_regions.sh", dry_run=dry_run)


def step_recombination(cfg: dict, dry_run: bool = False):
    print("\n" + "=" * 70)
    print("STEP 3: Recombination Analysis  (snakemake recombination_summary + single-sample plots + events + track plots)")
    print("=" * 70)

    samples = cfg.get("timepoint_samples") or []
    threads = cfg.get("threads", 56)
    day0_base = cfg.get("day0_base_name") or cfg.get("base_name")

    # Fallback: if no explicit timepoint_samples, use base_name as the sample
    if not samples:
        if cfg.get("base_name"):
            samples = [cfg["base_name"]]
        else:
            print("ERROR: No samples for recombination step. Use --samples or set timepoint_samples in config.")
            sys.exit(1)

    if not day0_base:
        print("ERROR: day0_base_name (or base_name) must be set so recombination rules know which reference to use.")
        sys.exit(1)

    print(f"  day0 reference sample: {day0_base}")
    print(f"  timepoint samples:     {samples}")
    print(f"  threads:               {threads}")

    for sample in samples:
        print(f"\n  --- Sample: {sample} ---")

        # Rewrite config.yaml so the Snakefile targets this time point
        # but keeps day0_base_name pinned to the reference sample.
        per_sample_cfg = dict(cfg)
        per_sample_cfg["base_name"] = sample
        per_sample_cfg["day0_base_name"] = day0_base

        # `rule all` aggregates every recombination-step output (summary,
        # single-sample plots, events extraction, per-chr_end track plots),
        # so a single `snakemake all` triggers all the rules that should
        # produce them. Snakemake skips anything already cached.
        cmd = ["snakemake", "-s", "_pipeline/Snakefile", "all", "-c", str(threads)]

        if dry_run:
            print(f"  [DRY RUN] Would write config.yaml with base_name={sample}, day0_base_name={day0_base}")
            print(f"  [DRY RUN] Would run: {' '.join(cmd)}")
            continue

        write_snakemake_config(per_sample_cfg)
        print(f"  Executing: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, cwd=PIPELINE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"\nERROR: snakemake recombination_summary/through_y_prime_analysis failed for {sample} (exit code {e.returncode})")
            sys.exit(e.returncode)


# ──────────────────────────────────────────────────────────────────────────────
# Sample config generator
# ──────────────────────────────────────────────────────────────────────────────

SAMPLE_CONFIG = """\
# pipeline_config.yaml — TeloTracker single-config file
# Edit the values below, then run:
#   python run_pipeline.py --config pipeline_config.yaml

# ── Required ──────────────────────────────────────────────────────────────────

# Base name of your day-0 sample (no file extension, no .bam)
base_name: "dorado_7302_day0_PromethION_no_tag_yes_rejection"

# Yeast strain identifier
strain: "7302"

# Directory containing input BAM/FASTQ files
bam_dir: "samples_dorado_basecalled"

# ── Dorado (required for create_ref) ─────────────────────────────────────────

# "docker" uses a Singularity .sif image; "local" uses dorado from PATH
dorado_mode: "docker"

# Full path to the Dorado Singularity image (only needed when dorado_mode: docker)
dorado_image: "/path/to/dorado.sif"

# Basecaller model that was used to generate the reads
dorado_model: "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"

# ── Optional / Defaults ───────────────────────────────────────────────────────

threads: 56
anchor_set: "telomerase_shutoff_anchors"

# BLAST parameter — minimum raw gapped alignment score
min_raw_gapped_score: 5000

# Y prime clustering: percent identity threshold for grouping
# Higher = more groups (finer resolution, risk of misassignment with ONT error)
# Lower  = fewer groups (coarser, more reliable)
# Only consulted when yprime_stop_mode is "threshold".
yprime_identity_threshold: 97.0

# Y prime clustering: hierarchical linkage method
#   complete = merge clusters only when EVERY cross-cluster pair is above threshold (strict)
#   average  = merge when the MEAN cross-cluster pair is above threshold (recommended default)
#   single   = merge if ANY cross-cluster pair is above threshold (most permissive)
yprime_linkage: "average"

# Y prime clustering: how to pick the cluster count
#   silhouette = data-driven; picks k with the highest silhouette score (recommended default)
#   threshold  = fixed cutoff from yprime_identity_threshold above
yprime_stop_mode: "silhouette"

# Override the automatically-extracted Y prime library with a curated one.
# Uncomment and set if you want to use a specific FASTA.
# y_prime_lib_override: "_pipeline/references/7302_features/repeatmasker_7302_all_y_primes.fasta"

# Uncomment to supply a custom reference FASTA for the label_regions step
# (defaults to the output of create_ref.sh)
# reference_fasta: "path/to/custom_reference.fasta"

# ── Recombination Analysis — Step 3 ──────────────────────────────────────────

# Day-0 sample that has already had create_ref + label_regions run against it.
# The recombination step uses its reference when analyzing later time points.
# Defaults to base_name if omitted.
# day0_base_name: "dorado_7302_day0_PromethION_no_tag_yes_rejection"

# Time-point samples to analyze against the day-0 reference.
# Each is run as: snakemake recombination_summary -c THREADS with the
# appropriate base_name set in config.yaml.
timepoint_samples:
  - "dorado_7302_day3_PromethION_no_tag_yes_rejection"
  - "dorado_7302_day6_PromethION_no_tag_yes_rejection"
"""


def init_config(output_path: str = "pipeline_config.yaml"):
    dest = Path(output_path)
    if dest.exists():
        ans = input(f"{dest} already exists. Overwrite? [y/N] ").strip().lower()
        if ans != "y":
            print("Aborted.")
            return
    dest.write_text(SAMPLE_CONFIG)
    print(f"Sample config written to: {dest}")
    print("Edit it, then run:")
    print(f"  python run_pipeline.py --config {dest}")


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="TeloTracker pipeline wrapper — one config to rule them all.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Meta
    parser.add_argument("--init-config", metavar="FILE", nargs="?",
                        const="pipeline_config.yaml",
                        help="Write a sample config file and exit "
                             "(default filename: pipeline_config.yaml)")
    parser.add_argument("--config", metavar="FILE",
                        help="Path to pipeline_config.yaml")
    parser.add_argument("--steps", nargs="+",
                        choices=["create_ref", "label_regions", "recombination", "day0", "all"],
                        default=["day0"],
                        metavar="STEP",
                        help="Steps to run: create_ref, label_regions, recombination, day0, all "
                             "(default: day0 = create_ref + label_regions). 'all' adds recombination. "
                             "Multiple individual steps can also be listed.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be executed without running anything")
    parser.add_argument("--no-organize", action="store_true",
                        help="Skip the post-step pass that creates sibling symlinks "
                             "(reference.fasta, labels.bed, ...) at the top of "
                             "results/<base>/ pointing into _pipeline/")

    # Core parameters
    grp = parser.add_argument_group("core parameters (override config file)")
    grp.add_argument("--base-name", metavar="NAME",
                     help="Day-0 sample base name (no extension)")
    grp.add_argument("--strain", metavar="ID",
                     help="Yeast strain number (e.g. 7302)")
    grp.add_argument("--bam-dir", metavar="DIR",
                     help="Directory containing input BAM/FASTQ files "
                          "(default: samples_dorado_basecalled)")
    grp.add_argument("--threads", type=int, metavar="N",
                     help="CPU threads (default: 56)")
    grp.add_argument("--anchor-set", metavar="NAME",
                     help="Anchor FASTA set name (default: test_anchors)")

    # Dorado
    grp2 = parser.add_argument_group("dorado parameters (required for create_ref)")
    grp2.add_argument("--dorado-mode", choices=["docker", "local"],
                      help="Dorado execution mode: 'docker' (Singularity) or 'local'")
    grp2.add_argument("--dorado-image", metavar="PATH",
                      help="Path to Dorado Singularity .sif file (docker mode)")
    grp2.add_argument("--dorado-model", metavar="MODEL",
                      help="Basecaller model name (e.g. dna_r10.4.1_e8.2_400bps_sup@v5.2.0)")

    # Label regions
    grp3 = parser.add_argument_group("label_regions parameters")
    grp3.add_argument("--reference-fasta", metavar="PATH",
                      help="Custom reference FASTA for label_regions "
                           "(default: output of create_ref.sh)")

    # Recombination
    grp4 = parser.add_argument_group("recombination parameters")
    grp4.add_argument("--samples", nargs="+", metavar="NAME",
                      help="Sample base names to analyze (runs snakemake recombination_summary "
                           "for each). Overrides timepoint_samples in config.")

    return parser.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # --init-config: write sample config and exit
    if args.init_config is not None:
        init_config(args.init_config)
        return

    # Load config file if given
    file_cfg = {}
    if args.config:
        if not os.path.exists(args.config):
            print(f"ERROR: config file not found: {args.config}")
            sys.exit(1)
        file_cfg = load_yaml_config(args.config)
        print(f"Loaded config from: {args.config}")

    cfg = merge_configs(file_cfg, args)

    # Resolve steps — expand shortcuts, preserve order, deduplicate
    steps = []
    for s in args.steps:
        if s == "all":
            steps.extend(["create_ref", "label_regions", "recombination"])
        elif s == "day0":
            steps.extend(["create_ref", "label_regions"])
        else:
            steps.append(s)
    seen = set()
    steps = [s for s in steps if not (s in seen or seen.add(s))]

    validate_cfg(cfg, steps)

    print("\nTeloTracker Pipeline Wrapper")
    print(f"Steps: {' → '.join(steps)}")
    if args.dry_run:
        print("[DRY RUN — nothing will be executed]")

    # Ensure we run from the pipeline root
    os.chdir(PIPELINE_DIR)

    organize_enabled = not args.no_organize
    day0_base = cfg.get("base_name", "")
    strain = cfg.get("strain", "")
    # create_ref.sh strips a trailing .bam, so strip here too for symlink-target resolution
    if day0_base.endswith(".bam"):
        day0_base = day0_base[:-4]

    if "create_ref" in steps:
        step_create_ref(cfg, dry_run=args.dry_run)
        organize_outputs(day0_base, strain, enabled=organize_enabled, dry_run=args.dry_run)

    if "label_regions" in steps:
        step_label_regions(cfg, dry_run=args.dry_run)
        organize_outputs(day0_base, strain, enabled=organize_enabled, dry_run=args.dry_run)

    if "recombination" in steps:
        step_recombination(cfg, dry_run=args.dry_run)
        # Recombination runs once per time-point sample; organize each one.
        for sample in (cfg.get("timepoint_samples") or [day0_base]):
            organize_outputs(sample, strain, enabled=organize_enabled, dry_run=args.dry_run)

    print("\n" + "=" * 70)
    print("Wrapper finished.")
    print("=" * 70)


if __name__ == "__main__":
    main()
