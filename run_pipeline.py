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

  # Read profiling only (telomere length, Y' count per end; no strain or day-0 reference):
  python run_pipeline.py --base-name my_sample --steps profile

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

# ------------------------------------------------------------------------------
# Config helpers
# ------------------------------------------------------------------------------

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
        "input_type":     args.input_type,
    }
    for key, val in cli_map.items():
        if val is not None:
            cfg[key] = val

    if args.no_adapter_check:
        cfg["adapter_check"] = False

    # Samples for the analyze step
    if args.samples:
        cfg["timepoint_samples"] = args.samples

    return cfg


def normalize_cfg(cfg: dict) -> dict:
    """Coerce config values once, before anything uses them.

    - strain and sample names become strings (YAML reads `strain: 7302` as an int);
    - a trailing .bam is dropped from base_name and the time-point samples (the Snakefile
      finds the input file's extension itself, so names never carry it);
    - yprime_grouping ("silhouette", or a percent identity such as 99.5 for a fixed cutoff)
      is resolved into yprime_stop_mode / yprime_identity_threshold, which older configs
      set directly and which are still accepted.
    """
    def _name(v):
        v = str(v)
        return v[:-4] if v.endswith(".bam") else v

    if cfg.get("strain") is not None:
        cfg["strain"] = str(cfg["strain"])
    if cfg.get("base_name"):
        cfg["base_name"] = _name(cfg["base_name"])
    if cfg.get("timepoint_samples"):
        cfg["timepoint_samples"] = [_name(x) for x in cfg["timepoint_samples"]]

    grouping = cfg.get("yprime_grouping")
    if grouping is not None:
        if str(grouping).strip().lower() == "silhouette":
            cfg["yprime_stop_mode"] = "silhouette"
        else:
            try:
                cfg["yprime_identity_threshold"] = float(grouping)
            except (TypeError, ValueError):
                print(f"ERROR: yprime_grouping must be 'silhouette' or a percent identity "
                      f"such as 99.5, not '{grouping}'")
                sys.exit(1)
            cfg["yprime_stop_mode"] = "threshold"
    return cfg


def pipeline_path(path) -> Path:
    """A config path: absolute, or relative to the pipeline directory."""
    path = Path(path)
    return path if path.is_absolute() else PIPELINE_DIR / path


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
        if mode not in ("docker", "local"):
            errors.append(f"dorado_mode must be 'docker' or 'local', not '{mode}'")
        elif mode == "docker":
            image = cfg.get("dorado_image")
            if not image:
                errors.append("dorado_mode is 'docker' but dorado_image is not set: point it at a "
                              "Singularity .sif with dorado")
            elif not pipeline_path(image).is_file():
                errors.append(f"dorado_image not found: {image}")

    if "label_regions" in steps:
        for key in ("label_yprime_fasta", "label_x_fasta"):
            if cfg.get(key) and not pipeline_path(cfg[key]).is_file():
                errors.append(f"{key} not found: {cfg[key]}")

    for step in ("profile", "recombination"):
        if step in steps and not step_samples(cfg):
            errors.append(f"At least one sample is required for the {step} step. "
                          "Use --samples, timepoint_samples or base_name in the config file.")
    if "recombination" in steps:
        if not cfg.get("strain"):
            errors.append("--strain (or strain in config) is required for recombination: "
                          "it names the day-0 reference, label and library files")

    if "label_regions" in steps:
        stop_mode = cfg.get("yprime_stop_mode", "silhouette")
        if stop_mode not in ("silhouette", "threshold"):
            errors.append(f"yprime_stop_mode must be 'silhouette' or 'threshold', not '{stop_mode}'")
        threshold = cfg.get("yprime_identity_threshold", 97.0)
        try:
            threshold = float(threshold)
            if not 50 <= threshold <= 100:
                errors.append(f"the Y' identity threshold is a percent between 50 and 100 "
                              f"(e.g. 99.5, not 0.995); got {threshold}")
        except (TypeError, ValueError):
            errors.append(f"the Y' identity threshold must be a number, not '{threshold}'")
        if (stop_mode == "silhouette" and "yprime_identity_threshold" in cfg
                and "yprime_grouping" not in cfg):
            print("WARNING: yprime_identity_threshold is set but Y' groups are chosen by silhouette, "
                  "which uses the threshold only when there are two Y' variants. For a fixed "
                  "cutoff set yprime_grouping: <percent> (e.g. 99.5).")

    try:
        resolve_yprime_reference(cfg, recombination=False)
    except ValueError as e:
        errors.append(str(e))

    # anchor_set now selects the anchor FASTA itself (not just output filenames), so a
    # missing file must fail here with a clear message rather than deep inside BLAST.
    anchor_set = cfg.get("anchor_set", "telomerase_shutoff_anchors")
    anchors_fa = PIPELINE_DIR / "_pipeline" / "references" / f"{anchor_set}.fasta"
    if not anchors_fa.is_file():
        errors.append(f"anchor_set '{anchor_set}' has no anchor file at {anchors_fa} "
                      f"(expected _pipeline/references/<anchor_set>.fasta)")
    input_type = cfg.get("input_type", "reads")
    if input_type not in ("reads", "contigs"):
        errors.append(f"input_type must be 'reads' or 'contigs', not '{input_type}'")
    elif input_type == "contigs":
        for step in ("create_ref", "recombination"):
            if step in steps:
                errors.append(f"{step} needs reads: input_type 'contigs' is for --steps profile only")
    if not isinstance(cfg.get("adapter_check", True), bool):
        errors.append(f"adapter_check must be true or false, not '{cfg.get('adapter_check')}'")
    for key, default in (("min_read_length", 2000), ("min_qscore", 10)):
        try:
            if float(cfg.get(key, default)) < 0:
                raise ValueError
        except (TypeError, ValueError):
            errors.append(f"{key} must be a number >= 0, not '{cfg.get(key)}'")

    if "create_ref" in steps and cfg.get("polish_scope", "genome") not in ("genome", "ends"):
        errors.append(f"polish_scope must be 'genome' or 'ends', not '{cfg.get('polish_scope')}'")
    if "create_ref" in steps and cfg.get("base_reference"):
        base_ref = PIPELINE_DIR / cfg["base_reference"]
        if not base_ref.is_file():
            errors.append(f"base_reference not found: {base_ref}")

    if errors:
        for e in errors:
            print(f"ERROR: {e}")
        sys.exit(1)


# ------------------------------------------------------------------------------
# Shell script patching
# ------------------------------------------------------------------------------

def substitute_var(content: str, var_name: str, new_value: str, script_name: str = "") -> str:
    """
    Replace a bash variable assignment line like:
        VAR_NAME="old"             -> VAR_NAME="new_value"
        VAR_NAME=old               -> VAR_NAME="new_value"
        VAR_NAME='old'             -> VAR_NAME="new_value"
        export VAR_NAME="old"      -> export VAR_NAME="new_value"
        VAR_NAME="old"  # comment  -> VAR_NAME="new_value"  # comment

    Handles quoted and unquoted values, leading indentation, `export`, trailing comments and
    CRLF line endings. Anchored to line-start so commented-out lines (#VAR_NAME=...) are ignored.
    The new value is written inside double quotes with double quotes and backslashes escaped;
    `$` still expands, as it always has. A variable that is not found is an error: the script
    would otherwise run with its own built-in default.
    """
    # Group 1: indentation, optional "export ", "VAR="
    # Group 2: the current value -- a double-quoted string, a single-quoted string or a bare word
    # Group 3: optional trailing whitespace + comment
    # Group 4: optional carriage return (CRLF files)
    pattern = (rf'^([ \t]*(?:export[ \t]+)?{re.escape(var_name)}=)'
               r'("(?:[^"\\]|\\.)*"|\'[^\']*\'|[^\s#]*)([ \t]*#[^\r\n]*)?(\r?)$')
    quoted = '"' + new_value.replace('\\', '\\\\').replace('"', '\\"') + '"'
    new_content, n = re.subn(pattern, lambda m: m.group(1) + quoted + (m.group(3) or '') + m.group(4),
                             content, flags=re.MULTILINE)
    if n == 0:
        sys.exit(f"ERROR: could not find '{var_name}=' in {script_name} to substitute; "
                 f"the script would run with its built-in default.")
    return new_content


def patch_script(script_path: Path, subs: dict) -> str:
    """Apply multiple variable substitutions to a script file."""
    content = script_path.read_text()
    for var, value in subs.items():
        content = substitute_var(content, var, str(value), script_path.name)
    return content


# ------------------------------------------------------------------------------
# Snakemake config.yaml writer
# ------------------------------------------------------------------------------

def write_snakemake_config(cfg: dict, recombination: bool = True, dest: Path = None, dry_run: bool = False):
    """Write the config.yaml consumed by the Snakemake pipeline (print it instead on a dry run)."""
    if dest is None:
        dest = PIPELINE_DIR / "_pipeline" / "config.yaml"
    text = snakemake_config_text(cfg, recombination)
    if dry_run:
        print(f"  [DRY RUN] Would write {dest}:")
        for line in text.splitlines():
            print(f"    {line}")
        return
    dest.write_text(text)
    print(f"  Written: {dest}")


def day0_files(cfg: dict) -> dict:
    """Where recombination finds the day-0 sample's files, as written into config.yaml. The three
    libraries keep a literal {strain}, which the Snakefile fills in."""
    strain = cfg["strain"]
    # Day-0 sample whose reference + labels are used as the read-against.
    # Internal-only key: run_pipeline.py sets this during the recombination
    # step so the per-timepoint config rewrite keeps `base_name` pointing at
    # the timepoint sample while the library paths still target day-0.
    day0_base = cfg.get("_day0_base", cfg["base_name"])
    labels = f"results/{day0_base}/_pipeline/pretelomeric_labels"
    return {
        # A supplied reference (reference_fasta, e.g. a curated strain genome labelled with
        # label_regions) is the one recombination must align to; otherwise create_ref's output.
        "day0_ref": cfg.get("reference_fasta")
                    or f"results/{day0_base}/_pipeline/assembly_{strain}/assembly_{strain}_dorado_reference.fasta",
        "day0_bed": f"{labels}/pretelomeric_regions_{strain}_simp.bed",
        "y_prime_lib": f"{labels}/extracted_yprimes_{{strain}}.fasta",
        "spacer_lib_dir": f"{labels}/pairings_for_spacers/{{strain}}_pairings",
        "x_element_lib": f"{labels}/clustered_x_elements_{{strain}}.fasta",
    }


def snakemake_config_text(cfg: dict, recombination: bool = True) -> str:
    """The config.yaml text: the key structure the current Snakefile reads at module load time.

    With recombination off (profile, create_ref) the Snakefile loads only the profiling rules,
    so no strain or day-0 settings are needed."""
    base_name    = cfg["base_name"]
    strain       = cfg.get("strain")
    bam_dir      = cfg.get("bam_dir", "samples_dorado_basecalled")
    anchor_set   = cfg.get("anchor_set", "telomerase_shutoff_anchors")
    min_gapped   = cfg.get("min_raw_gapped_score", 5000)
    yp_override  = cfg.get("y_prime_lib_override")

    lines = [
        "# config.yaml -- auto-generated by run_pipeline.py",
        f'base_name: "{base_name}"',
        f'anchor_set: "{anchor_set}"',
        f'bam_dir: "{bam_dir}"',
    ]
    if strain:
        lines.append(f'strain: "{strain}"')
    lines += [
        f"recombination: {'true' if recombination else 'false'}",
        "",
        "# BLAST parameter",
        f"min_raw_gapped_score: {min_gapped}",
        "",
        "# Profiling: input_type reads (basecalls) or contigs (assemblies); read filters; whether reads",
        "# must show the sequencing adapter after the telomere (off for adapter-trimmed or public reads)",
        f'input_type: "{cfg.get("input_type", "reads")}"',
        f"adapter_check: {'true' if adapter_check_on(cfg) else 'false'}",
        f"min_read_length: {cfg.get('min_read_length', 2000)}",
        f"min_qscore: {cfg.get('min_qscore', 10)}",
        "",
    ]
    if recombination:
        lines += [
            "# Recombination calling (the Snakefile's defaults when absent; written so the",
            "# recombination rules and the onion-skin time course use the same values)",
            f'attribution_mode: "{cfg.get("attribution_mode", "v2")}"',
            f'y_prime_id_level: "{cfg.get("y_prime_id_level", "family")}"',
            "",
        ]
    yp_path, yp_label = resolve_yprime_reference(cfg, recombination)
    lines += [
        "# Expected Y' copies per end: the zero line of the delta-Y' plots (yprime_reference;",
        "# an empty path means none, and the plots show absolute counts)",
        f'yprime_reference_path: "{yp_path}"',
        f'yprime_reference_label: "{yp_label}"',
        "",
        "# Paths to references",
        "references:",
        # Derived from anchor_set so a strain-specific set is actually used, not just
        # baked into output filenames. With the default this line is unchanged.
        f'  anchors: "_pipeline/references/{anchor_set}.fasta"',
        '  adapters: "_pipeline/references/nanopore_sqk-slk114_adapter_sequence_truncated.txt"',
        '  probe: "_pipeline/references/y_prime_probe.fasta"',
        "",
    ]
    if recombination:
        day0 = day0_files(cfg)
        lines += [
            f'  day0_ref: "{day0["day0_ref"]}"',
            f'  day0_bed: "{day0["day0_bed"]}"',
            f'  y_prime_lib: "{day0["y_prime_lib"]}"',
        ]
        if yp_override:
            lines.append(f'  y_prime_lib_override: "{yp_override}"')
        lines += [
            f'  spacer_lib_dir: "{day0["spacer_lib_dir"]}"',
            f'  x_element_lib: "{day0["x_element_lib"]}"',
            "",
        ]
    return "\n".join(lines)


def resolve_yprime_reference(cfg: dict, recombination: bool):
    """(file, label) of the expected Y' counts for cfg's yprime_reference (default auto); file is
    relative to the pipeline directory when inside it, '' for none. In recombination runs, auto
    falls back to the day-0 sample's label BED. Raises ValueError for an unknown name or file."""
    yr = _scripts_import("yprime_reference")
    value = cfg.get("yprime_reference", "auto")
    text = "" if value is None else str(value).strip()
    if text.lower() not in ("", "auto", "none") and ("/" in text or text.endswith((".bed", ".tsv"))):
        value = str(PIPELINE_DIR / text)          # files are given relative to the pipeline directory
    day0_bed = (str(PIPELINE_DIR / day0_files(cfg)["day0_bed"])
                if recombination and cfg.get("strain") and cfg.get("base_name") else None)
    path, label = yr.resolve(value, cfg.get("strain"), day0_bed, str(PIPELINE_DIR / "_pipeline" / "references"))
    if path:
        rel = os.path.relpath(path, PIPELINE_DIR)
        path = path if rel.startswith("..") else rel
    return path or "", label


def adapter_check_on(cfg: dict) -> bool:
    """Whether profiling checks for the sequencing adapter after the telomere (never for contigs)."""
    return bool(cfg.get("adapter_check", True)) and cfg.get("input_type", "reads") == "reads"


def missing_sample_inputs(cfg: dict, samples: list, root: Path = None) -> list:
    """One line per sample whose input file is not in bam_dir (by input_type's extensions)."""
    pre = _scripts_import("prerequisites")
    candidates = pre.INPUT_TYPES.get(cfg.get("input_type", "reads"), pre.READ_INPUTS)
    input_dir = str(Path(root or PIPELINE_DIR) / cfg.get("bam_dir", "samples_dorado_basecalled"))
    return [f"input for sample {s}: none of {pre.tried_inputs(input_dir, s, candidates)}"
            for s in samples if pre.find_sample_input(input_dir, s, candidates) is None]


def step_samples(cfg: dict) -> list:
    """The samples profile / recombination run on: timepoint_samples, else the day-0 sample."""
    if cfg.get("timepoint_samples"):
        return list(cfg["timepoint_samples"])
    return [cfg["base_name"]] if cfg.get("base_name") else []


def _scripts_import(module: str):
    """Import a helper module from _pipeline/scripts (shared with the Snakefile)."""
    scripts = str(PIPELINE_DIR / "_pipeline" / "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    return __import__(module)


def check_recombination_prereqs(cfg: dict, samples: list, root: Path = None) -> list:
    """Problems (one line each) that would stop recombination: the day-0 sample's reference,
    labels and libraries, and each sample's input file. Empty when everything is there.
    Paths are relative to root (the pipeline directory)."""
    pre = _scripts_import("prerequisites")
    root = Path(root or PIPELINE_DIR)
    day0 = {k: v.replace("{strain}", cfg["strain"]) for k, v in day0_files(cfg).items()}
    if cfg.get("y_prime_lib_override"):
        day0["y_prime_lib"] = cfg["y_prime_lib_override"]
    problems = [f"{what}: {path}  ({hint})" for what, path, hint in pre.missing_day0_inputs(
        *(str(root / day0[k]) for k in ("day0_ref", "day0_bed", "y_prime_lib",
                                         "spacer_lib_dir", "x_element_lib")))]
    return problems + missing_sample_inputs(cfg, samples, root)


# ------------------------------------------------------------------------------
# Script runner
# ------------------------------------------------------------------------------

# Prepended to every patched shell script before execution so that
# `conda activate` works even though subprocess bash is non-interactive
# (shell functions aren't inherited from the parent -- env vars are,
# so `CONDA_SHLVL` alone is NOT a reliable "conda is ready" signal).
# Check for the `conda` shell function itself. Runs above any `set -u`
# / `set -e` so a missing conda.sh falls through silently rather than
# aborting before we see the real error.
CONDA_INIT_PROLOGUE = (
    "# Auto-injected by run_pipeline.py -- initialize conda for non-interactive shells\n"
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
    os.unlink(tmp_path)


# ------------------------------------------------------------------------------
# Post-step organization
# ------------------------------------------------------------------------------

def organize_outputs(base_name: str, strain: str, enabled: bool = True, dry_run: bool = False):
    """Run _pipeline/scripts/organize_outputs.py to (re)create sibling symlinks + MANIFEST."""
    if not enabled or dry_run:
        return
    script = PIPELINE_DIR / "_pipeline" / "scripts" / "organize_outputs.py"
    # without a strain (profile runs) the script detects it, or links only strain-free outputs
    cmd = [sys.executable, str(script), base_name, "--quiet"] + (["--strain", strain] if strain else [])
    try:
        subprocess.run(cmd, check=True, cwd=PIPELINE_DIR)
    except subprocess.CalledProcessError as e:
        print(f"  WARNING: organize_outputs exited with code {e.returncode} -- symlinks may be incomplete")


# ------------------------------------------------------------------------------
# Pipeline steps
# ------------------------------------------------------------------------------

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

    print(f"  base_name:    {base_name}")
    print(f"  strain:       {strain}")
    print(f"  threads:      {threads}")
    print(f"  dorado_mode:  {dorado_mode}")
    if dorado_image:
        print(f"  dorado_image: {dorado_image}")
    print(f"  dorado_model: {dorado_model}")
    print(f"  polish_scope: {cfg.get('polish_scope', 'genome')}")

    # 1. Write Snakemake config.yaml (create_ref.sh runs only the profiling rules)
    print("\n  Writing config.yaml for Snakemake...")
    write_snakemake_config(cfg, recombination=False, dry_run=dry_run)

    # 2. Patch create_ref.sh
    # ANCHOR_SET must match the Snakemake config: create_ref.sh builds the
    # all_matches BLAST path from it for consensus scaffold selection, so a
    # mismatch here surfaces as a missing-file error partway through the run.
    subs = {
        "BASE_NAME":    base_name,
        "STRAIN_ID":    strain,
        "THREADS":      str(threads),
        "DORADO_MODE":  dorado_mode,
        "DORADO_MODEL": dorado_model,
        "ANCHOR_SET":   cfg.get("anchor_set", "telomerase_shutoff_anchors"),
        # The base reference must be truncated at the same anchors being used. Defaults
        # to the literal 6991 path so existing configs reproduce exactly; a new strain
        # sets base_reference to its own <strain>_only_to_anchors.fasta.
        "REFERENCE":    cfg.get("base_reference",
                                "_pipeline/references/6991_only_to_anchors.fasta"),
        "POLISH_SCOPE": cfg.get("polish_scope", "genome"),
        "POLISH_END_FLANK": str(cfg.get("polish_end_flank", 5000)),
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
        # Previously written only into the Snakemake config and never passed to the
        # clustering script, so threshold mode always cut at 97 whatever was set.
        "YPRIME_IDENTITY_THRESHOLD": cfg.get("yprime_identity_threshold", 97.0),
        "YPRIME_DEDUP_THRESHOLD":    cfg.get("yprime_dedup_threshold", 99.9),
        "YPRIME_BOUNDARY_TRIM":      "true" if cfg.get("yprime_boundary_trim", True) else "false",
    }
    # label_regions.sh hardcodes ANCHORS_FASTA=test_anchors.fasta, which is byte-identical
    # in content to telomerase_shutoff_anchors.fasta -- which is why nobody noticed it was
    # never patched. Only override for a non-default set, so existing runs reproduce the
    # patched script exactly.
    anchor_set = cfg.get("anchor_set", "telomerase_shutoff_anchors")
    if anchor_set != "telomerase_shutoff_anchors":
        subs["ANCHORS_FASTA"] = f"_pipeline/references/{anchor_set}.fasta"
    # The Y' and X sequences labelling searches the reference with (default: 6991-derived libraries).
    for key, var in (("label_yprime_fasta", "YPRIMES_FASTA"), ("label_x_fasta", "XPRIMES_FASTA")):
        if cfg.get(key):
            subs[var] = cfg[key]

    patched = patch_script(PIPELINE_DIR / "_pipeline" / "label_regions.sh", subs)

    # Also fix REFERENCE_DIR and REFERENCE_FASTA (script has hardcoded test values)
    if ref_dir:
        patched = substitute_var(patched, "REFERENCE_DIR",  ref_dir,   "label_regions.sh")
    if ref_fasta:
        patched = substitute_var(patched, "REFERENCE_FASTA", ref_fasta, "label_regions.sh")

    run_script(patched, "label_regions.sh", dry_run=dry_run)



def write_sample_report(sample: str, cfg: dict):
    """Render the HTML summary for one sample.

    Deliberately best-effort and never fatal: it runs after a failed snakemake too, because a
    partial report showing which stages produced output is more useful than none. The report
    itself audits its inputs and marks anything missing, so an incomplete run cannot be
    mistaken for a complete one.
    """
    if not cfg.get("make_report", True):
        return
    script = os.path.join(PIPELINE_DIR, "_pipeline", "scripts", "make_sample_report.py")
    pipeline_dir = os.path.join("results", sample, "_pipeline")
    if not os.path.isfile(script):
        print(f"  NOTE: {script} not found; skipping HTML report")
        return
    if not os.path.isdir(os.path.join(PIPELINE_DIR, pipeline_dir)):
        print(f"  NOTE: no {pipeline_dir}; skipping HTML report")
        return
    cmd = [sys.executable, script, "--pipeline-dir", pipeline_dir, "--base-name", sample]
    # How the day-0 Y' groups were made (label_regions output), shown in the report header.
    if cfg.get("base_name") and cfg.get("strain"):
        grouping = os.path.join("results", cfg["base_name"], "_pipeline", "pretelomeric_labels",
                                f"yprime_grouping_{cfg['strain']}.tsv")
        if os.path.exists(os.path.join(PIPELINE_DIR, grouping)):
            cmd += ["--yprime-grouping", grouping]
    try:
        subprocess.run(cmd, check=True, cwd=PIPELINE_DIR)
    except Exception as exc:                                        # noqa: BLE001
        # must never turn a good run into a failed one, nor mask a real error
        print(f"  WARNING: HTML report generation failed for {sample}: {exc}")


def onion_skin_timecourse(cfg: dict, samples: list, day0_base: str, dry_run: bool = False):
    """Follow Y' onion-skin layers across the time points (onion_skin_timecourse.py).

    Runs once, after every sample's recombination step, because it compares samples with each
    other. timepoint_samples must be listed in time order. Best-effort like the HTML report:
    it reads finished outputs and never turns a good run into a failed one."""
    if not cfg.get("onion_skin_timecourse", True) or len(samples) < 2:
        return
    strain = cfg["strain"]
    labels = f"results/{day0_base}/_pipeline/pretelomeric_labels"
    lib = cfg.get("y_prime_lib_override") or f"{labels}/extracted_yprimes_{strain}.fasta"
    out_dir = f"results/{day0_base}/_pipeline/onion_skin_timecourse"
    cmd = [sys.executable, os.path.join("_pipeline", "scripts", "onion_skin_timecourse.py"),
           "--results", "results", "--samples", *samples,
           "--day0-bed", f"{labels}/pretelomeric_regions_{strain}_simp.bed",
           "--y-prime-lib", lib,
           "--y-prime-id-level", cfg.get("y_prime_id_level", "family"),
           "--out-dir", out_dir, "--prefix", day0_base]
    print(f"\n  --- Onion skin across time points ({len(samples)} samples, in the order listed) ---")
    if dry_run:
        print(f"  [DRY RUN] Would run: {' '.join(cmd)}")
        return
    try:
        subprocess.run(cmd, check=True, cwd=PIPELINE_DIR)
    except Exception as exc:                                        # noqa: BLE001
        print(f"  WARNING: onion-skin timecourse failed: {exc}")


def step_profile(cfg: dict, dry_run: bool = False):
    """Read profiling only (anchors, chromosome ends, telomere length, Y' count, plots) for each
    sample: no strain and no day-0 reference needed."""
    print("\n" + "=" * 70)
    print("PROFILE: Read profiling  (snakemake profile: anchors, ends, telomere length, Y' count, plots)")
    print("=" * 70)

    samples = step_samples(cfg)
    threads = cfg.get("threads", 56)
    print(f"  samples: {samples}")
    print(f"  threads: {threads}")
    print(f"  input:   {cfg.get('input_type', 'reads')}, adapter check "
          f"{'on' if adapter_check_on(cfg) else 'off'}")

    missing = missing_sample_inputs(cfg, samples)
    if missing:
        print("\n  " + ("WARNING (dry run):" if dry_run else "ERROR: profiling cannot start; missing:"))
        for m in missing:
            print(f"    {m}")
        if not dry_run:
            sys.exit(1)

    for sample in samples:
        print(f"\n  --- Sample: {sample} ---")
        per_sample_cfg = dict(cfg)
        per_sample_cfg["base_name"] = sample
        cmd = ["snakemake", "-s", "_pipeline/Snakefile", "profile", "-c", str(threads)]
        if dry_run:
            write_snakemake_config(per_sample_cfg, recombination=False, dry_run=True)
            print(f"  [DRY RUN] Would run: {' '.join(cmd)}")
            continue
        write_snakemake_config(per_sample_cfg, recombination=False)
        print(f"  Executing: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, cwd=PIPELINE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"\nERROR: snakemake profile failed for {sample} (exit code {e.returncode})")
            write_sample_report(sample, cfg)
            sys.exit(e.returncode)
        write_sample_report(sample, cfg)


def step_recombination(cfg: dict, dry_run: bool = False):
    print("\n" + "=" * 70)
    print("STEP 3: Recombination Analysis  (snakemake all: summary, single-sample plots, events, track plots)")
    print("=" * 70)

    # timepoint_samples, or the day-0 sample itself (analysed against its own reference)
    samples = step_samples(cfg)
    threads = cfg.get("threads", 56)
    day0_base = cfg.get("base_name")

    if not samples:
        print("ERROR: No samples for recombination step. Use --samples or set timepoint_samples in config.")
        sys.exit(1)

    if not day0_base:
        print("ERROR: base_name must be set so recombination rules know which reference to use.")
        sys.exit(1)

    print(f"  day0 reference sample: {day0_base}")
    print(f"  timepoint samples:     {samples}")
    print(f"  threads:               {threads}")

    # Checked here rather than in validate_cfg so that `--steps label_regions recombination`
    # can first build the libraries this needs.
    problems = check_recombination_prereqs(cfg, samples)
    if problems:
        print("\n  " + ("WARNING (dry run; an earlier step would have to create these):" if dry_run
                        else "ERROR: recombination cannot start; these are missing:"))
        for p in problems:
            print(f"    {p}")
        if not dry_run:
            print("  Run label_regions on the day-0 sample (create_ref first, or set reference_fasta), "
                  "or check base_name, strain and bam_dir.")
            sys.exit(1)

    for sample in samples:
        print(f"\n  --- Sample: {sample} ---")

        # Rewrite config.yaml so the Snakefile targets this time point.
        # _day0_base is a private key consumed only by write_snakemake_config:
        # it makes the library/reference paths point back at the day-0 dir
        # while base_name points at the timepoint.
        per_sample_cfg = dict(cfg)
        per_sample_cfg["base_name"] = sample
        per_sample_cfg["_day0_base"] = day0_base

        # `rule all` aggregates every recombination-step output (summary,
        # single-sample plots, events extraction, per-chr_end track plots),
        # so a single `snakemake all` triggers all the rules that should
        # produce them. Snakemake skips anything already cached.
        cmd = ["snakemake", "-s", "_pipeline/Snakefile", "all", "-c", str(threads)]

        if dry_run:
            write_snakemake_config(per_sample_cfg, recombination=True, dry_run=True)
            print(f"  [DRY RUN] Would run: {' '.join(cmd)}")
            continue

        write_snakemake_config(per_sample_cfg, recombination=True)
        print(f"  Executing: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, cwd=PIPELINE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"\nERROR: snakemake all failed for {sample} (exit code {e.returncode})")
            # still write the report: it shows which stages did produce output, which is
            # exactly what you want when diagnosing the failure
            write_sample_report(sample, cfg)
            sys.exit(e.returncode)
        write_sample_report(sample, cfg)

    onion_skin_timecourse(cfg, samples, day0_base, dry_run=dry_run)


# ------------------------------------------------------------------------------
# Sample config generator
# ------------------------------------------------------------------------------

SAMPLE_CONFIG = """\
# pipeline_config.yaml -- TeloTracker single-config file
# Edit the values below, then run:
#   python run_pipeline.py --config pipeline_config.yaml

# -- Required ------------------------------------------------------------------

# Base name of your day-0 sample (no file extension, no .bam)
base_name: "dorado_7302_day0_PromethION_no_tag_yes_rejection"

# Yeast strain identifier (names the day-0 reference, label and library files; not needed for
# --steps profile). The curated day-0 strains ship as 6991_curated, 7172_curated, 7302_curated.
strain: "7302"

# Directory containing input BAM/FASTQ files
bam_dir: "samples_dorado_basecalled"

# -- Dorado (required for create_ref) -----------------------------------------

# "docker" uses a Singularity .sif image; "local" uses dorado from PATH
dorado_mode: "docker"

# Full path to the Dorado Singularity image (required when dorado_mode: docker; create_ref
# stops up front if it is missing)
dorado_image: "/path/to/dorado.sif"

# Basecaller model that was used to generate the reads
dorado_model: "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"

# -- Profiling (all steps; --steps profile runs only this part) -----------------

# reads   = BAM / FASTQ basecalls in bam_dir (default)
# contigs = an assembled genome, <base_name>.fasta / .fa / .fna (.gz) in bam_dir; profile only:
#           each contig is cut into two halves, one per chromosome end
input_type: "reads"

# Telomere reads count as "qualifying" (telomere plots, Y' plots, read summary) only when the
# sequencing adapter follows the telomere, which shows the read reached the molecule's end.
# Set false for reads without adapters (dorado's default trimming, SRA/ENA downloads): porechop
# is skipped and the telomere repeat alone decides. read_summary.tsv warns when reads look trimmed.
adapter_check: true

# Read filters (reads only): minimum length in bp, minimum mean qscore (qs:f tag; 0 = no filter)
min_read_length: 2000
min_qscore: 10

# Expected Y' copies per chromosome end, the zero line of the delta-Y' plots:
#   auto = the shipped reference named by strain (6991_curated, 7172_curated, 7302_curated, or a
#          table in _pipeline/references/yprime_expected_counts/ such as 7372), else the day-0
#          labels in recombination runs, else none
#   none = no reference: plots show absolute Y' copies per end
#   a shipped name (above, also 7021 7093 7154 7174 7250 7321 7323 7324 KRLT3_7637 tlc1),
#   or a path to a label BED / counts TSV (chr_end<TAB>expected_y_primes)
yprime_reference: "auto"

# -- Optional / Defaults -------------------------------------------------------

threads: 56

# Write a self-contained HTML summary (<sample>_report.html) into each sample's
# _pipeline/ directory at the end of the recombination step. Also written when a run
# fails part-way; it audits its own inputs so an incomplete run is obvious. false to skip.
make_report: true
# Anchor set: names _pipeline/references/<anchor_set>.fasta, the 32 anchors reads are
# BLASTed against. The default is cut from the 6991 genome; a different strain background
# needs its own set, and a matching truncated base reference for create_ref (below).
anchor_set: "telomerase_shutoff_anchors"
# base_reference: "_pipeline/references/6991_only_to_anchors.fasta"   # create_ref's base; default shown
# Polishing scope for create_ref: "genome" polishes the whole reference with every read;
# "ends" polishes each arm (extension + polish_end_flank bp of the base) with only the reads
# anchored to it and leaves the interiors as the base reference. Use "ends" when the library
# is end-enriched (interiors at a few x): genome-wide polishing at that depth corrupts them.
# polish_scope: "genome"
# polish_end_flank: 5000

# BLAST parameter -- minimum raw gapped alignment score
min_raw_gapped_score: 5000

# Y prime grouping: how Y' copies are grouped into IDs (ID1, ID2, ...).
#   silhouette = data-driven: picks the number of groups with the best silhouette score (default)
#   <percent>  = fixed identity cutoff, e.g. 97, 99 or 99.5. Higher = more, finer groups.
#                97 and 99 are the benchmarked cuts (~10 and ~12 groups on day-0 references at
#                the same false-mismatch rate); 99.5 is the cutoff of the published cut99.5
#                analyses (see _pipeline/references/day0/README.md).
# Group IDs depend on this setting, so compare samples grouped the same way.
# (Configs that set yprime_stop_mode / yprime_identity_threshold directly still work.)
yprime_grouping: "silhouette"

# Y prime clustering: variants at least this identical are merged before clustering
yprime_dedup_threshold: 99.9

# Correct Y' elements whose anchor-proximal end overhangs near-identical partner elements
# (chr16L_Y_Prime_1 is 77 bp too long in every 6991/7172/7302 reference, which splits it into
# its own group at a 99 % cut), and flag copies much shorter than their partners.
yprime_boundary_trim: true

# Y prime clustering: hierarchical linkage method
#   complete = merge clusters only when EVERY cross-cluster pair is above threshold (strict)
#   average  = merge when the MEAN cross-cluster pair is above threshold (recommended default)
#   single   = merge if ANY cross-cluster pair is above threshold (most permissive)
yprime_linkage: "average"

# Y' and X sequences label_regions searches a new reference with (default: 6991-derived).
# label_yprime_fasta: "_pipeline/references/repeatmasker_6991_all_y_primes.fasta"
# label_x_fasta: "_pipeline/references/6991_xprimes.fasta"

# Override the automatically-extracted Y prime library with a curated one.
# Uncomment and set if you want to use a specific FASTA.
# y_prime_lib_override: "_pipeline/references/7302_features/repeatmasker_7302_all_y_primes.fasta"

# Uncomment to supply your own reference FASTA (e.g. a curated strain genome) instead of
# building one with create_ref.sh: label_regions labels it, and the recombination step then
# aligns reads to it. Run --steps label_regions recombination (skip create_ref).
# reference_fasta: "path/to/custom_reference.fasta"

# -- Recombination Analysis -- Step 3 ------------------------------------------

# base_name above is treated as the day-0 reference sample for recombination:
# the step reads results/<base_name>/_pipeline/{assembly_<strain>,pretelomeric_labels}
# and analyzes each entry of timepoint_samples against it.

# Time-point samples to analyze against the day-0 reference, IN TIME ORDER.
# Each is run as: snakemake -s _pipeline/Snakefile all -c THREADS, with base_name in
# _pipeline/config.yaml set to that sample.
timepoint_samples:
  - "dorado_7302_day3_PromethION_no_tag_yes_rejection"
  - "dorado_7302_day6_PromethION_no_tag_yes_rejection"

# With two or more time points, compare them after the last one: arrays that gain further
# copies of the same Y' over time ("onion skin" layers). Reads timepoint_samples in the order
# listed. Output: results/<base_name>/_pipeline/onion_skin_timecourse/. false to skip.
onion_skin_timecourse: true
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


# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="TeloTracker pipeline wrapper -- one config to rule them all.",
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
                        choices=["profile", "create_ref", "label_regions", "recombination", "day0", "all"],
                        default=["day0"],
                        metavar="STEP",
                        help="Steps to run: profile, create_ref, label_regions, recombination, day0, all "
                             "(default: day0 = create_ref + label_regions). 'all' adds recombination. "
                             "'profile' runs read profiling only (anchors, ends, telomere length, Y' count) "
                             "and needs no strain or day-0 reference. "
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

    # Profiling
    grp5 = parser.add_argument_group("profiling parameters")
    grp5.add_argument("--input-type", choices=["reads", "contigs"],
                      help="reads (BAM/FASTQ basecalls, default) or contigs (assembly FASTA; profile only)")
    grp5.add_argument("--no-adapter-check", action="store_true",
                      help="reads lack the sequencing adapter (adapter-trimmed or public data): skip "
                           "porechop and count telomere reads by repeat length alone")

    # Recombination
    grp4 = parser.add_argument_group("recombination parameters")
    grp4.add_argument("--samples", nargs="+", metavar="NAME",
                      help="Sample base names to analyze (runs snakemake all for each). "
                           "Overrides timepoint_samples in config.")

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

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

    cfg = normalize_cfg(merge_configs(file_cfg, args))

    # Resolve steps -- expand shortcuts, preserve order, deduplicate
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
    print(f"Steps: {' -> '.join(steps)}")
    if args.dry_run:
        print("[DRY RUN -- nothing will be executed]")

    # Ensure we run from the pipeline root
    os.chdir(PIPELINE_DIR)

    organize_enabled = not args.no_organize
    day0_base = cfg.get("base_name", "")
    strain = cfg.get("strain", "")

    if "create_ref" in steps:
        step_create_ref(cfg, dry_run=args.dry_run)
        organize_outputs(day0_base, strain, enabled=organize_enabled, dry_run=args.dry_run)

    if "label_regions" in steps:
        step_label_regions(cfg, dry_run=args.dry_run)
        organize_outputs(day0_base, strain, enabled=organize_enabled, dry_run=args.dry_run)

    if "profile" in steps:
        step_profile(cfg, dry_run=args.dry_run)
        for sample in step_samples(cfg):
            organize_outputs(sample, strain, enabled=organize_enabled, dry_run=args.dry_run)

    if "recombination" in steps:
        step_recombination(cfg, dry_run=args.dry_run)
        # Recombination runs once per time-point sample; organize each one.
        for sample in step_samples(cfg):
            organize_outputs(sample, strain, enabled=organize_enabled, dry_run=args.dry_run)

    print("\n" + "=" * 70)
    print("Wrapper finished.")
    print("=" * 70)


if __name__ == "__main__":
    main()
