#!/usr/bin/env python3
"""
Organize teloTracker outputs.

Creates sibling symlinks at the base of results/<base_name>/ pointing to the
user-facing final outputs that live inside results/<base_name>/_pipeline/, and
writes a MANIFEST.md explaining the layout.

Idempotent: safe to call after every pipeline step. Missing final outputs are
skipped silently (e.g. recombination_summary.tsv isn't there yet while
create_ref.sh is still running). Stale symlinks are refreshed. Real files are
never overwritten.

Standalone:
    python scripts/organize_outputs.py <base_name> [--strain STRAIN] [--results-dir DIR]
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, Optional


FINALS: Dict[str, str] = {
    "reference.fasta":           "assembly_{strain}/assembly_{strain}_dorado_reference.fasta",
    "labels.bed":                "pretelomeric_labels/pretelomeric_regions_{strain}_simp.bed",
    "labels.gff3":               "pretelomeric_labels/pretelomeric_regions_{strain}.gff3",
    "quality_report.txt":        "pretelomeric_labels/pretelomeric_regions_{strain}_quality_report.txt",
    "recombination_summary.tsv": "recombination/{base}_recombination_summary.tsv",
    # Real subdirectory under _pipeline/, surfaced at the top of results/<base>/
    # via a directory symlink. ensure_symlink() handles both files and dirs.
    "recombination_events":      "recombination_events",
}


def detect_strain(pipeline_dir: Path) -> Optional[str]:
    if not pipeline_dir.is_dir():
        return None
    for sub in pipeline_dir.iterdir():
        m = re.match(r"^assembly_(.+)$", sub.name)
        if m and sub.is_dir():
            return m.group(1)
    return None


def ensure_symlink(link_path: Path, target_relative: str) -> str:
    """
    Create or refresh a symlink at link_path pointing to target_relative (relative
    to link_path's parent). Returns: "created" | "updated" | "skipped" |
    "blocked" (real file in the way) | "missing" (target doesn't exist).
    """
    target_full = link_path.parent / target_relative
    if not target_full.exists():
        return "missing"

    if link_path.is_symlink():
        if os.readlink(link_path) == target_relative:
            return "skipped"
        link_path.unlink()
        link_path.symlink_to(target_relative)
        return "updated"

    if link_path.exists():
        return "blocked"

    link_path.symlink_to(target_relative)
    return "created"


def write_manifest(base_dir: Path, present_finals) -> None:
    lines = [
        f"# {base_dir.name} — output layout",
        "",
        "## Final outputs (symlinks)",
        "",
        "The top-level entries below are pointers into `_pipeline/` for convenience.",
        "Deleting a symlink does not delete the underlying file.",
        "",
    ]
    if present_finals:
        lines.extend(f"- `{name}`" for name in present_finals)
    else:
        lines.append("_No final outputs produced yet — the pipeline may still be running._")
    lines += [
        "",
        "## Intermediates",
        "",
        "Every file produced by the pipeline lives inside `_pipeline/`. Useful for",
        "debugging and QC but not needed for downstream analysis. The internal",
        "structure (assembly/, blast/, porechop/, recombination/,",
        "pretelomeric_labels/, ...) mirrors the original pipeline layout.",
        "",
    ]
    (base_dir / "MANIFEST.md").write_text("\n".join(lines) + "\n")


def organize(base_name: str, strain: Optional[str] = None,
             results_dir: Path = Path("results"), verbose: bool = True) -> Dict[str, str]:
    base_dir = results_dir / base_name
    pipeline_dir = base_dir / "_pipeline"

    if not pipeline_dir.is_dir():
        if verbose:
            print(f"[organize] {pipeline_dir} does not exist yet — skipping")
        return {}

    if strain is None:
        strain = detect_strain(pipeline_dir)

    statuses: Dict[str, str] = {}
    present = []
    for name, target_template in FINALS.items():
        if strain is None and "{strain}" in target_template:
            statuses[name] = "missing"
            continue
        target_rel = "_pipeline/" + target_template.format(strain=strain, base=base_name)
        status = ensure_symlink(base_dir / name, target_rel)
        statuses[name] = status
        if status in ("created", "updated", "skipped"):
            present.append(name)
        if verbose and status not in ("skipped", "missing"):
            print(f"[organize] {name}: {status}")

    write_manifest(base_dir, present)
    return statuses


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1] if __doc__ else "")
    ap.add_argument("base_name", help="Project base name (directory under results/)")
    ap.add_argument("--strain", default=None, help="Strain ID (auto-detected if omitted)")
    ap.add_argument("--results-dir", default="results", type=Path)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    statuses = organize(
        base_name=args.base_name,
        strain=args.strain,
        results_dir=args.results_dir,
        verbose=not args.quiet,
    )
    blocked = [n for n, s in statuses.items() if s == "blocked"]
    if blocked:
        print(f"[organize] WARNING: real files exist at these paths; refusing to overwrite: {blocked}",
              file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
