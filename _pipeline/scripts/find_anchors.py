#!/usr/bin/env python3
"""
find_anchors.py

Find unique 5kb anchor sequences near each chromosome end.

Strategy:
  For each chromosome end (left/right), generate candidate 5kb windows
  sliding from the end toward the center in steps. BLAST all candidates
  against the full reference (excluding self-chromosome hits). Select the
  most distal window where no other chromosome has a hit with:
      pident >= identity_threshold  AND
      alignment_length / window_size >= coverage_threshold

  If no window passes, use the fallback: the window with the lowest
  maximum cross-chromosome match score.

Output:
  <output>_anchors.fasta   — FASTA of all anchor sequences
  <output>_anchors.tsv     — Summary table

Usage:
  python find_anchors.py reference.fasta -o anchors [options]
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------

def check_dependencies():
    for tool in ["makeblastdb", "blastn"]:
        if not shutil.which(tool):
            sys.exit(f"ERROR: '{tool}' not found in PATH. Install BLAST+.")
    try:
        from Bio import SeqIO          # noqa: F401
        from Bio.SeqRecord import SeqRecord  # noqa: F401
        from Bio.Seq import Seq        # noqa: F401
    except ImportError:
        sys.exit("ERROR: BioPython not installed. Run: pip install biopython")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Find unique 5kb anchor sequences near each chromosome end",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("reference",
                   help="Input FASTA reference")
    p.add_argument("-o", "--output", default="anchors",
                   help="Output file prefix")
    p.add_argument("--window-size", type=int, default=5040,
                   help="Anchor window size (bp). Must stay 5040: "
                        "filter_for_reads_with_anchors.py hardcodes that constant when "
                        "computing read_length_past_anchor on R arms.")
    p.add_argument("--step-size", type=int, default=500,
                   help="Sliding step size (bp)")
    p.add_argument("--max-slide", type=int, default=200_000,
                   help="Maximum distance to slide from the chromosome end (bp)")
    p.add_argument("--identity-threshold", type=float, default=85.0,
                   help="Minimum %% identity to consider a cross-chromosome hit disqualifying")
    p.add_argument("--coverage-threshold", type=float, default=35.0,
                   help="Minimum %% of window that must be covered by a hit to be disqualifying. "
                        "Calibrated on 6991: at 50 the derived anchors carried cross-matches up "
                        "to 2,380 bp (vs the curated set's 1,849 bp worst case); at 35 they carry "
                        "none, while anchor distances stay in the curated range (max 52.0 kb "
                        "vs the curated 53.5 kb).")
    p.add_argument("--exclude", nargs="+", default=["mito"],
                   help="Chromosome names to exclude (e.g. mitochondria)")
    p.add_argument("--threads", type=int, default=4,
                   help="BLAST threads")
    p.add_argument("--db-dir",
                   help="Directory for BLAST database (default: auto temp dir, deleted on exit)")
    p.add_argument("--keep-db", action="store_true",
                   help="Keep BLAST database after the run (only applies when --db-dir is set)")
    p.add_argument("--max-distance-from-end", type=int, default=60_000,
                   help="Fail if an anchor ends up further than this from its chromosome end. "
                        "Calibrated against the curated 6991 set, whose anchors span "
                        "1,839 bp (chr15L) to 53,546 bp (chr12R) from their ends.")
    p.add_argument("--allow-fallback", action="store_true",
                   help="Emit non-unique fallback anchors instead of failing. A fallback anchor "
                        "cross-matches another locus, and filter_for_reads_with_anchors.py DELETES "
                        "every read that matches two anchors, so this silently loses data.")
    p.add_argument("--base-reference-out",
                   help="Also write the companion <strain>_only_to_anchors.fasta that "
                        "create_ref.sh extends from: each chromosome truncated to "
                        "[L_anchor_start .. R_anchor_end].")
    return p.parse_args()


# ---------------------------------------------------------------------------
# BLAST helpers
# ---------------------------------------------------------------------------

def build_blast_db(fasta_path: str, db_path: str):
    print(f"Building BLAST database at {db_path} ...")
    r = subprocess.run(
        ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl", "-out", db_path],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        sys.exit(f"makeblastdb failed:\n{r.stderr}")
    print("  Done.\n")


def run_blast(query_fasta: str, db_path: str, threads: int) -> list[dict]:
    """
    BLAST a multi-FASTA query against the database.
    Returns a list of dicts: {qseqid, sseqid, pident, aln_len, qlen}
    """
    cmd = [
        "blastn",
        "-query",           query_fasta,
        "-db",              db_path,
        "-outfmt",          "6 qseqid sseqid pident length qlen sstart send",
        "-perc_identity",   "75",   # slightly below our threshold to capture near-misses
        "-num_threads",     str(threads),
        "-dust",            "no",   # do NOT mask — critical near telomeres / subtelomeric repeats
        "-task",            "blastn",
        "-max_target_seqs", "50",   # enough to cover all chromosomes many times over
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"blastn failed:\n{r.stderr}")

    hits = []
    for line in r.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        sstart, send = int(parts[5]), int(parts[6])
        hits.append({
            "qseqid": parts[0],
            "sseqid": parts[1],
            "pident": float(parts[2]),
            "aln_len": int(parts[3]),
            "qlen":    int(parts[4]),
            # normalise to 0-based half-open, orientation-independent
            "sbeg":    min(sstart, send) - 1,
            "send":    max(sstart, send),
        })
    return hits


# ---------------------------------------------------------------------------
# Window scoring
# ---------------------------------------------------------------------------

def evaluate_window(
    hits: list[dict],
    source_chrom: str,
    window_start: int,
    window_size: int,
    identity_thresh: float,
    coverage_thresh: float,
) -> tuple[bool, float]:
    """
    Evaluate whether a window is disqualified by similarity to any OTHER locus.

    Returns:
        (disqualified, max_match_score)

    max_match_score = max(pident × aln_coverage/100) across all non-self hits.
    A window is disqualified if ANY hit satisfies:
        pident >= identity_thresh  AND  aln_coverage >= coverage_thresh

    where aln_coverage = aln_len / window_size * 100.

    Only the window's own position counts as "self". Excluding the whole source
    chromosome would let a window that also matches the OTHER end of the same
    chromosome pass, and shared subtelomeric repeat is exactly where that happens.
    The chr12 rDNA array is the other case this catches.
    """
    max_score = 0.0
    disqualified = False
    win_beg, win_end = window_start, window_start + window_size

    for h in hits:
        # self = same sequence AND overlapping this window's own coordinates
        if h["sseqid"] == source_chrom and h["sbeg"] < win_end and h["send"] > win_beg:
            continue

        aln_cov = h["aln_len"] / window_size * 100.0
        # composite score: fraction of window covered at that identity
        score = h["pident"] * (aln_cov / 100.0)
        if score > max_score:
            max_score = score

        if h["pident"] >= identity_thresh and aln_cov >= coverage_thresh:
            disqualified = True

    return disqualified, max_score


# ---------------------------------------------------------------------------
# Core per-end logic
# ---------------------------------------------------------------------------

def find_anchor_for_end(
    chrom_name: str,
    chrom_seq: str,
    end: str,           # 'left' or 'right'
    db_path: str,
    args,
) -> dict | None:
    """
    Generate all candidate windows for one chromosome end, BLAST them all
    in a single batch, then pick the most distal passing window.
    """
    chrom_len = len(chrom_seq)
    window_size = args.window_size
    step_size = args.step_size

    # Cap max_slide so we never overlap into the other half of the chromosome
    max_slide = min(args.max_slide, chrom_len // 2 - window_size)
    if max_slide <= 0:
        max_slide = max(0, chrom_len - window_size)

    # --- Generate candidate windows (ordered: most distal first) ---
    windows = []   # list of (name, start, seq)
    for slide in range(0, max_slide + 1, step_size):
        start = slide if end == "left" else chrom_len - window_size - slide
        if start < 0 or start + window_size > chrom_len:
            break
        win_name = f"{chrom_name}_{end}_s{start}"
        win_seq  = chrom_seq[start : start + window_size]
        windows.append((win_name, start, win_seq))

    if not windows:
        print(f"    WARNING: no valid windows for {chrom_name} {end} end — skipping")
        return None

    # --- Write all windows to a temp FASTA and run BLAST once ---
    tmp = tempfile.NamedTemporaryFile(
        suffix=".fasta", mode="w", delete=False,
        prefix=f"anchor_{chrom_name}_{end}_",
    )
    for win_name, _, win_seq in windows:
        tmp.write(f">{win_name}\n{win_seq}\n")
    tmp.close()

    try:
        all_hits = run_blast(tmp.name, db_path, args.threads)
    finally:
        os.unlink(tmp.name)

    # --- Index hits by window name ---
    hits_by_window: dict[str, list[dict]] = {}
    for h in all_hits:
        hits_by_window.setdefault(h["qseqid"], []).append(h)

    # --- Find the most distal passing window; track best fallback ---
    best_fallback_idx   = 0
    best_fallback_score = float("inf")

    for idx, (win_name, start, win_seq) in enumerate(windows):
        disq, max_score = evaluate_window(
            hits_by_window.get(win_name, []),
            source_chrom=chrom_name,
            window_start=start,
            window_size=window_size,
            identity_thresh=args.identity_threshold,
            coverage_thresh=args.coverage_threshold,
        )

        if not disq:
            # Most distal passing window — return immediately
            dist = start if end == "left" else chrom_len - (start + window_size)
            return {
                "chrom":            chrom_name,
                "end":              end,
                "start":            start,
                "end_coord":        start + window_size,
                "seq":              win_seq,
                "max_match_pct":    max_score,
                "distance_from_end": dist,
                "fallback":         False,
                "n_windows_checked": idx + 1,
            }

        # Track the fallback (lowest max match seen so far)
        if max_score < best_fallback_score:
            best_fallback_score = max_score
            best_fallback_idx   = idx

    # --- No window passed threshold — use fallback ---
    win_name, start, win_seq = windows[best_fallback_idx]
    dist = start if end == "left" else chrom_len - (start + window_size)
    print(
        f"    WARNING: no window below threshold — "
        f"using fallback at pos={start} (max_match={best_fallback_score:.1f}%)"
    )
    return {
        "chrom":             chrom_name,
        "end":               end,
        "start":             start,
        "end_coord":         start + window_size,
        "seq":               win_seq,
        "max_match_pct":     best_fallback_score,
        "distance_from_end": dist,
        "fallback":          True,
        "n_windows_checked": len(windows),
    }


# ---------------------------------------------------------------------------
# Whole-set verification
# ---------------------------------------------------------------------------

def cross_match_lengths(fasta_path: str) -> list[tuple[str, str, float, int]]:
    """All-vs-all BLAST of the anchor set; returns non-self (q, s, pident, aln_len)."""
    cmd = [
        "blastn", "-query", fasta_path, "-subject", fasta_path,
        "-outfmt", "6 qseqid sseqid pident length", "-dust", "no", "-task", "blastn",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"blastn (self-comparison) failed:\n{r.stderr}")
    out = []
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 4 or p[0] == p[1]:
            continue
        out.append((p[0], p[1], float(p[2]), int(p[3])))
    return out


def verify_anchor_set(records, fasta_path: str, window_size: int) -> tuple[list, list]:
    """
    Check every constraint the rest of the pipeline silently assumes.
    Returns (problems, warnings); empty problems means the set is usable.

    On cross-matching: filter_for_reads_with_anchors.py:100 does
    drop_duplicates(subset=['read_id'], keep=False), so a read matching two anchors
    is DELETED. But line 86 only counts a hit when read_bp_used_for_match exceeds
    total_anchor_length/2, so the bar is shared sequence approaching that half-length,
    not mere shared k-mers. Calibration: the shipped, working 6991 set has a maximum
    cross-anchor alignment of 1,849 bp against a 2,520 bp filter — ~670 bp of headroom.
    A zero-shared-sequence bar would reject the known-good set.
    """
    problems, warnings = [], []
    expected = {f"chr{n}{arm}_anchor" for n in range(1, 17) for arm in ("L", "R")}
    got = [r.id for r in records]

    if len(got) != 32:
        problems.append(f"expected 32 anchors, got {len(got)}")
    dupes = {n for n in got if got.count(n) > 1}
    if dupes:
        problems.append(f"duplicate anchor ids: {sorted(dupes)}")
    if expected - set(got):
        problems.append(f"missing anchor ids: {sorted(expected - set(got))}")
    if set(got) - expected:
        problems.append(f"unexpected anchor ids: {sorted(set(got) - expected)}")

    for r in records:
        if len(r.seq) != window_size:
            problems.append(f"{r.id} is {len(r.seq)} bp, must be exactly {window_size}")
        n_count = str(r.seq).upper().count("N")
        if n_count:
            problems.append(f"{r.id} contains {n_count} N bases")

    fail_at = window_size // 2          # the downstream significance filter
    warn_at = int(fail_at * 0.8)        # headroom comparable to the shipped set
    worst = {}
    for q, s, pident, aln_len in cross_match_lengths(fasta_path):
        key = tuple(sorted((q, s)))
        if aln_len > worst.get(key, (0, 0))[0]:
            worst[key] = (aln_len, pident)
    for (a, b), (aln_len, pident) in sorted(worst.items(), key=lambda x: -x[1][0]):
        if aln_len >= fail_at:
            problems.append(
                f"{a} and {b} share a {aln_len} bp alignment at {pident:.1f}% "
                f"(>= the {fail_at} bp significance filter — reads matching both "
                "would be discarded)"
            )
        elif aln_len >= warn_at:
            warnings.append(
                f"{a} and {b} share {aln_len} bp at {pident:.1f}% "
                f"(below the {fail_at} bp filter, but less headroom than the "
                "shipped set's 1,849 bp worst case)"
            )
    return problems, warnings


def write_base_reference(path, chroms, summary_rows, anchor_records, window_size):
    """
    Write the companion <strain>_only_to_anchors.fasta.

    create_ref.sh:37 extends outward from a reference truncated at the anchors --
    extend_reference_multi relies on "the base reference is truncated at the anchors so
    the read's true distal home is largely absent from it"
    (subtelomere_reference_pipeline_utils.py:509). Verified against the shipped 6991 file:
    each contig starts at its L anchor's first base and ends at its R anchor's last.

    Contigs are named plainly chr1..chr16 because contig_for_arm()
    (subtelomere_reference_pipeline_utils.py:575-580) resolves an arm by digit string and
    raises unless exactly one contig matches.
    """
    by_chrom = {}
    for r in summary_rows:
        by_chrom.setdefault(r["chrom"], {})[r["end"]] = r

    anchors = {rec.id: str(rec.seq) for rec in anchor_records}
    problems = []
    written = []
    with open(path, "w") as fh:
        for n in range(1, 17):
            chrom = f"chr{n}"
            ends = by_chrom.get(chrom)
            if not ends or "left" not in ends or "right" not in ends:
                problems.append(f"{chrom}: missing an anchor, cannot truncate")
                continue
            lo = ends["left"]["start"]           # first base of the L anchor
            hi = ends["right"]["end_coord"]      # one past the last base of the R anchor
            seq = chroms[chrom].seq if hasattr(chroms[chrom], "seq") else chroms[chrom]
            sub = str(seq)[lo:hi]

            # self-check: the slice must begin with the L anchor and end with the R anchor
            if sub[:window_size] != anchors.get(f"{chrom}L_anchor"):
                problems.append(f"{chrom}: first {window_size} bp != {chrom}L_anchor")
            if sub[-window_size:] != anchors.get(f"{chrom}R_anchor"):
                problems.append(f"{chrom}: last {window_size} bp != {chrom}R_anchor")

            fh.write(f">{chrom}\n")
            for i in range(0, len(sub), 60):
                fh.write(sub[i:i + 60] + "\n")
            written.append((chrom, len(sub)))

    if problems:
        print(f"\n  ✗ base reference problems:")
        for p in problems:
            print(f"      - {p}")
        sys.exit(1)

    total = sum(l for _, l in written)
    print(f"  ✓ base reference: {len(written)} chromosomes, {total:,} bp  →  {path}")
    subprocess.run(["samtools", "faidx", path], check=False,
                   capture_output=True, text=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    check_dependencies()
    args = parse_args()

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    # --- Load reference ---
    print(f"Loading reference: {args.reference}")
    records = {rec.id: rec for rec in SeqIO.parse(args.reference, "fasta")}

    chroms = {k: v for k, v in records.items() if k not in args.exclude}
    excluded_names = [k for k in records if k in args.exclude]
    print(f"  {len(chroms)} chromosome(s) to process  |  excluded: {excluded_names}")
    print(f"  Expecting {len(chroms) * 2} anchors  (2 per chromosome)\n")

    # --- Build BLAST database ---
    cleanup_db = False
    if args.db_dir:
        db_dir = args.db_dir
        os.makedirs(db_dir, exist_ok=True)
    else:
        db_dir = tempfile.mkdtemp(prefix="anchor_blast_")
        cleanup_db = not args.keep_db

    db_path = os.path.join(db_dir, "ref_db")

    try:
        build_blast_db(args.reference, db_path)

        anchor_records = []
        summary_rows   = []

        for chrom_name in sorted(chroms.keys()):
            chrom_seq = str(chroms[chrom_name].seq)
            chrom_len = len(chrom_seq)
            print(f"[{chrom_name}]  length = {chrom_len:,} bp")

            # Infer which end to search from the sequence name:
            #   chrXL → telomere is at position 0        → search "left"
            #   chrXR → telomere is at the rightmost pos → search "right"
            # Fall back to searching both ends for unlabeled sequences.
            name_upper = chrom_name.upper()
            if name_upper.endswith("L"):
                ends_to_search = ("left",)
            elif name_upper.endswith("R"):
                ends_to_search = ("right",)
            else:
                ends_to_search = ("left", "right")

            for end in ends_to_search:
                print(f"  {end} end:")
                result = find_anchor_for_end(chrom_name, chrom_seq, end, db_path, args)

                if result is None:
                    print(f"    ERROR — could not find anchor\n")
                    continue

                flag = "[FALLBACK]" if result["fallback"] else "[OK]     "
                print(
                    f"    {flag}  pos {result['start']:,}–{result['end_coord']:,}  "
                    f"dist_from_end={result['distance_from_end']:,} bp  "
                    f"max_cross_match={result['max_match_pct']:.1f}%  "
                    f"({result['n_windows_checked']} windows checked)"
                )

                # Naming must match what filter_for_reads_with_anchors.py and
                # split_and_label_all_reads_include_anchor.py expect: chr1L_anchor etc.
                # A contig named chr1 yields two anchors, so the arm letter has to be
                # appended here or both records collide on the same id.
                arm = "L" if end == "left" else "R"
                anchor_id = (chrom_name if name_upper.endswith(("L", "R"))
                             else f"{chrom_name}{arm}") + "_anchor"
                rec = SeqRecord(
                    Seq(result["seq"]),
                    id=anchor_id,
                    description=(
                        f"pos={result['start']}-{result['end_coord']} "
                        f"dist_from_end={result['distance_from_end']} "
                        f"max_cross_chrom_match={result['max_match_pct']:.1f}% "
                        f"fallback={result['fallback']}"
                    ),
                )
                anchor_records.append(rec)
                summary_rows.append(result)

            print()

        # --- Write outputs ---
        out_fasta = args.output + "_anchors.fasta"
        out_tsv   = args.output + "_anchors.tsv"

        SeqIO.write(anchor_records, out_fasta, "fasta")

        with open(out_tsv, "w") as f:
            f.write(
                "chrom\tend\tstart\tend_coord\t"
                "max_match_pct\tdistance_from_end\tfallback\tn_windows_checked\n"
            )
            for r in summary_rows:
                f.write(
                    f"{r['chrom']}\t{r['end']}\t{r['start']}\t{r['end_coord']}\t"
                    f"{r['max_match_pct']:.2f}\t{r['distance_from_end']}\t"
                    f"{r['fallback']}\t{r['n_windows_checked']}\n"
                )

        n_fallback = sum(1 for r in summary_rows if r["fallback"])
        print("=" * 60)
        print(f"Anchors written : {len(anchor_records)}  →  {out_fasta}")
        print(f"Summary table   : {out_tsv}")
        print(f"Fallbacks used  : {n_fallback} / {len(summary_rows)}")

        # --- Verification: every constraint the rest of the pipeline assumes ---
        print("\nVerifying anchor set ...")
        problems, warnings = verify_anchor_set(anchor_records, out_fasta, args.window_size)

        too_far = [r for r in summary_rows
                   if r["distance_from_end"] > args.max_distance_from_end]
        for r in too_far:
            problems.append(
                f"{r['chrom']} {r['end']}: anchor sits {r['distance_from_end']:,} bp "
                f"from the chromosome end (limit {args.max_distance_from_end:,})"
            )

        fallbacks = [r for r in summary_rows if r["fallback"]]
        if fallbacks and not args.allow_fallback:
            for r in fallbacks:
                problems.append(
                    f"{r['chrom']} {r['end']}: no unique window found "
                    f"(best cross-match {r['max_match_pct']:.1f}%)"
                )

        for w in warnings:
            print(f"  ! {w}")

        if problems:
            print(f"\n  ✗ {len(problems)} problem(s):")
            for p in problems:
                print(f"      - {p}")
            print(
                "\nThis anchor set is NOT safe to use. Options: raise --max-slide, "
                "adjust --identity-threshold / --coverage-threshold, or improve the "
                "assembly at the affected ends."
            )
            sys.exit(1)

        print(f"  ✓ 32 anchors, all {args.window_size} bp, no disqualifying "
              "cross-matches, no fallbacks")

        # BLAST db in place — the FASTA path is used directly as the blastn -db prefix
        build_blast_db(out_fasta, out_fasta)
        print(f"  ✓ BLAST database built at {out_fasta}")

        # --- Companion base reference for create_ref.sh ---
        if args.base_reference_out:
            write_base_reference(args.base_reference_out, chroms, summary_rows,
                                 anchor_records, args.window_size)

    finally:
        if cleanup_db:
            shutil.rmtree(db_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
