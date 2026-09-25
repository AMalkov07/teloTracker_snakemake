#!/usr/bin/env python3
"""QC a de novo assembly before deriving anchors from it.

The decisive metric is TELOMERE COMPLETENESS at contig ends. Anchors are chosen as the
most distal unique window, so if contig ends stop short of the telomere the anchors land
further inward than the reference set's and stop being "just proximal to the subtelomeric
repeats". Everything else here is context for that judgement.

Usage: assembly_qc.py --assembly <fasta> [--assembly-info <flye assembly_info.txt>]
                      [--s288c <fasta>] [--threads N] --out-prefix <prefix>
"""
import argparse
import os
import re
import subprocess
import sys

TELO_WINDOW = 2000      # how far in from each contig end to look
MIN_TRACT = 24          # shortest run that counts as a real telomeric tract

# S. cerevisiae telomeric repeat is TG(1-3) read toward the chromosome end; the
# complement C(1-3)A is what appears at a left/5' end on the plus strand.
TG_RE = re.compile(r"(?:T{1,3}G{1,3}){4,}", re.I)
CA_RE = re.compile(r"(?:C{1,3}A{1,3}){4,}", re.I)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def best_tract(seq, pattern):
    """Longest match of `pattern` in seq -> (length, start, end) or (0, None, None)."""
    best = (0, None, None)
    for m in pattern.finditer(seq):
        if len(m.group()) > best[0]:
            best = (len(m.group()), m.start(), m.end())
    return best


def telomere_status(seq):
    """Per-contig telomere assessment at both ends."""
    left, right = seq[:TELO_WINDOW], seq[-TELO_WINDOW:]
    # a left end terminates in C(1-3)A on the plus strand; a right end in TG(1-3)
    l_len, l_start, _ = best_tract(left, CA_RE)
    r_len, _, r_end = best_tract(right, TG_RE)
    return {
        "left_tract_bp": l_len,
        # distance from the very first base to where the tract begins
        "left_gap_bp": l_start if l_len >= MIN_TRACT else None,
        "left_telomeric": l_len >= MIN_TRACT,
        "right_tract_bp": r_len,
        "right_gap_bp": (TELO_WINDOW - r_end) if r_len >= MIN_TRACT else None,
        "right_telomeric": r_len >= MIN_TRACT,
    }


def n50(lengths):
    s = sorted(lengths, reverse=True)
    half, run = sum(s) / 2, 0
    for x in s:
        run += x
        if run >= half:
            return x
    return 0


def s288c_coverage(assembly, s288c, threads, out_paf):
    """minimap2 assembly -> S288C; returns {target: aligned_bases} and the PAF path."""
    cmd = ["minimap2", "-x", "asm5", "-c", "--secondary=no", "-t", str(threads),
           s288c, assembly]
    with open(out_paf, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        print(f"  WARNING: minimap2 failed, skipping S288C comparison:\n{r.stderr[:500]}",
              file=sys.stderr)
        return None
    per_target, target_len = {}, {}
    for line in open(out_paf):
        f = line.split("\t")
        if len(f) < 11:
            continue
        tgt, tlen, matches = f[5], int(f[6]), int(f[9])
        per_target[tgt] = per_target.get(tgt, 0) + matches
        target_len[tgt] = tlen
    return per_target, target_len


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--assembly", required=True)
    p.add_argument("--assembly-info")
    p.add_argument("--s288c")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--expected-size", type=int, default=12_070_000)
    p.add_argument("--out-prefix", required=True)
    a = p.parse_args()

    contigs = read_fasta(a.assembly)
    lengths = {k: len(v) for k, v in contigs.items()}
    total = sum(lengths.values())

    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
    report = open(a.out_prefix + "_assembly_qc.txt", "w")

    def out(s=""):
        print(s)
        report.write(s + "\n")

    out("=" * 72)
    out("ASSEMBLY QC")
    out("=" * 72)
    out(f"assembly      : {a.assembly}")
    out(f"contigs       : {len(contigs)}")
    out(f"total length  : {total:,} bp  ({100.0 * total / a.expected_size:.1f}% of "
        f"{a.expected_size:,} expected)")
    out(f"N50           : {n50(lengths.values()):,} bp")
    out(f"largest       : {max(lengths.values()):,} bp")
    out()

    # ---- telomere completeness: the gate ----
    out("-" * 72)
    out("TELOMERE COMPLETENESS AT CONTIG ENDS")
    out("  A contig end is 'telomeric' if a >=%d bp TG(1-3)/C(1-3)A tract sits within "
        "the terminal %d bp." % (MIN_TRACT, TELO_WINDOW))
    out("  gap = bases between the contig terminus and the tract (0 is ideal).")
    out("-" * 72)
    out(f"{'contig':<24} {'len':>12}  {'L tract':>8} {'L gap':>7}  "
        f"{'R tract':>8} {'R gap':>7}")
    n_telo = 0
    for name in sorted(contigs, key=lambda k: -lengths[k]):
        t = telomere_status(contigs[name])
        n_telo += int(t["left_telomeric"]) + int(t["right_telomeric"])
        lg = "-" if t["left_gap_bp"] is None else f"{t['left_gap_bp']:,}"
        rg = "-" if t["right_gap_bp"] is None else f"{t['right_gap_bp']:,}"
        out(f"{name:<24} {lengths[name]:>12,}  {t['left_tract_bp']:>8} {lg:>7}  "
            f"{t['right_tract_bp']:>8} {rg:>7}")
    n_ends = 2 * len(contigs)
    out()
    out(f"telomeric ends: {n_telo} / {n_ends}  ({100.0 * n_telo / n_ends:.0f}%)")
    out("  For reference, a set of 32 chromosome ends is the target. Ends without a")
    out("  telomeric tract will place their anchor further inward than the 6991 set.")
    out()

    # ---- S288C comparison ----
    if a.s288c:
        out("-" * 72)
        out("ALIGNMENT TO S288C")
        out("-" * 72)
        res = s288c_coverage(a.assembly, a.s288c, a.threads, a.out_prefix + "_vs_s288c.paf")
        if res:
            per_target, target_len = res
            covered = 0
            for tgt in sorted(per_target, key=lambda t: -per_target[t]):
                frac = 100.0 * per_target[tgt] / target_len[tgt]
                covered += 1
                out(f"  {tgt:<16} {per_target[tgt]:>10,} matching bp  "
                    f"({frac:5.1f}% of {target_len[tgt]:,})")
            out(f"\n  S288C sequences hit: {covered}")
        out()

    # ---- Flye per-contig coverage ----
    if a.assembly_info and os.path.exists(a.assembly_info):
        out("-" * 72)
        out("FLYE assembly_info.txt")
        out("-" * 72)
        for line in open(a.assembly_info):
            out("  " + line.rstrip())
        out()

    report.close()
    print(f"\nwritten {a.out_prefix}_assembly_qc.txt")


if __name__ == "__main__":
    main()
