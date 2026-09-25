#!/usr/bin/env python3
"""Assign de novo assembly contigs to chr1..chr16 and put them on the plus strand.

Flye emits `contig_1`, `contig_2`... in arbitrary orientation. Downstream code needs
records named exactly `chr1`..`chr16` on the forward strand:

  * find_anchors.py slices windows forward and names them chr<N><L|R>_anchor
  * subtelomere_reference_pipeline_utils.py:575-580 (contig_for_arm) resolves an arm
    like "12R" by matching the digit string and RAISES unless exactly one contig matches

Assignment is by whole-genome alignment to S288C. Deliberately does NOT splice any S288C
sequence into the output -- strain divergence from the reference is the entire reason
this tooling exists.

Usage: assign_contigs_to_chromosomes.py --assembly <fasta> --s288c <fasta>
                                        --out <fasta> [--threads N]
                                        [--min-fraction 0.7] [--allow-gapped]
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict

# RefSeq accession -> arabic chromosome. Explicit, rather than parsing Roman numerals
# out of free-text descriptions.
S288C_CHROMS = {
    "NC_001133.9": "chr1",  "NC_001134.8": "chr2",  "NC_001135.5": "chr3",
    "NC_001136.10": "chr4", "NC_001137.3": "chr5",  "NC_001138.5": "chr6",
    "NC_001139.9": "chr7",  "NC_001140.6": "chr8",  "NC_001141.2": "chr9",
    "NC_001142.9": "chr10", "NC_001143.9": "chr11", "NC_001144.5": "chr12",
    "NC_001145.3": "chr13", "NC_001146.8": "chr14", "NC_001147.6": "chr15",
    "NC_001148.4": "chr16",
}
MITO = "NC_001224.1"

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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


def revcomp(s):
    return s.translate(COMP)[::-1]


def align(assembly, s288c, threads, paf_path):
    cmd = ["minimap2", "-x", "asm5", "-c", "--secondary=no", "-t", str(threads),
           s288c, assembly]
    with open(paf_path, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit(f"minimap2 failed:\n{r.stderr}")
    rows = []
    for line in open(paf_path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        rows.append({
            "qname": f[0], "qlen": int(f[1]), "qstart": int(f[2]), "qend": int(f[3]),
            "strand": f[4], "tname": f[5], "tlen": int(f[6]),
            "tstart": int(f[7]), "tend": int(f[8]), "matches": int(f[9]),
        })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--assembly", required=True)
    p.add_argument("--s288c", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--min-fraction", type=float, default=0.7,
                   help="Fraction of a contig's aligned bases that must land on its "
                        "assigned chromosome")
    p.add_argument("--allow-gapped", action="store_true",
                   help="Permit joining two contigs for one chromosome with an N gap. "
                        "Untested downstream -- prefer improving the assembly.")
    p.add_argument("--min-contig", type=int, default=20_000,
                   help="Ignore contigs shorter than this (plasmids, debris)")
    a = p.parse_args()

    contigs = read_fasta(a.assembly)
    paf = os.path.splitext(a.out)[0] + "_vs_s288c.paf"
    rows = align(a.assembly, a.s288c, a.threads, paf)

    # ---- per contig: matching bases per target, and strand vote ----
    per_contig = defaultdict(lambda: defaultdict(int))
    strand_vote = defaultdict(lambda: defaultdict(int))
    for r in rows:
        chrom = S288C_CHROMS.get(r["tname"])
        if chrom is None:
            continue  # mito, 2-micron, unplaced
        per_contig[r["qname"]][chrom] += r["matches"]
        strand_vote[r["qname"]][r["strand"]] += r["matches"]

    assignments, problems = {}, []
    for cname, seq in contigs.items():
        if len(seq) < a.min_contig:
            print(f"  skip {cname} ({len(seq):,} bp < --min-contig)")
            continue
        targets = per_contig.get(cname)
        if not targets:
            problems.append(f"{cname} ({len(seq):,} bp) has no alignment to S288C")
            continue
        total = sum(targets.values())
        chrom, best = max(targets.items(), key=lambda kv: kv[1])
        frac = best / total
        if frac < a.min_fraction:
            problems.append(
                f"{cname} is ambiguous: only {frac:.0%} of aligned bases on {chrom} "
                f"(next: {sorted(targets.items(), key=lambda kv: -kv[1])[1:3]})")
            continue
        strand = max(strand_vote[cname].items(), key=lambda kv: kv[1])[0]
        # position along the target, for ordering when a chromosome is fragmented
        tstarts = [r["tstart"] for r in rows
                   if r["qname"] == cname and S288C_CHROMS.get(r["tname"]) == chrom]
        assignments.setdefault(chrom, []).append({
            "contig": cname, "strand": strand, "frac": frac,
            "len": len(seq), "tstart": min(tstarts) if tstarts else 0,
        })

    # ---- report and build output ----
    print("\n" + "=" * 72)
    print("CONTIG -> CHROMOSOME ASSIGNMENT")
    print("=" * 72)
    out_records = {}
    for n in range(1, 17):
        chrom = f"chr{n}"
        parts = sorted(assignments.get(chrom, []), key=lambda d: d["tstart"])
        if not parts:
            problems.append(f"{chrom}: no contig assigned")
            print(f"  {chrom:<7} MISSING")
            continue
        desc = ", ".join(f"{p['contig']}({p['strand']},{p['len']:,}bp,{p['frac']:.0%})"
                         for p in parts)
        print(f"  {chrom:<7} {len(parts)} contig(s): {desc}")

        seqs = [contigs[p["contig"]] if p["strand"] == "+" else revcomp(contigs[p["contig"]])
                for p in parts]
        if len(seqs) == 1:
            out_records[chrom] = seqs[0]
        elif a.allow_gapped:
            out_records[chrom] = ("N" * 100).join(seqs)
            problems.append(
                f"{chrom}: joined {len(seqs)} contigs with an N gap (--allow-gapped). "
                "The N run will propagate into the anchors' base reference.")
        else:
            problems.append(
                f"{chrom}: fragmented into {len(parts)} contigs. Re-run with "
                "--allow-gapped to join them with an N gap, or improve the assembly.")

    if problems:
        print("\n" + "-" * 72)
        print(f"{len(problems)} PROBLEM(S):")
        for pr in problems:
            print(f"  - {pr}")

    if len(out_records) != 16:
        print(f"\nERROR: {len(out_records)}/16 chromosomes resolved; not writing output.")
        sys.exit(1)

    with open(a.out, "w") as fh:
        for n in range(1, 17):
            chrom = f"chr{n}"
            fh.write(f">{chrom}\n")
            s = out_records[chrom]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    print(f"\nwritten {a.out}  (16 chromosomes, plus strand)")
    if problems:
        sys.exit(1)


if __name__ == "__main__":
    main()
