#!/usr/bin/env python3
"""Check that every anchor in a set actually recruits reads, and that none cross-match.

This is the test that would have caught the 6212 failure immediately: 11 of 32 anchors
recruited ZERO reads, which no existing check reports. filter_for_reads_with_anchors.py
prints anchor value-counts, but an anchor with no reads simply doesn't appear.

Reads the raw per-read anchor BLAST table, i.e.
  results/<base>/_pipeline/blast/all_matches_<base>_blasted_<anchor_set>.tsv

Usage: check_anchor_recruitment.py --all-matches <tsv> [--min-reads 50]
                                   [--min-identity 98.5] [--strict]
"""
import argparse
import csv
import statistics
import sys
from collections import defaultdict

EXPECTED = [f"chr{n}{arm}_anchor" for n in range(1, 17) for arm in ("L", "R")]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--all-matches", required=True)
    p.add_argument("--min-reads", type=int, default=50)
    p.add_argument("--min-identity", type=float, default=98.5)
    p.add_argument("--strict", action="store_true",
                   help="Exit non-zero if any criterion fails")
    a = p.parse_args()

    per_anchor = defaultdict(list)
    reads_to_anchors = defaultdict(set)
    with open(a.all_matches) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            anchor = row["anchor_name"]
            per_anchor[anchor].append(float(row["pident"]))
            reads_to_anchors[row["read_id"]].add(anchor)

    print("=" * 72)
    print("ANCHOR RECRUITMENT")
    print("=" * 72)
    print(f"{'anchor':<20} {'reads':>8} {'median id':>10} {'min id':>8}")

    dead, low_reads, low_id = [], [], []
    for anchor in EXPECTED:
        ids = per_anchor.get(anchor, [])
        if not ids:
            dead.append(anchor)
            print(f"{anchor:<20} {0:>8} {'-':>10} {'-':>8}   <-- NO READS")
            continue
        med, lo = statistics.median(ids), min(ids)
        flag = ""
        if len(ids) < a.min_reads:
            low_reads.append(anchor)
            flag = "   <-- few reads"
        if med < a.min_identity:
            low_id.append(anchor)
            flag += "   <-- low identity"
        print(f"{anchor:<20} {len(ids):>8} {med:>10.2f} {lo:>8.2f}{flag}")

    # reads matching >1 anchor are silently discarded downstream
    # (filter_for_reads_with_anchors.py:100, drop_duplicates(keep=False))
    multi = {r: s for r, s in reads_to_anchors.items() if len(s) > 1}

    unexpected = sorted(set(per_anchor) - set(EXPECTED))

    print()
    print("-" * 72)
    print(f"anchors with zero reads      : {len(dead)}" +
          (f"  {dead}" if dead else "   OK"))
    print(f"anchors below {a.min_reads} reads      : {len(low_reads)}" +
          (f"  {low_reads}" if low_reads else "   OK"))
    print(f"anchors below {a.min_identity}% identity : {len(low_id)}" +
          (f"  {low_id}" if low_id else "   OK"))
    print(f"reads matching >1 anchor     : {len(multi)}" +
          ("   (these are DISCARDED downstream)" if multi else "   OK"))
    if unexpected:
        print(f"unexpected anchor names      : {unexpected}")

    failed = bool(dead or multi or unexpected)
    print()
    print("RESULT:", "FAIL" if failed else "PASS")
    if a.strict and failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
