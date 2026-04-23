#!/usr/bin/env python3
"""
Cluster per-chr-end x-element sequences so that chr_ends whose x-elements
are indistinguishable share a cluster ID.

Mirrors scripts/cluster_yprimes_paper_method.py — same hierarchical clustering,
same 97% identity threshold, same homopolymer condensation — with two
differences:

1. Header format: `>X_Element_chrXL;chrYR#IDn_Color` (semicolon-joined chr_end
   members, like the Y′ library's `origin` body).
2. Similarity formula: `pident × alen / min(qlen, slen)` ("coverage of the
   shorter sequence") instead of the Y′ paper formula that penalizes the
   longer sequence's unaligned tail. Necessary because strain-specific
   x-elements can legitimately differ in total length while sharing an
   identical core — e.g. in 6991, chr16L (608 bp) is a truncation of
   chr15R (697 bp) with 99.83% identity over the overlap.

Usage:
    python scripts/cluster_x_elements.py <input_fasta> \
        --output-fasta <clustered_fasta> \
        --output-dir <workdir>
"""
from __future__ import annotations

import argparse
import os
import re
import sys

# Reuse helpers from the Y′ clustering script
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import numpy as np
import pandas as pd

from cluster_yprimes_paper_method import (
    ID_COLOR_NAMES,
    assign_ids_and_colors,
    condense_homopolymers,
    deduplicate_yprimes as _dedup_generic,
    find_clusters_threshold,
    load_fasta,
    run_blast,
    write_fasta,
)


CHR_END_RE = re.compile(r'(chr\d+[LR])')


def calc_similarity_matrix_coverage(blast_df, seq_names):
    """
    Like cluster_yprimes_paper_method.calc_similarity_matrix, but uses the
    coverage-of-shorter-sequence formula:

        sim = pident × alignment_length / min(qlen, slen)

    Rationale: when one sequence is a truncation of the other, BLAST still
    aligns the overlap at high identity; the "unaligned tail" term in the
    paper formula penalizes that (drops 99.83% to 86.94%) and leaves the
    pair unclustered even though biologically they're the same element.
    """
    n = len(seq_names)
    name_to_idx = {}
    for i, name in enumerate(seq_names):
        name_to_idx[name] = i
        first_word = name.split()[0]
        name_to_idx.setdefault(first_word, i)

    sim_matrix = np.zeros((n, n))
    np.fill_diagonal(sim_matrix, 100.0)
    if blast_df.empty:
        return sim_matrix

    for (qseqid, sseqid), group in blast_df.groupby(['qseqid', 'sseqid']):
        if qseqid == sseqid:
            continue
        qi = name_to_idx.get(qseqid)
        si = name_to_idx.get(sseqid)
        if qi is None or si is None:
            continue

        best = group.loc[group['length'].idxmax()]
        pident = float(best['pident'])
        aln_len = int(best['length'])
        qlen = int(best['qlen'])
        slen = int(best['slen'])

        denom = max(min(qlen, slen), 1)
        similarity = (pident / 100.0) * aln_len / denom * 100.0
        similarity = min(similarity, 100.0)

        sim_matrix[qi, si] = max(sim_matrix[qi, si], similarity)
        sim_matrix[si, qi] = max(sim_matrix[si, qi], similarity)

    return sim_matrix


def header_chr_end(header):
    """Extract the chr_end (chr1L / chr12R / ...) from an x-element header
    like `chr10L_x_ends` or the semicolon-joined body `chr9L;chr10L`."""
    m = CHR_END_RE.search(header.split()[0])
    return m.group(1) if m else header.split('_')[0]


def rewrite_x_element_header(members, id_str, color):
    """Build a clustered header.

    members is a list of chr_end strings (e.g. ['chr9L', 'chr10L']).
    Output: `X_Element_chr9L;chr10L#IDn_Color` (sorted alphabetically).
    """
    body = ';'.join(sorted(set(members)))
    return f'X_Element_{body}#{id_str}_{color}'


def write_clustered_fasta_x(original_entries, deduped_names, cluster_labels,
                            original_to_rep, output_path):
    """One FASTA entry per cluster. Body sequence is the representative's
    (uncondensed) sequence; the header lists every chr_end member and tags
    the cluster ID + color.
    """
    id_assignments = assign_ids_and_colors(cluster_labels, deduped_names)

    # Map deduped representative name -> (id, color)
    rep_to_id = {
        deduped_names[i]: id_assignments[i] for i in range(len(deduped_names))
    }
    # Map cluster_label (1-based from scipy) -> (id, color)
    label_to_id = {}
    for i, cl in enumerate(cluster_labels):
        label_to_id.setdefault(cl, id_assignments[i])

    # Group every ORIGINAL chr_end by its cluster label
    from collections import defaultdict
    cluster_members = defaultdict(list)   # cluster_label -> [chr_end, ...]
    cluster_rep_seq = {}                  # cluster_label -> (header, seq) of rep

    rep_to_label = {deduped_names[i]: cl for i, cl in enumerate(cluster_labels)}
    rep_to_entry = {h: (h, s) for h, s in original_entries}

    for orig_header, _ in original_entries:
        rep_name = original_to_rep.get(orig_header, orig_header)
        cl = rep_to_label.get(rep_name)
        if cl is None:
            continue
        cluster_members[cl].append(header_chr_end(orig_header))
        if cl not in cluster_rep_seq:
            cluster_rep_seq[cl] = rep_to_entry.get(rep_name, (orig_header, ''))

    output_entries = []
    for cl in sorted(cluster_members.keys()):
        id_str, color = label_to_id.get(cl, ('Unknown', 'Gray'))
        members = cluster_members[cl]
        _, rep_seq = cluster_rep_seq[cl]
        new_header = rewrite_x_element_header(members, id_str, color)
        output_entries.append((new_header, rep_seq))

    write_fasta(output_entries, output_path)
    print(f"  Clustered x-element FASTA: {output_path} ({len(output_entries)} clusters)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1] if __doc__ else '')
    p.add_argument('input_fasta',
                   help='Per-chr-end x-element FASTA (output of make_x_element_sequences.py)')
    p.add_argument('--output-fasta', required=True,
                   help='Path to write the clustered x-element FASTA')
    p.add_argument('--output-dir', default='.',
                   help='Directory for cluster-assignment TSV, similarity matrix, etc.')
    p.add_argument('--identity-threshold', type=float, default=97.0)
    p.add_argument('--dedup-threshold', type=float, default=99.9)
    p.add_argument('--linkage', default='complete',
                   choices=['average', 'complete', 'single'])
    p.add_argument('--no-condense-homopolymers', action='store_true')
    p.add_argument('--threads', type=int, default=4)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print("X-element clustering")
    print("=" * 70)
    print(f"  Input FASTA:         {args.input_fasta}")
    print(f"  Output FASTA:        {args.output_fasta}")
    print(f"  Working dir:         {args.output_dir}")
    print(f"  Similarity formula:  pident * alen / min(qlen, slen)  (coverage-of-shorter)")
    print(f"  Identity threshold:  {args.identity_threshold}%")
    print(f"  Dedup threshold:     {args.dedup_threshold}%")
    print(f"  Linkage:             {args.linkage}")

    entries = load_fasta(args.input_fasta)
    print(f"\nLoaded {len(entries)} x-element sequences")

    # Step 1: Homopolymer condensation
    if not args.no_condense_homopolymers:
        print("\n--- Step 1: Condense homopolymers >4bp ---")
        condensed_entries = [(h, condense_homopolymers(s.upper())) for h, s in entries]
        n_condensed = sum(1 for (_, s), (_, cs) in zip(entries, condensed_entries)
                          if len(s) != len(cs))
        print(f"  {n_condensed}/{len(entries)} sequences had homopolymers condensed")
    else:
        print("\n--- Step 1: Homopolymer condensation SKIPPED ---")
        condensed_entries = [(h, s.upper()) for h, s in entries]

    working_fasta = os.path.join(args.output_dir, 'x_elements_working.fasta')
    write_fasta(condensed_entries, working_fasta)

    seq_names = [h for h, _ in condensed_entries]

    # Step 2: All-vs-all BLAST + coverage-of-shorter similarity
    print("\n--- Step 2: All-vs-all BLAST ---")
    blast_df = run_blast(working_fasta, threads=args.threads)
    print(f"  BLAST produced {len(blast_df)} hits")
    sim_matrix = calc_similarity_matrix_coverage(blast_df, seq_names)

    n = len(seq_names)
    off_diag = sim_matrix[np.triu_indices(n, k=1)]
    if len(off_diag):
        print(f"  Similarity range: {off_diag.min():.1f}% - {off_diag.max():.1f}% "
              f"(mean {off_diag.mean():.1f}%)")

    # Step 3: Deduplicate near-identical pairs (default ≥99.9% → merged)
    original_to_rep = {name: name for name in seq_names}
    if args.dedup_threshold < 100.0:
        print(f"\n--- Step 3: Deduplicate at {args.dedup_threshold}% similarity ---")
        deduped_entries, deduped_names, original_to_rep = _dedup_generic(
            condensed_entries, sim_matrix, seq_names, threshold=args.dedup_threshold)
        print(f"  {len(condensed_entries)} sequences -> {len(deduped_entries)} unique variants")

        if len(deduped_entries) < len(condensed_entries):
            deduped_fasta = os.path.join(args.output_dir, 'x_elements_deduped.fasta')
            write_fasta(deduped_entries, deduped_fasta)
            blast_df2 = run_blast(deduped_fasta, threads=args.threads)
            sim_matrix = calc_similarity_matrix_coverage(blast_df2, deduped_names)
    else:
        print("\n--- Step 3: Deduplication SKIPPED (threshold=100%) ---")
        deduped_entries = condensed_entries
        deduped_names = seq_names

    # Step 4: Hierarchical clustering at the identity threshold
    distance_matrix = 100.0 - sim_matrix
    distance_matrix = np.maximum(distance_matrix, 0.0)
    np.fill_diagonal(distance_matrix, 0.0)
    distance_matrix = (distance_matrix + distance_matrix.T) / 2.0

    print(f"\n--- Step 4: Hierarchical clustering "
          f"({args.linkage}-linkage, threshold={args.identity_threshold}%) ---")
    n_clusters, _Z, scores, labels = find_clusters_threshold(
        distance_matrix, linkage_method=args.linkage, threshold=args.identity_threshold)
    print(f"  {n_clusters} clusters from {len(deduped_names)} unique variants")
    if scores and n_clusters in scores:
        print(f"  Silhouette score at k={n_clusters}: {scores[n_clusters]:.4f}")

    # Step 5: Report cluster membership
    print("\n--- Step 5: Clusters ---")
    from collections import defaultdict
    by_label = defaultdict(list)
    for i, cl in enumerate(labels):
        by_label[cl].append(deduped_names[i])
    for cl in sorted(by_label):
        members = by_label[cl]
        chr_ends = sorted({header_chr_end(m) for m in members})
        print(f"  Cluster {cl}: {', '.join(chr_ends)}")

    # Step 6: Write the clustered FASTA with rewritten headers (input to analyze_features.py)
    print("\n--- Step 6: Writing clustered FASTA ---")
    write_clustered_fasta_x(entries, deduped_names, labels, original_to_rep,
                            args.output_fasta)

    # Cluster assignments TSV (one row per original chr_end)
    assignments_path = os.path.join(args.output_dir, 'x_element_cluster_assignments.tsv')
    with open(assignments_path, 'w') as fh:
        fh.write('chr_end\tcluster\tlength\n')
        rep_to_label = {deduped_names[i]: cl for i, cl in enumerate(labels)}
        for header, seq in entries:
            rep = original_to_rep.get(header, header)
            cl = rep_to_label.get(rep, '?')
            fh.write(f'{header_chr_end(header)}\t{cl}\t{len(seq)}\n')
    print(f"  Assignments TSV: {assignments_path}")

    # Similarity matrix
    sim_path = os.path.join(args.output_dir, 'x_element_similarity_matrix.tsv')
    short = [header_chr_end(n) for n in deduped_names]
    pd.DataFrame(sim_matrix, index=short, columns=short).to_csv(sim_path, sep='\t')
    print(f"  Similarity matrix: {sim_path}")

    print("\nDone.")


if __name__ == '__main__':
    main()
