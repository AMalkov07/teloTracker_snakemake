#!/usr/bin/env python3
"""
Cluster Y prime sequences with configurable methods.

Two modes:
  --paper-defaults    Use hardcoded paper settings (Dudragne et al. 2026)
  (no flag)           Use the CUSTOM_* variables defined below

Usage:
    # Paper defaults (ignores CUSTOM_ variables)
    python scripts/cluster_yprimes_paper_method.py \
        --input-fasta references/extracted_yprimes_7302.fasta \
        --output-dir output_clustering/ \
        --paper-defaults

    # Custom settings (reads from CUSTOM_ variables below)
    python scripts/cluster_yprimes_paper_method.py \
        --input-fasta references/extracted_yprimes_7302.fasta \
        --output-dir output_clustering/

    # Override any setting via command line (takes priority over CUSTOM_ variables)
    python scripts/cluster_yprimes_paper_method.py \
        --input-fasta references/extracted_yprimes_7302.fasta \
        --output-dir output_clustering/ \
        --linkage complete --stop-mode threshold
"""

# =============================================================================
# CUSTOM SETTINGS — Edit these to set your preferred defaults.
# These are used when --paper-defaults is NOT specified.
# Command-line flags override these if provided.
# =============================================================================

CUSTOM_CONDENSE_HOMOPOLYMERS = True     # True = condense, False = skip
CUSTOM_DEDUP_THRESHOLD       = 99.9    # 99.9 = paper method, 100.0 = no dedup
CUSTOM_LINKAGE               = 'average'  # 'average', 'complete', or 'single'
CUSTOM_STOP_MODE             = 'threshold'  # 'silhouette' or 'threshold'
CUSTOM_IDENTITY_THRESHOLD    = 97.0    # Only used when CUSTOM_STOP_MODE = 'threshold'

# =============================================================================

import argparse
import os
import re
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score


# Colors for the heatmap figure
CLUSTER_COLORS = [
    '#7FC97F', '#1F78B4', '#B2ABD2', '#33A02C', '#8B4513',
    '#FF7F00', '#FFB6C1', '#6A3D9A', '#E31A1C', '#B15928',
    '#FDBF6F', '#CAB2D6', '#A6CEE3', '#FB9A99', '#00CED1',
    '#FFD700', '#808080', '#556B2F', '#DC143C', '#4169E1',
]

# Color names for FASTA header output (IDn_ColorName)
ID_COLOR_NAMES = [
    'Gray', 'Red', 'Green', 'Orange', 'Purple',
    'Blue', 'Yellow', 'Cyan', 'Magenta', 'Brown',
    'Pink', 'Teal', 'Olive', 'Navy', 'Coral',
    'Lavender', 'Maroon', 'Gold', 'Lime', 'Slate',
]


def parse_args():
    p = argparse.ArgumentParser(
        description='Cluster Y prime sequences (configurable method)')

    p.add_argument('input_fasta',
                   help='Input FASTA of extracted Y primes')
    p.add_argument('--output-dir', default='.',
                   help='Output directory for results and figures (default: current directory)')

    # Mode selection
    p.add_argument('--paper-defaults', action='store_true',
                   help='Use hardcoded paper settings (overrides CUSTOM_ variables and all flags below)')

    # Configurable steps (defaults come from CUSTOM_ variables above)
    p.add_argument('--no-condense-homopolymers', action='store_true', default=None,
                   help='Skip homopolymer condensation')
    p.add_argument('--condense-homopolymers', action='store_true', default=None,
                   help='Enable homopolymer condensation')
    p.add_argument('--dedup-threshold', type=float, default=None,
                   help='Similarity threshold for deduplication (100 = skip)')
    p.add_argument('--linkage', choices=['average', 'complete', 'single'], default=None,
                   help='Linkage method for hierarchical clustering')
    p.add_argument('--stop-mode', choices=['silhouette', 'threshold'], default=None,
                   help='How to determine cluster count')
    p.add_argument('--identity-threshold', type=float, default=None,
                   help='Identity threshold for --stop-mode threshold')

    p.add_argument('--output-fasta', default=None,
                   help='Output FASTA with rewritten IDn_Color headers (for pipeline integration)')
    p.add_argument('--label', default=None,
                   help='Label for output filenames (auto-generated if not set)')
    p.add_argument('--threads', type=int, default=4)

    args = p.parse_args()

    # Resolve settings: --paper-defaults > CLI flags > CUSTOM_ variables
    if args.paper_defaults:
        args.condense = True
        args.dedup_threshold = 99.9
        args.linkage = 'average'
        args.stop_mode = 'silhouette'
        args.identity_threshold = 97.0
    else:
        # Condense: check explicit flags, fall back to CUSTOM_
        if args.no_condense_homopolymers:
            args.condense = False
        elif args.condense_homopolymers:
            args.condense = True
        else:
            args.condense = CUSTOM_CONDENSE_HOMOPOLYMERS

        if args.dedup_threshold is None:
            args.dedup_threshold = CUSTOM_DEDUP_THRESHOLD
        if args.linkage is None:
            args.linkage = CUSTOM_LINKAGE
        if args.stop_mode is None:
            args.stop_mode = CUSTOM_STOP_MODE
        if args.identity_threshold is None:
            args.identity_threshold = CUSTOM_IDENTITY_THRESHOLD

    # Auto-generate label
    if args.label is None:
        if args.paper_defaults:
            args.label = 'paper_defaults'
        else:
            parts = []
            parts.append('condense' if args.condense else 'nocondense')
            parts.append(f'dedup{args.dedup_threshold:g}')
            parts.append(args.linkage)
            if args.stop_mode == 'silhouette':
                parts.append('silhouette')
            else:
                parts.append(f'thresh{args.identity_threshold:g}')
            args.label = '_'.join(parts)

    return args


# =============================================================================
# FASTA I/O
# =============================================================================

def load_fasta(fasta_path):
    entries = []
    current_header = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    entries.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)
    if current_header is not None:
        entries.append((current_header, ''.join(current_seq)))
    return entries


def write_fasta(entries, output_path):
    with open(output_path, 'w') as f:
        for header, seq in entries:
            f.write(f'>{header}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def rewrite_header(header, new_id, new_color):
    """Rewrite the IDn_Color portion of a Y prime header.

    Preserves the location, Size, and Type fields.
    Format: Y_Prime_locations#Size/Type/IDn_Color
    """
    if '#' not in header:
        return f'{header}#Unknown/Unknown/{new_id}_{new_color}'

    prefix, metadata = header.split('#', 1)
    parts = metadata.split('/')
    size = parts[0] if len(parts) >= 1 else 'Unknown'
    array_type = parts[1] if len(parts) >= 2 else 'Unknown'

    return f'{prefix}#{size}/{array_type}/{new_id}_{new_color}'


def assign_ids_and_colors(cluster_labels, seq_names):
    """Assign IDn and color name to each cluster.

    Clusters ordered by size (largest first), then alphabetically.
    Returns dict: seq_index -> (id_str, color_name)
    """
    from collections import defaultdict

    cluster_info = defaultdict(lambda: {'count': 0, 'first_name': None, 'members': []})
    for i, cl in enumerate(cluster_labels):
        cluster_info[cl]['count'] += 1
        cluster_info[cl]['members'].append(i)
        if cluster_info[cl]['first_name'] is None:
            cluster_info[cl]['first_name'] = seq_names[i]

    sorted_clusters = sorted(
        cluster_info.keys(),
        key=lambda c: (-cluster_info[c]['count'], cluster_info[c]['first_name'])
    )

    assignments = {}
    for new_id_num, old_label in enumerate(sorted_clusters, start=1):
        id_str = f'ID{new_id_num}'
        color = ID_COLOR_NAMES[(new_id_num - 1) % len(ID_COLOR_NAMES)]
        for idx in cluster_info[old_label]['members']:
            assignments[idx] = (id_str, color)

    return assignments


def write_clustered_fasta(original_entries, deduped_names, cluster_labels,
                          original_to_rep, output_path):
    """Write a FASTA with rewritten IDn_Color headers based on cluster assignments.

    Uses the ORIGINAL (uncondensed) sequences, with headers rewritten to
    reflect the new cluster IDs. This ensures the output is compatible with
    the rest of the pipeline regardless of clustering settings.

    original_to_rep maps every original header to its deduped representative.
    """
    id_assignments = assign_ids_and_colors(cluster_labels, deduped_names)

    # Map deduped names to their ID assignments
    name_to_id = {}
    for i, name in enumerate(deduped_names):
        name_to_id[name] = id_assignments[i]

    output_entries = []
    for header, seq in original_entries:
        # Find the deduped representative for this entry
        rep_name = original_to_rep.get(header, header)
        id_color = name_to_id.get(rep_name)

        if id_color is None:
            new_id, new_color = 'Unknown', 'Gray'
        else:
            new_id, new_color = id_color

        new_header = rewrite_header(header, new_id, new_color)
        output_entries.append((new_header, seq))

    write_fasta(output_entries, output_path)
    print(f"  Clustered FASTA: {output_path}")


# =============================================================================
# Step 1: Homopolymer condensation
# =============================================================================

def condense_homopolymers(seq, max_len=4):
    """Condense homopolymers longer than max_len to max_len."""
    if not seq:
        return seq
    result = [seq[0]]
    run_len = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            run_len += 1
            if run_len <= max_len:
                result.append(seq[i])
        else:
            run_len = 1
            result.append(seq[i])
    return ''.join(result)


# =============================================================================
# BLAST and similarity
# =============================================================================

def run_blast(fasta_path, threads=4):
    """Run all-vs-all BLAST. Returns DataFrame."""
    with tempfile.TemporaryDirectory(prefix='yprime_cluster_') as tmp_dir:
        db_path = os.path.join(tmp_dir, 'db')
        blast_out = os.path.join(tmp_dir, 'blast.txt')

        subprocess.run(
            ['makeblastdb', '-in', fasta_path, '-dbtype', 'nucl', '-out', db_path],
            check=True, capture_output=True)

        subprocess.run(
            ['blastn', '-query', fasta_path, '-db', db_path,
             '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend sstart send',
             '-evalue', '1e-10', '-max_target_seqs', '100',
             '-num_threads', str(threads), '-out', blast_out],
            check=True, capture_output=True)

        columns = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
                    'qstart', 'qend', 'sstart', 'send']
        df = pd.read_csv(blast_out, sep='\t', names=columns, header=None)

    return df


def calc_similarity_matrix(blast_df, seq_names, seq_lengths):
    """Compute similarity using the paper's formula:
    sim = %identity * alignment_length / alignment_span
    where alignment_span = alignment_length + max(unaligned_query, unaligned_subject)
    """
    n = len(seq_names)
    # Map both full header and first-word to index (BLAST uses first word only)
    name_to_idx = {}
    for i, name in enumerate(seq_names):
        name_to_idx[name] = i
        first_word = name.split()[0]
        if first_word not in name_to_idx:
            name_to_idx[first_word] = i
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
        pident = best['pident']
        aln_len = best['length']
        qlen = best['qlen']
        slen = best['slen']

        q_aligned = abs(best['qend'] - best['qstart']) + 1
        s_aligned = abs(best['send'] - best['sstart']) + 1
        q_unaligned = qlen - q_aligned
        s_unaligned = slen - s_aligned

        alignment_span = aln_len + max(q_unaligned, s_unaligned)
        similarity = (pident / 100.0) * aln_len / alignment_span * 100.0

        sim_matrix[qi, si] = max(sim_matrix[qi, si], similarity)
        sim_matrix[si, qi] = max(sim_matrix[si, qi], similarity)

    return sim_matrix


# =============================================================================
# Deduplication
# =============================================================================

def deduplicate_yprimes(entries, sim_matrix, seq_names, threshold=99.9):
    """Merge Y primes with similarity > threshold."""
    n = len(entries)
    merged_into = list(range(n))

    for i in range(n):
        for j in range(i+1, n):
            if sim_matrix[i, j] >= threshold:
                rep_i = i
                while merged_into[rep_i] != rep_i:
                    rep_i = merged_into[rep_i]
                rep_j = j
                while merged_into[rep_j] != rep_j:
                    rep_j = merged_into[rep_j]
                if rep_i != rep_j:
                    merged_into[rep_j] = rep_i

    for i in range(n):
        rep = i
        while merged_into[rep] != rep:
            rep = merged_into[rep]
        merged_into[i] = rep

    groups = {}
    for i in range(n):
        rep = merged_into[i]
        if rep not in groups:
            groups[rep] = []
        groups[rep].append(i)

    merged_entries = []
    merged_names = []
    original_to_rep = {}  # maps every original seq_name to its representative deduped name

    for new_idx, (rep, members) in enumerate(sorted(groups.items())):
        header = entries[rep][0]
        seq = entries[rep][1]
        rep_name = seq_names[rep]
        if len(members) > 1:
            member_names = [seq_names[m] for m in members]
            print(f"  Merged {len(members)} variants: {', '.join(m[:30] for m in member_names)}")
        merged_entries.append((header, seq))
        merged_names.append(rep_name)
        for m in members:
            original_to_rep[seq_names[m]] = rep_name

    return merged_entries, merged_names, original_to_rep


# =============================================================================
# Clustering
# =============================================================================

def find_clusters_silhouette(distance_matrix, linkage_method='average', max_clusters=None):
    """Find optimal cluster count using silhouette score."""
    n = distance_matrix.shape[0]
    if max_clusters is None:
        max_clusters = min(n - 1, 20)

    condensed = squareform(distance_matrix)
    Z = linkage(condensed, method=linkage_method)

    scores = {}
    for k in range(2, max_clusters + 1):
        labels = fcluster(Z, k, criterion='maxclust')
        if len(set(labels)) < 2:
            continue
        score = silhouette_score(distance_matrix, labels, metric='precomputed')
        scores[k] = score

    if not scores:
        return 2, Z, scores, fcluster(Z, 2, criterion='maxclust')

    best_k = max(scores, key=scores.get)
    labels = fcluster(Z, best_k, criterion='maxclust')
    return best_k, Z, scores, labels


def find_clusters_threshold(distance_matrix, linkage_method='complete', threshold=97.0):
    """Find clusters using a fixed identity threshold."""
    condensed = squareform(distance_matrix)
    Z = linkage(condensed, method=linkage_method)
    labels = fcluster(Z, t=(100.0 - threshold), criterion='distance')
    n_clusters = len(set(labels))

    # Also compute silhouette scores for informational purposes
    n = distance_matrix.shape[0]
    max_k = min(n - 1, 20)
    scores = {}
    for k in range(2, max_k + 1):
        k_labels = fcluster(Z, k, criterion='maxclust')
        if len(set(k_labels)) < 2:
            continue
        scores[k] = silhouette_score(distance_matrix, k_labels, metric='precomputed')

    return n_clusters, Z, scores, labels


# =============================================================================
# Figure 3A: Similarity matrix heatmap with dendrogram
# =============================================================================

def create_figure_3a(sim_matrix, seq_names, cluster_labels, lengths, linkage_method,
                     output_path):
    """Create a Figure 3A-style clustermap."""
    n_clusters = len(set(cluster_labels))
    color_map = {}
    for i, cl in enumerate(sorted(set(cluster_labels))):
        color_map[cl] = CLUSTER_COLORS[i % len(CLUSTER_COLORS)]

    row_colors = [color_map[cluster_labels[i]] for i in range(len(seq_names))]

    short_names = []
    for name in seq_names:
        m = re.match(r'Y_Prime_(.*?)#', name)
        short_names.append(m.group(1) if m else name[:20])

    df_sim = pd.DataFrame(sim_matrix, index=short_names, columns=short_names)
    cmap = sns.color_palette("YlOrRd", as_cmap=True)

    g = sns.clustermap(
        df_sim,
        method=linkage_method,
        metric='euclidean',
        row_colors=row_colors,
        col_colors=row_colors,
        cmap=cmap,
        vmin=0, vmax=100,
        figsize=(14, 12),
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        xticklabels=True,
        yticklabels=True,
        linewidths=0.5,
    )

    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.tick_params(labelsize=7)

    # Length bar on the right
    reordered_idx = g.dendrogram_row.reordered_ind
    ax_len = g.fig.add_axes([0.92, 0.05, 0.06, 0.75])
    reordered_lengths = [lengths[i] for i in reordered_idx]
    ax_len.barh(range(len(reordered_lengths)), reordered_lengths,
                color=[row_colors[i] for i in reordered_idx],
                edgecolor='black', linewidth=0.3, height=0.8)
    ax_len.set_ylim(-0.5, len(reordered_lengths) - 0.5)
    ax_len.set_xlabel('Length (bp)', fontsize=8)
    ax_len.tick_params(labelsize=6)
    ax_len.set_yticks([])
    ax_len.invert_yaxis()

    g.fig.suptitle(f'Hierarchical Clustering of Y\' Elements ({n_clusters} clusters)',
                   fontsize=14, fontweight='bold', y=1.02)

    # Cluster legend
    legend_handles = []
    for cl in sorted(set(cluster_labels)):
        count = sum(1 for l in cluster_labels if l == cl)
        legend_handles.append(plt.Line2D([0], [0], marker='s', color='w',
                              markerfacecolor=color_map[cl], markersize=10,
                              label=f'Cluster {cl} ({count})'))
    g.ax_heatmap.legend(handles=legend_handles, loc='upper left',
                        bbox_to_anchor=(1.15, 1.0), fontsize=7,
                        title='Clusters', title_fontsize=8)

    g.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"  Figure 3A saved to: {output_path}")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    label = args.label

    print("=" * 80)
    print("Y Prime Clustering")
    print("=" * 80)
    print(f"  Mode:                     {'Paper defaults' if args.paper_defaults else 'Custom'}")
    print(f"  Homopolymer condensation: {'ON' if args.condense else 'OFF'}")
    print(f"  Deduplication threshold:  {args.dedup_threshold}%")
    print(f"  Linkage method:           {args.linkage}")
    print(f"  Stop mode:                {args.stop_mode}", end='')
    if args.stop_mode == 'threshold':
        print(f" ({args.identity_threshold}%)")
    else:
        print()
    print(f"  Output label:             {label}")

    # Load sequences
    entries = load_fasta(args.input_fasta)
    print(f"\nLoaded {len(entries)} Y prime sequences")

    # Step 1: Homopolymer condensation
    if args.condense:
        print("\n--- Step 1: Condense homopolymers >4bp ---")
        condensed_entries = []
        n_condensed = 0
        for header, seq in entries:
            condensed = condense_homopolymers(seq.upper())
            condensed_entries.append((header, condensed))
            if len(seq) != len(condensed):
                n_condensed += 1
                name = header.split('#')[0] if '#' in header else header[:30]
                print(f"  {name}: {len(seq)} -> {len(condensed)} bp "
                      f"({len(seq) - len(condensed)} bp condensed)")
        print(f"  {n_condensed}/{len(entries)} sequences had homopolymers condensed")
    else:
        print("\n--- Step 1: Homopolymer condensation SKIPPED ---")
        condensed_entries = [(h, s.upper()) for h, s in entries]

    # Write working FASTA
    working_fasta = os.path.join(args.output_dir, f'{label}_working_yprimes.fasta')
    write_fasta(condensed_entries, working_fasta)

    seq_names = [h for h, s in condensed_entries]
    seq_lengths = {h: len(s) for h, s in condensed_entries}

    # Step 2: All-vs-all BLAST
    print("\n--- Step 2: All-vs-all BLAST ---")
    blast_df = run_blast(working_fasta, threads=args.threads)
    print(f"  BLAST produced {len(blast_df)} hits")

    sim_matrix = calc_similarity_matrix(blast_df, seq_names, seq_lengths)

    # Print similarity range
    n = len(seq_names)
    off_diag = sim_matrix[np.triu_indices(n, k=1)]
    if len(off_diag) > 0:
        print(f"  Similarity range: {off_diag.min():.1f}% - {off_diag.max():.1f}% "
              f"(mean {off_diag.mean():.1f}%)")

    # Step 3: Deduplication
    original_to_rep = {name: name for name in seq_names}  # identity mapping by default

    if args.dedup_threshold < 100.0:
        print(f"\n--- Step 3: Deduplicate at {args.dedup_threshold}% similarity ---")
        n_before = len(condensed_entries)
        deduped_entries, deduped_names, original_to_rep = deduplicate_yprimes(
            condensed_entries, sim_matrix, seq_names, threshold=args.dedup_threshold)
        n_after = len(deduped_entries)
        print(f"  {n_before} sequences -> {n_after} unique variants")

        if n_after < n_before:
            deduped_fasta = os.path.join(args.output_dir, f'{label}_deduped_yprimes.fasta')
            write_fasta(deduped_entries, deduped_fasta)
            blast_df2 = run_blast(deduped_fasta, threads=args.threads)
            deduped_lengths = {h: len(s) for h, s in deduped_entries}
            sim_matrix = calc_similarity_matrix(blast_df2, deduped_names, deduped_lengths)
        else:
            deduped_names = seq_names
    else:
        print(f"\n--- Step 3: Deduplication SKIPPED (threshold=100%) ---")
        deduped_entries = condensed_entries
        deduped_names = seq_names

    # Step 4: Hierarchical clustering
    distance_matrix = 100.0 - sim_matrix
    distance_matrix = np.maximum(distance_matrix, 0.0)
    np.fill_diagonal(distance_matrix, 0.0)
    distance_matrix = (distance_matrix + distance_matrix.T) / 2.0

    if args.stop_mode == 'silhouette':
        print(f"\n--- Step 4: Hierarchical clustering ({args.linkage}-linkage + silhouette score) ---")
        best_k, Z, scores, labels = find_clusters_silhouette(
            distance_matrix, linkage_method=args.linkage)
        print(f"\n  Silhouette scores by cluster count:")
        for k in sorted(scores.keys()):
            marker = " <-- best" if k == best_k else ""
            print(f"    k={k:>2}: {scores[k]:.4f}{marker}")
        print(f"\n  Optimal number of clusters: {best_k}")
    else:
        print(f"\n--- Step 4: Hierarchical clustering ({args.linkage}-linkage, "
              f"threshold={args.identity_threshold}%) ---")
        best_k, Z, scores, labels = find_clusters_threshold(
            distance_matrix, linkage_method=args.linkage, threshold=args.identity_threshold)
        print(f"\n  Number of clusters at {args.identity_threshold}% threshold: {best_k}")
        if scores:
            sil_at_k = scores.get(best_k, None)
            if sil_at_k is not None:
                print(f"  Silhouette score at k={best_k}: {sil_at_k:.4f}")

    # Step 5: Output results
    print(f"\n--- Step 5: Results ---")
    print(f"\n  {best_k} clusters from {len(deduped_names)} unique variants:")

    cluster_members = {}
    for i, cl in enumerate(labels):
        if cl not in cluster_members:
            cluster_members[cl] = []
        cluster_members[cl].append(deduped_names[i])

    for cl in sorted(cluster_members.keys()):
        members = cluster_members[cl]
        member_lengths = []
        for name in members:
            for h, s in entries:
                if h == name:
                    member_lengths.append(len(s))
                    break
        size_class = "Long" if member_lengths and np.mean(member_lengths) > 6400 else "Short"
        print(f"\n  Cluster {cl} ({len(members)} members, {size_class}):")
        for j, name in enumerate(members):
            short = name.split('#')[0] if '#' in name else name[:40]
            length = member_lengths[j] if j < len(member_lengths) else '?'
            print(f"    {short}  ({length} bp)")

    # Write cluster assignments TSV
    assignments_path = os.path.join(args.output_dir, f'{label}_cluster_assignments.tsv')
    with open(assignments_path, 'w') as f:
        f.write('sequence_name\tcluster\toriginal_length\n')
        for i, name in enumerate(deduped_names):
            orig_len = 0
            for h, s in entries:
                if h == name:
                    orig_len = len(s)
                    break
            f.write(f'{name}\t{labels[i]}\t{orig_len}\n')
    print(f"\n  Cluster assignments: {assignments_path}")

    # Write similarity matrix
    sim_path = os.path.join(args.output_dir, f'{label}_similarity_matrix.tsv')
    short_names = []
    for name in deduped_names:
        m = re.match(r'Y_Prime_(.*?)#', name)
        short_names.append(m.group(1) if m else name[:20])
    pd.DataFrame(sim_matrix, index=short_names, columns=short_names).to_csv(sim_path, sep='\t')
    print(f"  Similarity matrix: {sim_path}")

    # Write silhouette scores
    if scores:
        scores_path = os.path.join(args.output_dir, f'{label}_silhouette_scores.tsv')
        with open(scores_path, 'w') as f:
            f.write('n_clusters\tsilhouette_score\n')
            for k in sorted(scores.keys()):
                f.write(f'{k}\t{scores[k]:.6f}\n')
        print(f"  Silhouette scores: {scores_path}")

    # Write clustered FASTA (with original uncondensed sequences, rewritten headers)
    if args.output_fasta:
        print("\n--- Writing clustered FASTA ---")
        write_clustered_fasta(entries, deduped_names, labels, original_to_rep,
                              args.output_fasta)

    # Create Figure 3A
    print("\n--- Creating figures ---")
    orig_lengths = []
    for name in deduped_names:
        for h, s in entries:
            if h == name:
                orig_lengths.append(len(s))
                break

    fig3a_path = os.path.join(args.output_dir, f'{label}_clustermap.png')
    create_figure_3a(sim_matrix, deduped_names, labels, orig_lengths,
                     args.linkage, fig3a_path)

    # Silhouette score plot
    if scores:
        fig_sil, ax_sil = plt.subplots(figsize=(8, 4))
        ks = sorted(scores.keys())
        sil_vals = [scores[k] for k in ks]
        ax_sil.plot(ks, sil_vals, 'o-', color='#1F78B4', linewidth=2)
        ax_sil.axvline(x=best_k, color='red', linestyle='--', alpha=0.7,
                       label=f'Selected k={best_k}')
        ax_sil.set_xlabel('Number of clusters', fontsize=11)
        ax_sil.set_ylabel('Silhouette score', fontsize=11)
        ax_sil.set_title('Silhouette Score vs Number of Clusters',
                         fontsize=13, fontweight='bold')
        ax_sil.legend(fontsize=10)
        ax_sil.set_xticks(ks)
        ax_sil.grid(True, alpha=0.3)
        sil_path = os.path.join(args.output_dir, f'{label}_silhouette_plot.png')
        fig_sil.savefig(sil_path, dpi=150, bbox_inches='tight')
        plt.close(fig_sil)
        print(f"  Silhouette plot: {sil_path}")

    print("\nDone.")


if __name__ == '__main__':
    main()
