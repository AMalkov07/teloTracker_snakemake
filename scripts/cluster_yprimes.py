#!/usr/bin/env python3
"""
Cluster Y prime sequences by pairwise global identity and assign group IDs.

Takes extracted Y prime sequences (from extract_yprime_fasta.py), computes
all-vs-all BLAST, builds a global identity matrix, performs complete-linkage
hierarchical clustering, and rewrites FASTA headers with computed IDn_Color
group assignments.

Optionally produces a diagnostic variability report via MAFFT multiple
sequence alignment, showing exactly which regions differentiate the groups.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from math import log2

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


# Color palette for group assignment
COLOR_PALETTE = [
    'Gray', 'Red', 'Green', 'Orange', 'Purple',
    'Blue', 'Yellow', 'Cyan', 'Magenta', 'Brown',
    'Pink', 'Teal', 'Olive', 'Navy', 'Coral',
    'Lavender', 'Maroon', 'Gold', 'Lime', 'Slate',
]


def parse_args():
    p = argparse.ArgumentParser(
        description='Cluster Y prime sequences by pairwise identity and assign group IDs'
    )
    p.add_argument('--input-fasta', required=True,
                   help='Input FASTA of extracted Y primes')
    p.add_argument('--output-fasta', required=True,
                   help='Output FASTA with updated IDn_Color headers')
    p.add_argument('--output-tsv', required=True,
                   help='Output TSV mapping file')
    p.add_argument('--identity-threshold', type=float, default=97.0,
                   help='Pairwise identity threshold for grouping (default: 97.0)')
    p.add_argument('--threads', type=int, default=4,
                   help='Threads for BLAST/MAFFT')
    p.add_argument('--diagnostic-report', default=None,
                   help='Output directory for MAFFT-based diagnostic variability report')
    return p.parse_args()


# =============================================================================
# FASTA I/O
# =============================================================================

def load_fasta(fasta_path):
    """Load FASTA file. Returns list of (header, sequence) tuples."""
    entries = []
    current_header = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    entries.append((current_header, ''.join(current_seq)))
                current_header = line[1:]  # strip '>'
                current_seq = []
            elif line:
                current_seq.append(line)
    if current_header is not None:
        entries.append((current_header, ''.join(current_seq)))
    return entries


def write_fasta(entries, output_path):
    """Write list of (header, sequence) tuples to FASTA."""
    with open(output_path, 'w') as f:
        for header, seq in entries:
            f.write(f'>{header}\n')
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def parse_header(header):
    """Parse Y prime FASTA header.

    Format: Y_Prime_chr2L1#Short/Solo/ID4_Green-Light
    Returns dict with keys: locations, size, array_type, id, color, prefix
    """
    result = {'raw': header}

    if '#' not in header:
        result['prefix'] = header
        result['size'] = 'Unknown'
        result['array_type'] = 'Unknown'
        result['id'] = 'Unknown'
        result['color'] = 'Unknown'
        return result

    prefix, metadata = header.split('#', 1)
    result['prefix'] = prefix

    parts = metadata.split('/')
    result['size'] = parts[0] if len(parts) >= 1 else 'Unknown'
    result['array_type'] = parts[1] if len(parts) >= 2 else 'Unknown'

    if len(parts) >= 3:
        id_color = parts[2]
        if '_' in id_color:
            underscore_idx = id_color.index('_')
            result['id'] = id_color[:underscore_idx]
            result['color'] = id_color[underscore_idx + 1:]
        else:
            result['id'] = id_color
            result['color'] = 'Unknown'
    else:
        result['id'] = 'Unknown'
        result['color'] = 'Unknown'

    return result


def rewrite_header(header, new_id, new_color):
    """Rewrite the IDn_Color portion of a Y prime header."""
    if '#' not in header:
        return f'{header}#Unknown/Unknown/{new_id}_{new_color}'

    prefix, metadata = header.split('#', 1)
    parts = metadata.split('/')

    size = parts[0] if len(parts) >= 1 else 'Unknown'
    array_type = parts[1] if len(parts) >= 2 else 'Unknown'

    return f'{prefix}#{size}/{array_type}/{new_id}_{new_color}'


# =============================================================================
# All-vs-all BLAST
# =============================================================================

def run_allvsall_blast(fasta_path, threads=4):
    """Run BLASTn all-vs-all. Returns DataFrame."""
    with tempfile.TemporaryDirectory(prefix='yprime_cluster_') as tmp_dir:
        db_path = os.path.join(tmp_dir, 'yprime_db')
        blast_out = os.path.join(tmp_dir, 'blast_results.txt')

        # Build BLAST database
        subprocess.run(
            ['makeblastdb', '-in', fasta_path, '-dbtype', 'nucl', '-out', db_path],
            check=True, capture_output=True
        )

        # Run BLAST
        subprocess.run(
            ['blastn', '-query', fasta_path, '-db', db_path,
             '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend sstart send nident',
             '-evalue', '1e-10', '-max_target_seqs', '100',
             '-num_threads', str(threads), '-out', blast_out],
            check=True, capture_output=True
        )

        # Parse results
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
                    'qstart', 'qend', 'sstart', 'send', 'nident']
        df = pd.read_csv(blast_out, sep='\t', names=columns, header=None)

    return df


# =============================================================================
# Identity matrix
# =============================================================================

def compute_identity_matrix(blast_df, seq_names, seq_lengths):
    """Compute global pairwise identity from BLAST alignments.

    Global identity = nident / max(qlen, slen) * 100
    This penalizes partial alignments — unaligned regions count as mismatches.

    For pairs with multiple HSPs, takes the single best (longest) HSP to avoid
    double-counting. Self-hits are set to 100%.
    """
    n = len(seq_names)
    name_to_idx = {name: i for i, name in enumerate(seq_names)}
    identity_matrix = np.zeros((n, n))

    # Fill diagonal
    np.fill_diagonal(identity_matrix, 100.0)

    if blast_df.empty:
        return identity_matrix

    # For each pair, take the best HSP (highest nident)
    for (qseqid, sseqid), group in blast_df.groupby(['qseqid', 'sseqid']):
        if qseqid == sseqid:
            continue

        qi = name_to_idx.get(qseqid)
        si = name_to_idx.get(sseqid)
        if qi is None or si is None:
            continue

        # Take best HSP by nident
        best = group.loc[group['nident'].idxmax()]
        max_len = max(best['qlen'], best['slen'])
        global_ident = (best['nident'] / max_len) * 100.0

        # Keep the maximum across both directions (q->s and s->q)
        identity_matrix[qi, si] = max(identity_matrix[qi, si], global_ident)
        identity_matrix[si, qi] = max(identity_matrix[si, qi], global_ident)

    return identity_matrix


# =============================================================================
# Clustering
# =============================================================================

def cluster_sequences(identity_matrix, threshold=97.0):
    """Complete-linkage hierarchical clustering on identity matrix.

    Returns array of integer cluster labels (1-indexed from fcluster).
    """
    n = identity_matrix.shape[0]

    if n == 1:
        return np.array([1])

    # Convert to distance matrix
    distance_matrix = 100.0 - identity_matrix
    # Ensure non-negative and symmetric
    distance_matrix = np.maximum(distance_matrix, 0.0)
    np.fill_diagonal(distance_matrix, 0.0)
    distance_matrix = (distance_matrix + distance_matrix.T) / 2.0

    # Condensed form for scipy
    condensed = squareform(distance_matrix)

    # Complete linkage: all pairs within cluster must be within threshold
    Z = linkage(condensed, method='complete')

    # Cut at threshold distance
    labels = fcluster(Z, t=(100.0 - threshold), criterion='distance')

    return labels


def assign_ids_and_colors(cluster_labels, seq_names):
    """Assign IDn and Color to each cluster.

    Clusters ordered by size (largest first), then alphabetically.
    Returns dict: seq_index -> (id_str, color_str)
    """
    # Count cluster sizes and find first member name
    cluster_info = defaultdict(lambda: {'count': 0, 'first_name': None, 'members': []})
    for i, label in enumerate(cluster_labels):
        cluster_info[label]['count'] += 1
        cluster_info[label]['members'].append(i)
        if cluster_info[label]['first_name'] is None:
            cluster_info[label]['first_name'] = seq_names[i]

    # Sort clusters: largest first, then alphabetically by first member
    sorted_clusters = sorted(
        cluster_info.keys(),
        key=lambda c: (-cluster_info[c]['count'], cluster_info[c]['first_name'])
    )

    assignments = {}
    for new_id_num, old_label in enumerate(sorted_clusters, start=1):
        id_str = f'ID{new_id_num}'
        color = COLOR_PALETTE[(new_id_num - 1) % len(COLOR_PALETTE)]
        for idx in cluster_info[old_label]['members']:
            assignments[idx] = (id_str, color)

    return assignments


# =============================================================================
# Diagnostic variability report (MAFFT MSA)
# =============================================================================

def run_mafft(fasta_path, output_path, threads=4):
    """Run MAFFT multiple sequence alignment."""
    result = subprocess.run(
        ['mafft', '--auto', '--thread', str(threads), fasta_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"MAFFT error: {result.stderr}", file=sys.stderr)
        return False

    with open(output_path, 'w') as f:
        f.write(result.stdout)
    return True


def compute_variability_profile(alignment_entries, cluster_assignments, seq_names):
    """Compute per-position variability from MSA.

    Returns list of dicts with keys:
        position, entropy, consensus_base, base_counts, groups_at_position
    """
    if not alignment_entries:
        return []

    # Map aligned sequence names to cluster IDs
    name_to_idx = {name: i for i, name in enumerate(seq_names)}

    aln_length = len(alignment_entries[0][1])
    profile = []

    for pos in range(aln_length):
        base_counts = defaultdict(int)
        group_bases = defaultdict(lambda: defaultdict(int))
        total_non_gap = 0

        for header, seq in alignment_entries:
            # Match aligned header to original name
            aln_name = header.split('#')[0] if '#' in header else header
            aln_name_full = header

            idx = name_to_idx.get(aln_name_full)
            if idx is None:
                # Try matching by prefix
                for orig_name, orig_idx in name_to_idx.items():
                    if orig_name.startswith(aln_name) or aln_name.startswith(orig_name):
                        idx = orig_idx
                        break

            base = seq[pos].upper()
            base_counts[base] += 1

            if idx is not None and base != '-':
                group_id = cluster_assignments.get(idx, ('Unknown', 'Unknown'))[0]
                group_bases[group_id][base] += 1

            if base != '-':
                total_non_gap += 1

        # Shannon entropy
        entropy = 0.0
        if total_non_gap > 0:
            for base, count in base_counts.items():
                if base == '-':
                    continue
                p = count / total_non_gap
                if p > 0:
                    entropy -= p * log2(p)

        # Consensus base (most frequent non-gap)
        non_gap_counts = {b: c for b, c in base_counts.items() if b != '-'}
        consensus = max(non_gap_counts, key=non_gap_counts.get) if non_gap_counts else '-'

        # Fraction differing from consensus (among non-gap)
        consensus_count = non_gap_counts.get(consensus, 0)
        frac_differ = 1.0 - (consensus_count / total_non_gap) if total_non_gap > 0 else 0.0

        profile.append({
            'position': pos + 1,  # 1-based
            'entropy': entropy,
            'frac_differ': frac_differ,
            'consensus_base': consensus,
            'total_non_gap': total_non_gap,
            'total_seqs': len(alignment_entries),
            'gap_count': base_counts.get('-', 0),
            'base_counts': dict(base_counts),
            'group_bases': {k: dict(v) for k, v in group_bases.items()},
        })

    return profile


def identify_diagnostic_windows(profile, assignments, min_entropy=0.1, merge_gap=10):
    """Identify contiguous variable regions and annotate which groups they separate.

    Returns list of dicts: start, end, length, avg_entropy, max_entropy, separates_groups
    """
    # Find variable positions
    variable_positions = [p for p in profile if p['entropy'] > min_entropy]

    if not variable_positions:
        return []

    # Merge adjacent positions into windows
    windows = []
    current_start = variable_positions[0]['position']
    current_end = variable_positions[0]['position']
    current_entries = [variable_positions[0]]

    for p in variable_positions[1:]:
        if p['position'] - current_end <= merge_gap:
            current_end = p['position']
            current_entries.append(p)
        else:
            windows.append({
                'start': current_start,
                'end': current_end,
                'entries': current_entries
            })
            current_start = p['position']
            current_end = p['position']
            current_entries = [p]

    windows.append({
        'start': current_start,
        'end': current_end,
        'entries': current_entries
    })

    # Annotate each window
    all_groups = sorted(set(v[0] for v in assignments.values()))
    result = []

    for w in windows:
        entropies = [e['entropy'] for e in w['entries']]

        # Determine which groups differ in this window
        # For each variable position, check which groups have different consensus bases
        groups_that_differ = set()
        for entry in w['entries']:
            group_bases = entry.get('group_bases', {})
            # Find the majority base for each group at this position
            group_consensus = {}
            for gid, bases in group_bases.items():
                if bases:
                    group_consensus[gid] = max(bases, key=bases.get)

            # Find groups with different consensus
            consensus_values = list(group_consensus.values())
            if len(set(consensus_values)) > 1:
                # Multiple different bases across groups
                for gid, base in group_consensus.items():
                    groups_that_differ.add(gid)

        # Build a description of which groups separate
        if len(groups_that_differ) >= 2:
            separates = sorted(groups_that_differ)
        else:
            separates = all_groups  # Gap-based differentiation

        result.append({
            'start': w['start'],
            'end': w['end'],
            'length': w['end'] - w['start'] + 1,
            'avg_entropy': np.mean(entropies),
            'max_entropy': max(entropies),
            'num_variable_positions': len(w['entries']),
            'separates_groups': separates,
        })

    return result


def generate_size_class_annotation(entries, seq_names, cluster_assignments, output_dir, threads=4):
    """Generate annotated consensus for Long and Short Y prime size classes.

    For each size class, aligns all members via MAFFT, then annotates each
    position as conserved (identical across all members) or variable
    (with details on which groups differ and what bases they have).

    Outputs per size class:
        - <size>_yprime_alignment.fasta — the MAFFT alignment
        - <size>_yprime_annotated_regions.tsv — conserved/variable regions
        - <size>_yprime_consensus.fasta — consensus with variable sites as N
        - stdout summary
    """
    from collections import Counter

    # Split entries by size class
    size_classes = defaultdict(list)  # size -> [(header, seq, idx)]
    for i, (header, seq) in enumerate(entries):
        parsed = parse_header(header)
        size = parsed.get('size', 'Unknown')
        if size in ('Long', 'Short'):
            size_classes[size].append((header, seq, i))
        else:
            size_classes['Unknown'].append((header, seq, i))

    for size_class, members in sorted(size_classes.items()):
        if len(members) < 2:
            print(f"\n--- {size_class} Y primes: only {len(members)} sequence, skipping annotation ---")
            continue

        print(f"\n--- {size_class} Y prime annotation ({len(members)} sequences) ---")

        # Write temp FASTA for MAFFT
        size_fasta = os.path.join(output_dir, f'{size_class.lower()}_yprime_input.fasta')
        write_fasta([(h, s) for h, s, _ in members], size_fasta)

        # Run MAFFT
        aln_path = os.path.join(output_dir, f'{size_class.lower()}_yprime_alignment.fasta')
        if not run_mafft(size_fasta, aln_path, threads):
            print(f"  MAFFT failed for {size_class}, skipping")
            continue

        aln_entries = load_fasta(aln_path)
        if not aln_entries:
            continue

        aln_length = len(aln_entries[0][1])

        # Map aligned names to cluster IDs
        member_idx_map = {}  # aln position -> original idx
        for header, seq, idx in members:
            member_idx_map[header] = idx

        # Analyze each position
        positions = []
        consensus_seq = []
        for pos in range(aln_length):
            bases_by_group = defaultdict(list)
            all_bases = []

            for aln_header, aln_seq in aln_entries:
                base = aln_seq[pos].upper()
                # Find the original index for this aligned entry
                orig_idx = member_idx_map.get(aln_header)
                if orig_idx is None:
                    for h, s, idx in members:
                        if h.startswith(aln_header) or aln_header.startswith(h):
                            orig_idx = idx
                            break

                all_bases.append(base)
                if orig_idx is not None:
                    gid = cluster_assignments.get(orig_idx, ('?', '?'))[0]
                    bases_by_group[gid].append(base)

            non_gap = [b for b in all_bases if b != '-']

            if not non_gap:
                positions.append({
                    'pos': pos + 1, 'status': 'gap', 'consensus': '-',
                    'detail': '', 'groups': {}
                })
                consensus_seq.append('-')
                continue

            # Check if all non-gap bases are identical
            unique_bases = set(non_gap)
            consensus_base = Counter(non_gap).most_common(1)[0][0]
            num_gaps = all_bases.count('-')

            if len(unique_bases) == 1 and num_gaps == 0:
                status = 'conserved'
                detail = ''
            elif len(unique_bases) == 1 and num_gaps > 0:
                status = 'conserved_with_gaps'
                detail = f'{num_gaps}/{len(all_bases)} gaps'
            else:
                status = 'variable'
                # Build group detail: which group has which base
                group_summary = {}
                for gid, gbases in sorted(bases_by_group.items()):
                    non_gap_gbases = [b for b in gbases if b != '-']
                    if non_gap_gbases:
                        group_summary[gid] = Counter(non_gap_gbases).most_common(1)[0][0]
                    else:
                        group_summary[gid] = '-'
                detail = ' '.join(f'{gid}={b}' for gid, b in sorted(group_summary.items()))

            positions.append({
                'pos': pos + 1, 'status': status, 'consensus': consensus_base,
                'detail': detail, 'groups': {gid: Counter(non_gap_gbases).most_common(1)[0][0]
                                             for gid, gbases in bases_by_group.items()
                                             for non_gap_gbases in [[b for b in gbases if b != '-']]
                                             if non_gap_gbases}
            })
            consensus_seq.append(consensus_base if status == 'conserved' else 'N')

        # Write consensus FASTA (N at variable sites)
        consensus_path = os.path.join(output_dir, f'{size_class.lower()}_yprime_consensus.fasta')
        ungapped_consensus = ''.join(b for b in consensus_seq if b != '-')
        write_fasta([(f'{size_class}_Y_Prime_consensus', ungapped_consensus)], consensus_path)

        # Merge positions into regions
        regions = []
        current_status = positions[0]['status'] if positions else None
        region_start = 1
        region_details = []

        for p in positions:
            # Simplify: treat conserved_with_gaps as conserved for region merging
            p_status = 'conserved' if p['status'] in ('conserved', 'conserved_with_gaps', 'gap') else 'variable'
            if p_status != current_status:
                regions.append({
                    'start': region_start, 'end': p['pos'] - 1,
                    'status': current_status, 'details': region_details
                })
                region_start = p['pos']
                current_status = p_status
                region_details = []
            if p_status == 'variable':
                region_details.append(p)

        regions.append({
            'start': region_start, 'end': positions[-1]['pos'],
            'status': current_status, 'details': region_details
        })

        # Write annotated regions TSV
        regions_path = os.path.join(output_dir, f'{size_class.lower()}_yprime_annotated_regions.tsv')
        with open(regions_path, 'w') as f:
            f.write('start\tend\tlength\tstatus\tnum_variable_sites\tgroup_differences\n')
            for r in regions:
                length = r['end'] - r['start'] + 1
                num_var = len(r['details'])
                # Summarize group differences across variable sites in this region
                if r['details']:
                    # Collect all group->base mappings across variable positions
                    group_bases_across = defaultdict(list)
                    for d in r['details']:
                        for gid, base in d.get('groups', {}).items():
                            group_bases_across[gid].append(base)
                    # Create compact summary
                    group_summary_parts = []
                    for gid in sorted(group_bases_across.keys()):
                        bases = group_bases_across[gid]
                        summary = ''.join(bases[:10])
                        if len(bases) > 10:
                            summary += f'...({len(bases)})'
                        group_summary_parts.append(f'{gid}:{summary}')
                    group_diff = ' | '.join(group_summary_parts)
                else:
                    group_diff = ''
                f.write(f'{r["start"]}\t{r["end"]}\t{length}\t{r["status"]}\t{num_var}\t{group_diff}\n')

        print(f"  Alignment: {aln_path}")
        print(f"  Consensus: {consensus_path}")
        print(f"  Annotated regions: {regions_path}")

        # Print visual summary
        print(f"\n  {size_class} Y prime layout ({aln_length} alignment columns):")
        # Get group IDs present in this size class
        groups_in_class = sorted(set(
            cluster_assignments[idx][0] for _, _, idx in members
        ))
        print(f"  Groups: {', '.join(groups_in_class)}")

        for r in regions:
            length = r['end'] - r['start'] + 1
            if length < 5 and r['status'] == 'conserved':
                continue  # Skip tiny conserved spacers for readability
            num_var = len(r['details'])
            if r['status'] == 'variable':
                # Show which groups have which bases at this region
                group_bases_across = defaultdict(list)
                for d in r['details']:
                    for gid, base in d.get('groups', {}).items():
                        group_bases_across[gid].append(base)
                # Compact: group groups with same base pattern
                pattern_to_groups = defaultdict(list)
                for gid in sorted(group_bases_across.keys()):
                    pattern = ''.join(group_bases_across[gid])
                    pattern_to_groups[pattern].append(gid)

                diff_desc_parts = []
                for pattern, gids in sorted(pattern_to_groups.items(), key=lambda x: x[1][0]):
                    gids_str = ','.join(gids)
                    if len(pattern) <= 6:
                        diff_desc_parts.append(f"{gids_str}=[{pattern}]")
                    else:
                        diff_desc_parts.append(f"{gids_str}=[{pattern[:6]}...]")
                diff_desc = '  vs  '.join(diff_desc_parts)

                print(f"    {r['start']:>6}-{r['end']:<6}  ({length:>5} bp)  VARIABLE  "
                      f"{num_var} sites — {diff_desc}")
            else:
                print(f"    {r['start']:>6}-{r['end']:<6}  ({length:>5} bp)  conserved")

        # Clean up temp input file
        os.remove(size_fasta)


def generate_diagnostic_report(input_fasta, output_dir, seq_names, cluster_assignments, threads=4):
    """Generate MAFFT-based diagnostic variability report."""
    os.makedirs(output_dir, exist_ok=True)

    # Check MAFFT availability
    try:
        subprocess.run(['mafft', '--version'], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("WARNING: MAFFT not found, skipping diagnostic report", file=sys.stderr)
        return

    print("\n--- Diagnostic Variability Report ---")

    # Run MAFFT
    alignment_path = os.path.join(output_dir, 'alignment.fasta')
    print("Running MAFFT multiple sequence alignment...")
    if not run_mafft(input_fasta, alignment_path, threads):
        print("WARNING: MAFFT failed, skipping diagnostic report", file=sys.stderr)
        return

    # Load alignment
    alignment_entries = load_fasta(alignment_path)
    if not alignment_entries:
        print("WARNING: Empty alignment, skipping diagnostic report", file=sys.stderr)
        return

    print(f"Alignment: {len(alignment_entries)} sequences, {len(alignment_entries[0][1])} columns")

    # Compute per-position variability
    profile = compute_variability_profile(alignment_entries, cluster_assignments, seq_names)

    # Write variability profile TSV
    profile_path = os.path.join(output_dir, 'variability_profile.tsv')
    with open(profile_path, 'w') as f:
        f.write('position\tentropy\tfrac_differ\tconsensus\ttotal_non_gap\tgap_count\n')
        for p in profile:
            f.write(f"{p['position']}\t{p['entropy']:.4f}\t{p['frac_differ']:.4f}\t"
                    f"{p['consensus_base']}\t{p['total_non_gap']}\t{p['gap_count']}\n")
    print(f"Per-position variability: {profile_path}")

    # Identify diagnostic windows
    windows = identify_diagnostic_windows(profile, cluster_assignments)

    # Write diagnostic regions TSV
    regions_path = os.path.join(output_dir, 'diagnostic_regions.tsv')
    with open(regions_path, 'w') as f:
        f.write('start\tend\tlength\tavg_entropy\tmax_entropy\tnum_variable_positions\tseparates_groups\n')
        for w in windows:
            groups_str = ','.join(w['separates_groups'])
            f.write(f"{w['start']}\t{w['end']}\t{w['length']}\t"
                    f"{w['avg_entropy']:.4f}\t{w['max_entropy']:.4f}\t"
                    f"{w['num_variable_positions']}\t{groups_str}\n")
    print(f"Diagnostic regions: {regions_path}")

    # Print summary to stdout
    aln_length = len(alignment_entries[0][1])
    total_variable = sum(1 for p in profile if p['entropy'] > 0.1)
    total_conserved = sum(1 for p in profile if p['entropy'] <= 0.1 and p['total_non_gap'] > 0)

    print(f"\nAlignment length: {aln_length} columns")
    print(f"Conserved positions (entropy <= 0.1): {total_conserved}")
    print(f"Variable positions (entropy > 0.1): {total_variable}")
    print(f"Gap-only positions: {aln_length - total_variable - total_conserved}")

    if windows:
        print(f"\nDiagnostic regions ({len(windows)} found):")
        for w in windows:
            groups_str = ', '.join(w['separates_groups'])
            print(f"  Positions {w['start']}-{w['end']} "
                  f"({w['length']} bp, {w['num_variable_positions']} variable sites, "
                  f"max entropy {w['max_entropy']:.2f}) "
                  f"— separates: {groups_str}")
    else:
        print("\nNo diagnostic regions found (all positions are conserved or gap-only)")

    # Pairwise group comparison table
    all_groups = sorted(set(v[0] for v in cluster_assignments.values()))
    if len(all_groups) >= 2:
        # Build per-group consensus at each position
        group_consensus_seqs = {}
        for gid in all_groups:
            consensus = []
            for p in profile:
                gb = p.get('group_bases', {}).get(gid, {})
                if gb:
                    consensus.append(max(gb, key=gb.get))
                else:
                    consensus.append('-')
            group_consensus_seqs[gid] = consensus

        # Pairwise comparison
        from itertools import combinations
        comparison_path = os.path.join(output_dir, 'group_pairwise_comparison.tsv')
        with open(comparison_path, 'w') as f:
            f.write('group_A\tgroup_B\tlen_A\tlen_B\tshared_positions\tidentical_positions\t'
                    'different_positions\tglobal_pct_identity\tshared_pct_identity\tkey_diff_regions\n')

            for ga, gb in combinations(all_groups, 2):
                seq_a = group_consensus_seqs[ga]
                seq_b = group_consensus_seqs[gb]
                # Sequence lengths (non-gap bases)
                len_a = sum(1 for b in seq_a if b != '-')
                len_b = sum(1 for b in seq_b if b != '-')
                shared = 0
                identical = 0
                diff_positions = []
                for i, (a, b) in enumerate(zip(seq_a, seq_b)):
                    if a != '-' and b != '-':
                        shared += 1
                        if a == b:
                            identical += 1
                        else:
                            diff_positions.append(i + 1)  # 1-based

                # Global identity: identical / max(len_A, len_B)
                # Same metric used by the clustering algorithm
                max_len = max(len_a, len_b)
                global_pct = (identical / max_len * 100) if max_len > 0 else 0.0
                # Shared identity: identical / shared (for reference)
                shared_pct = (identical / shared * 100) if shared > 0 else 0.0

                # Summarize diff positions into ranges
                diff_ranges = []
                if diff_positions:
                    rng_start = diff_positions[0]
                    rng_end = diff_positions[0]
                    for dp in diff_positions[1:]:
                        if dp - rng_end <= 5:  # merge nearby diffs
                            rng_end = dp
                        else:
                            diff_ranges.append(f"{rng_start}-{rng_end}" if rng_start != rng_end else str(rng_start))
                            rng_start = dp
                            rng_end = dp
                    diff_ranges.append(f"{rng_start}-{rng_end}" if rng_start != rng_end else str(rng_start))

                key_regions = '; '.join(diff_ranges[:15])  # top 15 regions
                if len(diff_ranges) > 15:
                    key_regions += f' ... (+{len(diff_ranges) - 15} more)'

                f.write(f'{ga}\t{gb}\t{len_a}\t{len_b}\t{shared}\t{identical}\t{len(diff_positions)}\t'
                        f'{global_pct:.2f}\t{shared_pct:.2f}\t{key_regions}\n')

        print(f"Pairwise group comparison: {comparison_path}")

        # Print summary table to stdout
        print("\n--- Pairwise Group Comparison (global identity = identical / max(len_A, len_B)) ---")
        print(f"{'Group A':<8} {'Group B':<8} {'Len A':<8} {'Len B':<8} {'Identical':<10} {'Different':<10} {'Global %':<10} {'Shared %':<10}")
        for ga, gb in combinations(all_groups, 2):
            seq_a = group_consensus_seqs[ga]
            seq_b = group_consensus_seqs[gb]
            len_a = sum(1 for b in seq_a if b != '-')
            len_b = sum(1 for b in seq_b if b != '-')
            shared = sum(1 for a, b in zip(seq_a, seq_b) if a != '-' and b != '-')
            identical = sum(1 for a, b in zip(seq_a, seq_b) if a != '-' and b != '-' and a == b)
            different = shared - identical
            max_len = max(len_a, len_b)
            global_pct = (identical / max_len * 100) if max_len > 0 else 0.0
            shared_pct = (identical / shared * 100) if shared > 0 else 0.0
            print(f"{ga:<8} {gb:<8} {len_a:<8} {len_b:<8} {identical:<10} {different:<10} {global_pct:<10.2f} {shared_pct:<10.2f}")

    # Summary of conserved vs variable by region
    print("\n--- Region Summary ---")
    if profile:
        # Walk through profile and describe conserved/variable stretches
        in_variable = False
        region_start = 1
        regions_summary = []

        for p in profile:
            is_var = p['entropy'] > 0.1
            if is_var != in_variable:
                if region_start < p['position']:
                    region_type = "Variable" if in_variable else "Conserved"
                    regions_summary.append((region_start, p['position'] - 1, region_type))
                region_start = p['position']
                in_variable = is_var

        # Final region
        region_type = "Variable" if in_variable else "Conserved"
        regions_summary.append((region_start, profile[-1]['position'], region_type))

        # Merge small regions for readability (< 5bp variable regions into conserved)
        merged = []
        for start, end, rtype in regions_summary:
            length = end - start + 1
            if merged and rtype == merged[-1][2]:
                merged[-1] = (merged[-1][0], end, rtype)
            elif rtype == 'Variable' and length < 3 and merged and merged[-1][2] == 'Conserved':
                merged[-1] = (merged[-1][0], end, 'Conserved')
            else:
                merged.append((start, end, rtype))

        for start, end, rtype in merged:
            length = end - start + 1
            if length >= 10:  # Only print regions >= 10bp
                print(f"  {start:>6}-{end:<6}  ({length:>5} bp)  {rtype}")

    # Size-class annotation (Long vs Short)
    entries = load_fasta(input_fasta)
    generate_size_class_annotation(entries, seq_names, cluster_assignments, output_dir, threads)


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    print(f"Clustering Y prime sequences from: {args.input_fasta}")
    print(f"Identity threshold: {args.identity_threshold}%")

    # Step 1: Load sequences
    entries = load_fasta(args.input_fasta)
    if not entries:
        print("ERROR: No sequences found in input FASTA", file=sys.stderr)
        sys.exit(1)

    seq_names = [h for h, s in entries]
    seq_lengths = {h: len(s) for h, s in entries}
    print(f"Loaded {len(entries)} Y prime sequences")

    # Step 2: All-vs-all BLAST
    print("Running all-vs-all BLAST...")
    blast_df = run_allvsall_blast(args.input_fasta, threads=args.threads)
    print(f"BLAST produced {len(blast_df)} hits")

    # Step 3: Compute global identity matrix
    identity_matrix = compute_identity_matrix(blast_df, seq_names, seq_lengths)

    # Print identity matrix summary
    n = len(seq_names)
    off_diag = identity_matrix[np.triu_indices(n, k=1)]
    if len(off_diag) > 0:
        print(f"Pairwise identity range: {off_diag.min():.1f}% - {off_diag.max():.1f}% "
              f"(mean {off_diag.mean():.1f}%)")

    # Step 4: Hierarchical clustering
    labels = cluster_sequences(identity_matrix, threshold=args.identity_threshold)
    n_clusters = len(set(labels))
    print(f"Clustering produced {n_clusters} groups at {args.identity_threshold}% threshold")

    # Step 5: Assign IDs and colors
    assignments = assign_ids_and_colors(labels, seq_names)

    # Print cluster summary
    cluster_members = defaultdict(list)
    for idx, (id_str, color) in assignments.items():
        cluster_members[(id_str, color)].append(seq_names[idx])

    print("\nCluster assignments:")
    for (id_str, color), members in sorted(cluster_members.items()):
        print(f"  {id_str}_{color} ({len(members)} members):")
        for m in members:
            parsed = parse_header(m)
            print(f"    {parsed['prefix']}  [{parsed['size']}]")

    # Step 6: Rewrite FASTA headers
    output_entries = []
    for i, (header, seq) in enumerate(entries):
        new_id, new_color = assignments[i]
        new_header = rewrite_header(header, new_id, new_color)
        output_entries.append((new_header, seq))

    write_fasta(output_entries, args.output_fasta)
    print(f"\nClustered FASTA written to: {args.output_fasta}")

    # Step 7: Write TSV mapping
    with open(args.output_tsv, 'w') as f:
        f.write('original_header\tnew_header\tcluster_id\tcluster_color\t'
                'cluster_size\tseq_length\tmin_within_identity\tmax_between_identity\n')

        for i, (orig_header, seq) in enumerate(entries):
            new_id, new_color = assignments[i]
            new_header = rewrite_header(orig_header, new_id, new_color)

            # Compute min within-cluster identity
            same_cluster = [j for j in range(n) if assignments[j][0] == new_id and j != i]
            if same_cluster:
                min_within = min(identity_matrix[i, j] for j in same_cluster)
            else:
                min_within = 100.0  # Singleton

            # Compute max between-cluster identity
            diff_cluster = [j for j in range(n) if assignments[j][0] != new_id]
            if diff_cluster:
                max_between = max(identity_matrix[i, j] for j in diff_cluster)
            else:
                max_between = 0.0

            cluster_size = len(same_cluster) + 1

            f.write(f'{orig_header}\t{new_header}\t{new_id}\t{new_color}\t'
                    f'{cluster_size}\t{len(seq)}\t{min_within:.2f}\t{max_between:.2f}\n')

    print(f"Cluster mapping written to: {args.output_tsv}")

    # Step 8: Diagnostic report (optional)
    if args.diagnostic_report:
        generate_diagnostic_report(
            args.input_fasta, args.diagnostic_report,
            seq_names, assignments, threads=args.threads
        )

    print("\nDone.")


if __name__ == '__main__':
    main()
