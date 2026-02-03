#!/usr/bin/env python3
"""
Adjust feature boundaries to maximize ITS (Interstitial Telomeric Sequence) regions.

ITS sequences are telomere-like:
- L arm (left side): AC-rich sequences
- R arm (right side): TG-rich sequences

Algorithm:
1. Plan which features to use (anchors, X primes, Y primes)
2. For each feature boundary, check if there are non-telomeric bases within 5bp
3. If yes: EXPAND the boundary until no non-telomeric bases within 5bp (or hit next feature)
4. If no expansion occurred: TRIM telomeric bases to maximize ITS size
"""

import argparse
import os
from Bio import SeqIO
from typing import Dict, List, Tuple, Optional


def get_telomeric_bases(arm: str) -> set:
    """
    Get the bases that are considered telomeric for this arm.

    Args:
        arm: 'L' or 'R'

    Returns:
        Set of telomeric bases
    """
    if arm == 'R':
        return {'T', 'G'}  # TG repeats on R arm (toward high coordinates/telomere)
    else:  # L arm
        return {'A', 'C'}  # AC repeats on L arm (toward low coordinates/telomere)


def get_non_telomeric_bases(arm: str) -> set:
    """Get bases that are NOT telomeric for this arm."""
    if arm == 'R':
        return {'A', 'C'}
    else:
        return {'T', 'G'}


def has_non_telomeric_in_window(seq: str, arm: str, window: int = 5) -> bool:
    """
    Check if there are non-telomeric bases in the sequence window.

    Args:
        seq: Sequence to check
        arm: 'L' or 'R'
        window: Window size to check

    Returns:
        True if non-telomeric bases exist in window
    """
    telomeric = get_telomeric_bases(arm)
    check_seq = seq[:window].upper() if len(seq) >= window else seq.upper()

    for base in check_seq:
        if base in 'ACGT' and base not in telomeric:
            return True
    return False


def expand_start_boundary(ref_seq: str, current_start: int, arm: str,
                          min_boundary: int, window: int = 5) -> Tuple[int, bool]:
    """
    Expand the start boundary backwards until no non-telomeric bases within window.

    Args:
        ref_seq: Reference sequence
        current_start: Current start position (0-based)
        arm: 'L' or 'R'
        min_boundary: Minimum position we can expand to (e.g., end of previous feature)
        window: Window size to check

    Returns:
        Tuple of (new_start, did_expand)
    """
    telomeric = get_telomeric_bases(arm)
    new_start = current_start
    did_expand = False

    while new_start > min_boundary:
        # Check window at current position
        check_seq = ref_seq[new_start:new_start + window].upper()

        # Check if all bases in window are telomeric
        has_non_telo = any(b in 'ACGT' and b not in telomeric for b in check_seq)

        if not has_non_telo:
            # No non-telomeric bases - stop expanding
            break

        # Found non-telomeric base - expand backwards
        new_start -= 1
        did_expand = True

    return new_start, did_expand


def expand_end_boundary(ref_seq: str, current_end: int, arm: str,
                        max_boundary: int, window: int = 5) -> Tuple[int, bool]:
    """
    Expand the end boundary forwards until no non-telomeric bases within window.

    Args:
        ref_seq: Reference sequence
        current_end: Current end position (0-based, exclusive)
        arm: 'L' or 'R'
        max_boundary: Maximum position we can expand to (e.g., start of next feature)
        window: Window size to check

    Returns:
        Tuple of (new_end, did_expand)
    """
    telomeric = get_telomeric_bases(arm)
    new_end = current_end
    did_expand = False

    while new_end < max_boundary:
        # Check window ending at current position
        check_start = max(0, new_end - window)
        check_seq = ref_seq[check_start:new_end].upper()

        # Check if all bases in window are telomeric
        has_non_telo = any(b in 'ACGT' and b not in telomeric for b in check_seq)

        if not has_non_telo:
            # No non-telomeric bases - stop expanding
            break

        # Found non-telomeric base - expand forwards
        new_end += 1
        did_expand = True

    return new_end, did_expand


def trim_start_boundary(ref_seq: str, current_start: int, current_end: int, arm: str) -> int:
    """
    Trim telomeric bases from start to maximize ITS.

    Args:
        ref_seq: Reference sequence
        current_start: Current start position (0-based)
        current_end: Current end position (0-based, exclusive)
        arm: 'L' or 'R'

    Returns:
        New start position
    """
    telomeric = get_telomeric_bases(arm)
    new_start = current_start

    # Trim while we see telomeric bases
    while new_start < current_end - 10:  # Keep at least 10bp of feature
        base = ref_seq[new_start].upper()
        if base in telomeric:
            new_start += 1
        else:
            break

    return new_start


def trim_end_boundary(ref_seq: str, current_start: int, current_end: int, arm: str) -> int:
    """
    Trim telomeric bases from end to maximize ITS.

    Args:
        ref_seq: Reference sequence
        current_start: Current start position (0-based)
        current_end: Current end position (0-based, exclusive)
        arm: 'L' or 'R'

    Returns:
        New end position
    """
    telomeric = get_telomeric_bases(arm)
    new_end = current_end

    # Trim while we see telomeric bases
    while new_end > current_start + 10:  # Keep at least 10bp of feature
        base = ref_seq[new_end - 1].upper()
        if base in telomeric:
            new_end -= 1
        else:
            break

    return new_end


def get_feature_order(features: List[Dict], arm: str) -> List[Dict]:
    """
    Order features from anchor toward telomere.

    Args:
        features: List of feature dicts with 'start', 'end', 'type' keys
        arm: 'L' or 'R'

    Returns:
        Sorted list of features
    """
    if arm == 'R':
        # R arm: anchor at low coords, telomere at high coords
        # Order: anchor -> x_core -> x_variable -> y_primes (sorted by position)
        return sorted(features, key=lambda x: x['start'])
    else:
        # L arm: telomere at low coords, anchor at high coords
        # Order from anchor toward telomere: anchor -> x_core -> x_variable -> y_primes
        # But in coordinate order, y_primes come first (low coords)
        return sorted(features, key=lambda x: x['start'], reverse=True)


def adjust_feature_boundaries(ref_seq: str, features: List[Dict], arm: str,
                              chr_start: int = 0, chr_end: int = None,
                              window: int = 5) -> List[Dict]:
    """
    Adjust all feature boundaries to maximize ITS regions.

    Args:
        ref_seq: Reference sequence for this chromosome
        features: List of feature dicts with 'start', 'end', 'type', 'name' keys
        arm: 'L' or 'R'
        chr_start: Chromosome start (for boundary limits)
        chr_end: Chromosome end (for boundary limits)
        window: Window size for checking non-telomeric bases

    Returns:
        List of adjusted feature dicts
    """
    if chr_end is None:
        chr_end = len(ref_seq)

    if not features:
        return features

    # Sort features by coordinate
    sorted_features = sorted(features, key=lambda x: x['start'])
    adjusted_features = []

    for i, feature in enumerate(sorted_features):
        # Determine boundaries for expansion (next/previous feature or chromosome ends)
        if i == 0:
            if arm == 'R':
                # First feature on R arm - previous boundary is anchor end or chr start
                prev_boundary = chr_start
            else:
                # First feature on L arm (sorted, so lowest coord) - previous is chr start (telomere)
                prev_boundary = chr_start
        else:
            prev_boundary = sorted_features[i - 1]['end']

        if i == len(sorted_features) - 1:
            if arm == 'R':
                # Last feature on R arm - next boundary is chr end (telomere)
                next_boundary = chr_end
            else:
                # Last feature on L arm - next boundary is anchor or chr end
                next_boundary = chr_end
        else:
            next_boundary = sorted_features[i + 1]['start']

        # Get current boundaries
        current_start = feature['start']
        current_end = feature['end']

        # Check and adjust start boundary
        start_seq = ref_seq[current_start:current_start + window]
        needs_start_expand = has_non_telomeric_in_window(start_seq, arm, window)

        if needs_start_expand:
            new_start, _ = expand_start_boundary(ref_seq, current_start, arm, prev_boundary, window)
        else:
            # No expansion needed - trim to maximize ITS
            new_start = trim_start_boundary(ref_seq, current_start, current_end, arm)

        # Check and adjust end boundary
        end_seq = ref_seq[max(0, current_end - window):current_end]
        needs_end_expand = has_non_telomeric_in_window(end_seq[::-1], arm, window)  # Reverse to check from end

        if needs_end_expand:
            new_end, _ = expand_end_boundary(ref_seq, current_end, arm, next_boundary, window)
        else:
            # No expansion needed - trim to maximize ITS
            new_end = trim_end_boundary(ref_seq, current_start, current_end, arm)

        # Create adjusted feature
        adjusted_feature = feature.copy()
        adjusted_feature['original_start'] = current_start
        adjusted_feature['original_end'] = current_end
        adjusted_feature['start'] = new_start
        adjusted_feature['end'] = new_end
        adjusted_feature['start_adjusted'] = new_start != current_start
        adjusted_feature['end_adjusted'] = new_end != current_end

        adjusted_features.append(adjusted_feature)

    return adjusted_features


def adjust_chr_end_features(ref_sequences: Dict[str, str],
                            chr_end_regions: Dict[str, Dict],
                            window: int = 5) -> Dict[str, Dict]:
    """
    Adjust feature boundaries for all chromosome ends.

    Args:
        ref_sequences: Dict mapping chr name to sequence
        chr_end_regions: Dict mapping chr_end to {'anchor': [...], 'x_prime': [...], 'y_prime': [...]}
        window: Window size for boundary checking

    Returns:
        Adjusted chr_end_regions
    """
    adjusted_regions = {}

    for chr_end, regions in chr_end_regions.items():
        arm = chr_end[-1]  # L or R

        # Get chromosome name and sequence
        # Try to find matching chromosome in ref_sequences
        chr_name = None
        ref_seq = None

        for ref_chr, seq in ref_sequences.items():
            # Match chr1, chr1_extended, chr1L -> chr1
            ref_chr_base = ref_chr.replace('_extended', '').lower()
            chr_end_base = chr_end[:-1].lower()  # Remove L/R

            if ref_chr_base == chr_end_base or ref_chr_base.startswith(chr_end_base):
                chr_name = ref_chr
                ref_seq = seq
                break

        if ref_seq is None:
            print(f"Warning: Could not find reference sequence for {chr_end}")
            adjusted_regions[chr_end] = regions
            continue

        chr_len = len(ref_seq)

        # Collect all features for this chr_end
        all_features = []

        for anchor in regions.get('anchor', []):
            all_features.append({**anchor, 'feature_type': 'anchor'})

        for xprime in regions.get('x_prime', []):
            all_features.append({**xprime, 'feature_type': 'x_prime'})

        for i, yprime in enumerate(regions.get('y_prime', [])):
            all_features.append({**yprime, 'feature_type': 'y_prime', 'yprime_index': i})

        # Adjust boundaries
        adjusted_features = adjust_feature_boundaries(
            ref_seq, all_features, arm,
            chr_start=0, chr_end=chr_len,
            window=window
        )

        # Reconstruct chr_end_regions structure
        adjusted_regions[chr_end] = {
            'anchor': [],
            'x_prime': [],
            'y_prime': []
        }

        for feat in adjusted_features:
            feat_type = feat.pop('feature_type')
            if 'yprime_index' in feat:
                del feat['yprime_index']
            adjusted_regions[chr_end][feat_type].append(feat)

    return adjusted_regions


def load_reference_sequences(fasta_path: str) -> Dict[str, str]:
    """Load reference sequences from FASTA file."""
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def main():
    parser = argparse.ArgumentParser(
        description='Adjust feature boundaries to maximize ITS regions'
    )
    parser.add_argument(
        '--reference', '-r',
        required=True,
        help='Reference FASTA file'
    )
    parser.add_argument(
        '--input-bed', '-i',
        required=True,
        help='Input BED file with features'
    )
    parser.add_argument(
        '--output-bed', '-o',
        required=True,
        help='Output BED file with adjusted features'
    )
    parser.add_argument(
        '--window', '-w',
        type=int,
        default=5,
        help='Window size for checking non-telomeric bases (default: 5)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print verbose output'
    )

    args = parser.parse_args()

    print(f"Loading reference from {args.reference}...")
    ref_sequences = load_reference_sequences(args.reference)
    print(f"Loaded {len(ref_sequences)} chromosomes")

    # Parse input BED file
    print(f"Loading features from {args.input_bed}...")
    features_by_chr = {}

    with open(args.input_bed, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            chr_name = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            strand = parts[4] if len(parts) > 4 else '+'

            if chr_name not in features_by_chr:
                features_by_chr[chr_name] = []

            features_by_chr[chr_name].append({
                'chr': chr_name,
                'start': start,
                'end': end,
                'name': name,
                'strand': strand
            })

    print(f"Loaded features for {len(features_by_chr)} chromosomes")

    # Adjust boundaries
    print(f"Adjusting feature boundaries (window={args.window})...")
    adjusted_count = 0

    with open(args.output_bed, 'w') as out:
        for chr_name, features in sorted(features_by_chr.items()):
            # Get reference sequence
            ref_seq = None
            for ref_chr, seq in ref_sequences.items():
                if ref_chr == chr_name or ref_chr.startswith(chr_name):
                    ref_seq = seq
                    break

            if ref_seq is None:
                print(f"Warning: No reference for {chr_name}, keeping original")
                for feat in features:
                    length = feat['end'] - feat['start']
                    out.write(f"{feat['chr']}\t{feat['start']}\t{feat['end']}\t"
                             f"{feat['name']}\t{feat['strand']}\t{length}\n")
                continue

            # Determine arm from feature names
            for feat in features:
                name = feat['name']
                if 'L_' in name or name.endswith('L'):
                    arm = 'L'
                elif 'R_' in name or name.endswith('R'):
                    arm = 'R'
                else:
                    # Try to infer from position
                    arm = 'L' if feat['start'] < len(ref_seq) / 2 else 'R'

                feat['arm'] = arm

            # Group by arm and adjust
            for arm in ['L', 'R']:
                arm_features = [f for f in features if f.get('arm') == arm]
                if not arm_features:
                    continue

                adjusted = adjust_feature_boundaries(
                    ref_seq, arm_features, arm,
                    window=args.window
                )

                for feat in adjusted:
                    if feat.get('start_adjusted') or feat.get('end_adjusted'):
                        adjusted_count += 1
                        if args.verbose:
                            print(f"  {feat['name']}: {feat.get('original_start', feat['start'])}-"
                                  f"{feat.get('original_end', feat['end'])} -> "
                                  f"{feat['start']}-{feat['end']}")

                    length = feat['end'] - feat['start']
                    out.write(f"{feat['chr']}\t{feat['start']}\t{feat['end']}\t"
                             f"{feat['name']}\t{feat['strand']}\t{length}\n")

    print(f"Adjusted {adjusted_count} feature boundaries")
    print(f"Output written to {args.output_bed}")


if __name__ == '__main__':
    main()
