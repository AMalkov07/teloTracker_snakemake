#!/usr/bin/env python3
"""
Y prime detection with fragment merging logic.

Y primes are highly variable subtelomeric elements that can occur in tandem repeats.
This module handles:
1. Detection of Y prime fragments
2. Merging adjacent fragments into complete Y primes
3. Classification by size (Short: 5000-6400bp, Long: 6400-6900bp)
"""

import pandas as pd
from typing import List, Dict, Tuple


# Y prime size thresholds (based on repeatmasker_6991_all_y_primes.fasta analysis)
YPRIME_SHORT_MIN = 5000   # Minimum short Y prime size
YPRIME_SHORT_MAX = 6400   # Maximum short Y prime size
YPRIME_LONG_MIN = 6400    # Minimum long Y prime size
YPRIME_LONG_MAX = 6900    # Maximum long Y prime size

# Fragment merging parameters
MAX_FRAGMENT_GAP = 2000   # Maximum gap between fragments to consider merging
MIN_FRAGMENT_SIZE = 1000  # Minimum size to consider as a potential fragment
MIN_YPRIME_IDENTITY = 80.0  # Minimum percent identity for Y prime hits
MIN_YPRIME_COVERAGE = 80.0  # Minimum query coverage for Y prime hits


def classify_yprime_size(length: int) -> str:
    """
    Classify Y prime by size.

    Args:
        length: Y prime length in bp

    Returns:
        'short', 'long', or 'unknown'
    """
    if YPRIME_SHORT_MIN <= length < YPRIME_SHORT_MAX:
        return 'short'
    elif YPRIME_LONG_MIN <= length <= YPRIME_LONG_MAX:
        return 'long'
    else:
        return 'unknown'


def calculate_overlap(start1: int, end1: int, start2: int, end2: int) -> int:
    """
    Calculate overlap between two intervals.

    Args:
        start1, end1: First interval
        start2, end2: Second interval

    Returns:
        Overlap length (0 if no overlap)
    """
    return max(0, min(end1, end2) - max(start1, start2))


def calculate_gap(start1: int, end1: int, start2: int, end2: int) -> int:
    """
    Calculate gap between two intervals.

    Args:
        start1, end1: First interval (should be before second)
        start2, end2: Second interval

    Returns:
        Gap length (0 if overlapping or adjacent)
    """
    if end1 >= start2:
        return 0  # Overlapping or reversed
    return start2 - end1


def merge_yprime_fragments(fragments: List[Dict], chr_name: str) -> List[Dict]:
    """
    Merge adjacent Y prime fragments into complete Y primes.

    Strategy:
    1. Group fragments by chromosome and proximity
    2. Merge fragments that are within MAX_FRAGMENT_GAP bp of each other
    3. Check if merged size matches expected Y prime sizes
    4. Keep individual fragments if they're already valid Y primes
    5. Only merge if result is a valid Y prime size

    Args:
        fragments: List of Y prime hit dictionaries
        chr_name: Chromosome name (for grouping)

    Returns:
        List of merged Y prime regions
    """
    if not fragments:
        return []

    # Sort fragments by position on chromosome
    sorted_frags = sorted(fragments, key=lambda x: (x['chr'], x['start']))

    merged = []
    current_group = [sorted_frags[0]]

    for i in range(1, len(sorted_frags)):
        prev = current_group[-1]
        curr = sorted_frags[i]

        # Check if on same chromosome
        if prev['chr'] != curr['chr']:
            # Process previous group
            merged.extend(process_fragment_group(current_group))
            current_group = [curr]
            continue

        # Check gap between fragments
        gap = calculate_gap(prev['start'], prev['end'], curr['start'], curr['end'])

        if gap <= MAX_FRAGMENT_GAP:
            # Close enough to potentially merge
            current_group.append(curr)
        else:
            # Too far apart - process previous group and start new one
            merged.extend(process_fragment_group(current_group))
            current_group = [curr]

    # Process final group
    merged.extend(process_fragment_group(current_group))

    return merged


def process_fragment_group(group: List[Dict]) -> List[Dict]:
    """
    Process a group of potentially mergeable fragments.

    Strategy:
    1. If single fragment and valid size → keep as is
    2. If single fragment and too small → check if it's a valid partial match
    3. If multiple fragments → try merging
    4. If merged size is valid Y prime → create merged region
    5. Otherwise keep individual fragments that pass thresholds

    Args:
        group: List of fragment dictionaries

    Returns:
        List of processed Y prime regions
    """
    if not group:
        return []

    if len(group) == 1:
        # Single fragment
        frag = group[0]
        length = frag['end'] - frag['start']

        # Check if it's a valid Y prime size
        if length >= MIN_FRAGMENT_SIZE:
            size_class = classify_yprime_size(length)
            if size_class != 'unknown':
                # Valid Y prime
                frag['yprime_class'] = size_class
                frag['merged'] = False
                frag['fragment_count'] = 1
                return [frag]
            elif frag['query_coverage'] >= MIN_YPRIME_COVERAGE:
                # Small but high coverage - might be a complete short Y prime
                frag['yprime_class'] = 'short' if length < YPRIME_SHORT_MAX else 'long'
                frag['merged'] = False
                frag['fragment_count'] = 1
                frag['note'] = f'Small ({length}bp) but high coverage ({frag["query_coverage"]:.1f}%)'
                return [frag]

        # Fragment too small or unknown size
        if frag['pident'] >= MIN_YPRIME_IDENTITY:
            frag['yprime_class'] = 'fragment'
            frag['merged'] = False
            frag['fragment_count'] = 1
            frag['note'] = f'Fragment ({length}bp) - too small for complete Y prime'
            return [frag]

        return []  # Filter out low-quality small fragments

    # Multiple fragments - try merging
    merged_start = min(f['start'] for f in group)
    merged_end = max(f['end'] for f in group)
    merged_length = merged_end - merged_start

    # Calculate merged statistics
    total_aligned_length = sum(f['end'] - f['start'] for f in group)
    avg_pident = sum(f['pident'] * (f['end'] - f['start']) for f in group) / total_aligned_length
    max_bitscore = max(f['bitscore'] for f in group)
    sources = [f['source'] for f in group]

    # Check if merged size is valid
    size_class = classify_yprime_size(merged_length)

    if size_class != 'unknown':
        # Valid merged Y prime!
        merged_region = {
            'chr': group[0]['chr'],
            'start': merged_start,
            'end': merged_end,
            'strand': group[0]['strand'],  # Use first fragment's strand
            'type': 'y_prime',
            'source': f"MERGED_{len(group)}_fragments",
            'pident': avg_pident,
            'bitscore': max_bitscore,
            'query_coverage': total_aligned_length / merged_length * 100,
            'yprime_class': size_class,
            'merged': True,
            'fragment_count': len(group),
            'fragment_sources': sources,
            'note': f'Merged from {len(group)} fragments ({total_aligned_length}bp aligned, {merged_length}bp total)'
        }
        return [merged_region]

    # Merged size not valid - return individual fragments that pass thresholds
    result = []
    for frag in group:
        length = frag['end'] - frag['start']
        if length >= MIN_FRAGMENT_SIZE and frag['pident'] >= MIN_YPRIME_IDENTITY:
            frag['yprime_class'] = 'fragment'
            frag['merged'] = False
            frag['fragment_count'] = 1
            frag['note'] = f'Fragment ({length}bp) - merge failed (total: {merged_length}bp)'
            result.append(frag)

    return result


def filter_overlapping_yprimes(yprimes: List[Dict]) -> List[Dict]:
    """
    Filter overlapping Y primes, keeping the best one.

    When multiple Y primes overlap significantly, keep the one with:
    1. Highest bitscore
    2. Best identity
    3. Longest length

    Args:
        yprimes: List of Y prime dictionaries

    Returns:
        Filtered list with non-overlapping Y primes
    """
    if len(yprimes) <= 1:
        return yprimes

    # Sort by position
    sorted_yp = sorted(yprimes, key=lambda x: (x['chr'], x['start']))

    filtered = []
    skip_indices = set()

    for i in range(len(sorted_yp)):
        if i in skip_indices:
            continue

        current = sorted_yp[i]
        overlapping = [current]

        # Find all overlapping Y primes
        for j in range(i + 1, len(sorted_yp)):
            if sorted_yp[j]['chr'] != current['chr']:
                break

            overlap = calculate_overlap(
                current['start'], current['end'],
                sorted_yp[j]['start'], sorted_yp[j]['end']
            )

            # If overlap is >50% of either region, consider them overlapping
            curr_len = current['end'] - current['start']
            other_len = sorted_yp[j]['end'] - sorted_yp[j]['start']

            if overlap > 0.5 * min(curr_len, other_len):
                overlapping.append(sorted_yp[j])
                skip_indices.add(j)

        # If multiple overlapping, keep the best one
        if len(overlapping) > 1:
            # Sort by: merged (prefer merged), bitscore, identity, length
            best = max(overlapping, key=lambda x: (
                x.get('merged', False),
                x['bitscore'],
                x['pident'],
                x['end'] - x['start']
            ))
            filtered.append(best)
        else:
            filtered.append(current)

    return filtered


def assign_yprimes_to_chr_ends(yprime_df: pd.DataFrame, anchor_regions: Dict) -> Tuple[Dict, List[Dict]]:
    """
    Assign Y primes to chromosome ends based on proximity to anchors.

    Strategy:
    1. For each anchor region, search for Y primes on the same chromosome
    2. Y primes should be telomere-proximal (beyond the anchor)
    3. Multiple Y primes can exist in tandem
    4. No chromosome constraint - any Y prime can be on any chr end

    Args:
        yprime_df: DataFrame with Y prime BLAST results
        anchor_regions: Dict of anchor regions by chr_end

    Returns:
        Tuple of (chr_end_yprimes dict, quality_report list)
    """
    chr_end_yprimes = {}
    quality_report = []

    # Convert Y prime hits to dictionaries
    yprime_hits = []
    for _, row in yprime_df.iterrows():
        hit = {
            'chr': row['sseqid'],
            'start': int(row['start']),
            'end': int(row['end']),
            'strand': row['strand'],
            'type': 'y_prime',
            'source': row['qseqid'],
            'pident': row['pident'],
            'bitscore': row['bitscore'],
            'query_coverage': row['query_coverage']
        }
        yprime_hits.append(hit)

    # Group Y primes by chromosome
    yprimes_by_chr = {}
    for hit in yprime_hits:
        chr_name = hit['chr']
        if chr_name not in yprimes_by_chr:
            yprimes_by_chr[chr_name] = []
        yprimes_by_chr[chr_name].append(hit)

    # Merge fragments for each chromosome
    merged_yprimes_by_chr = {}
    for chr_name, hits in yprimes_by_chr.items():
        merged = merge_yprime_fragments(hits, chr_name)
        filtered = filter_overlapping_yprimes(merged)
        merged_yprimes_by_chr[chr_name] = filtered

    # Assign Y primes to chromosome ends based on anchor positions
    for chr_end, regions in anchor_regions.items():
        anchors = regions.get('anchor', [])
        if not anchors:
            continue

        # Get anchor position
        anchor = anchors[0]  # Should only be one after filtering
        anchor_chr = anchor['chr']
        anchor_start = anchor['start']
        anchor_end = anchor['end']

        # Get Y primes on same chromosome
        chr_yprimes = merged_yprimes_by_chr.get(anchor_chr, [])

        # Determine which arm (L or R) this chr_end is
        arm = chr_end[-1]  # Last character (L or R)

        # Assign Y primes based on position relative to anchor
        assigned_yprimes = []
        for yp in chr_yprimes:
            # For L arm: Y primes should be before (left of) anchor
            # For R arm: Y primes should be after (right of) anchor
            if arm == 'L':
                # Y prime should be telomere-proximal (smaller coordinates)
                if yp['end'] <= anchor_start:
                    assigned_yprimes.append(yp)
            else:  # R arm
                # Y prime should be telomere-proximal (larger coordinates)
                if yp['start'] >= anchor_end:
                    assigned_yprimes.append(yp)

        chr_end_yprimes[chr_end] = assigned_yprimes

        # Report Y primes found
        if assigned_yprimes:
            for yp in assigned_yprimes:
                length = yp['end'] - yp['start']
                quality_report.append({
                    'chr_end': chr_end,
                    'chr': yp['chr'],
                    'start': yp['start'],
                    'end': yp['end'],
                    'length': length,
                    'pident': yp['pident'],
                    'bitscore': yp['bitscore'],
                    'query_coverage': yp.get('query_coverage', 0),
                    'class': yp.get('yprime_class', 'unknown'),
                    'merged': yp.get('merged', False),
                    'fragment_count': yp.get('fragment_count', 1),
                    'note': yp.get('note', '')
                })

    return chr_end_yprimes, quality_report
