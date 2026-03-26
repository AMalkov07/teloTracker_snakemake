#!/usr/bin/env python3
"""
Y prime detection with intelligent fragment merging.

Strategy:
1. Filter for high-quality, substantial BLAST hits (not tiny repetitive matches)
2. Group hits by chromosome and position
3. Merge only 2-3 adjacent large fragments into complete Y primes
4. Classify by size and avoid over-merging repetitive regions
5. Resolve overlaps intelligently: truncate small overlaps, find alternatives for large ones
"""

import pandas as pd
from typing import List, Dict, Tuple, Optional
from copy import deepcopy


# Y prime size thresholds (based on repeatmasker_6991_all_y_primes.fasta analysis)
YPRIME_SHORT_MIN = 5000   # Minimum short Y prime size
YPRIME_SHORT_MAX = 6400   # Maximum short Y prime size
YPRIME_LONG_MIN = 6400    # Minimum long Y prime size
YPRIME_LONG_MAX = 6900    # Maximum long Y prime size

# Detection parameters
MIN_YPRIME_IDENTITY = 80.0  # Minimum percent identity for Y prime hits
MIN_YPRIME_COVERAGE = 50.0  # Minimum query coverage (allow partial hits for fragment detection)
MIN_HIT_LENGTH = 2500       # Minimum hit length to be considered a real fragment (not just repetitive match)
MAX_FRAGMENT_GAP = 500      # Maximum gap between fragments (for telomeric repeats between Y primes)
MAX_FRAGMENTS_TO_MERGE = 3  # Maximum fragments to merge (2-3 real fragments, not hundreds of tiny hits)

# Overlap resolution parameters (configurable)
MAX_TRUNCATION_BP = 100     # Maximum bp to truncate from an alignment to resolve overlap
SMALL_OVERLAP_THRESHOLD = 10  # Overlaps <= this are considered "small"
SMALL_OVERLAP_QUALITY_DIFF = 0.5  # For small overlaps, prefer truncation if quality diff > this %
ALTERNATIVE_QUALITY_THRESHOLD = 2.0  # Consider alternatives within this % identity of original

# Adjusted identity: penalize unaligned bases as if they were mismatches
# This prevents fragments from having artificially high identity scores
USE_ADJUSTED_IDENTITY = True  # Set to False to use standard BLAST identity


def calculate_adjusted_pident(pident: float, alignment_length: int, query_length: int) -> float:
    """
    Calculate adjusted percent identity that penalizes partial alignments.

    Standard BLAST pident = matches / alignment_length
    Adjusted pident = matches / query_length

    This treats unaligned bases as mismatches, penalizing fragments that only
    align part of the query Y prime sequence.

    Example:
        - 6kb Y prime query
        - 3kb alignment with 99.9% identity
        - Standard pident: 99.9%
        - Matches: 0.999 * 3000 = 2997
        - Adjusted pident: 2997 / 6000 = 49.95%

    Args:
        pident: BLAST percent identity (0-100)
        alignment_length: Length of the alignment in bp
        query_length: Total length of the query sequence in bp

    Returns:
        Adjusted percent identity (0-100)
    """
    if query_length <= 0:
        return pident

    # Calculate number of matches from BLAST pident
    matches = (pident / 100.0) * alignment_length

    # Calculate adjusted identity using full query length
    adjusted = (matches / query_length) * 100.0

    # Cap at 100% - alignment gaps can cause adjusted > 100% when
    # alignment_length > query_length, which would incorrectly favor
    # shorter Y prime references over longer ones
    return min(adjusted, 100.0)


def classify_yprime_size(length: int) -> str:
    """Classify Y prime by size."""
    if YPRIME_SHORT_MIN <= length < YPRIME_SHORT_MAX:
        return 'short'
    elif YPRIME_LONG_MIN <= length <= YPRIME_LONG_MAX:
        return 'long'
    else:
        return 'unknown'


def calculate_gap(end1: int, start2: int) -> int:
    """Calculate gap between two intervals."""
    return max(0, start2 - end1)


def filter_high_quality_hits(yprime_df: pd.DataFrame) -> List[Dict]:
    """
    Filter Y prime BLAST hits for substantial, high-quality matches.

    Removes:
    - Tiny repetitive matches (<2.5kb)
    - Low coverage fragments
    - Low identity hits (using adjusted identity if enabled)

    Keeps:
    - Complete Y primes (high coverage, correct size)
    - Large fragments (>2.5kb) that could be merged

    Note: When USE_ADJUSTED_IDENTITY is True, identity is calculated as if
    unaligned bases were mismatches, penalizing partial alignments.
    """
    high_quality_hits = []

    for _, row in yprime_df.iterrows():
        length = row['end'] - row['start']
        pident = row['pident']
        coverage = row['query_coverage']
        qlen = row['qlen']

        # Calculate adjusted identity if enabled
        if USE_ADJUSTED_IDENTITY:
            adjusted_pident = calculate_adjusted_pident(pident, length, qlen)
        else:
            adjusted_pident = pident

        # Skip low quality hits (using adjusted identity)
        if adjusted_pident < MIN_YPRIME_IDENTITY:
            continue

        # Skip tiny hits (just repetitive DNA, not real Y prime fragments)
        if length < MIN_HIT_LENGTH:
            continue

        # Skip very low coverage (unless it's a large fragment)
        if coverage < MIN_YPRIME_COVERAGE and length < 4000:
            continue

        hit = {
            'chr': row['sseqid'],
            'start': int(row['start']),
            'end': int(row['end']),
            'length': length,
            'strand': row['strand'],
            'type': 'y_prime',
            'source': row['qseqid'],
            'pident': pident,  # Keep original BLAST identity
            'adjusted_pident': adjusted_pident,  # Store adjusted identity
            'bitscore': row['bitscore'],
            'query_coverage': coverage,
            'qlen': qlen
        }
        high_quality_hits.append(hit)

    return high_quality_hits


def merge_adjacent_fragments(hits: List[Dict]) -> List[Dict]:
    """
    Merge 2-3 adjacent large fragments into complete Y primes.

    Strategy:
    - Group hits that are close together (<500bp apart)
    - Only merge if result makes a valid Y prime size
    - Limit to 2-3 fragments max (not hundreds)
    - Keep complete Y primes as-is
    """
    if not hits:
        return []

    # Sort by position
    sorted_hits = sorted(hits, key=lambda x: (x['chr'], x['start']))

    merged_results = []
    i = 0

    while i < len(sorted_hits):
        current = sorted_hits[i]
        current_length = current['end'] - current['start']

        # Check if this hit is already a valid Y prime size
        size_class = classify_yprime_size(current_length)

        if size_class != 'unknown' or current['query_coverage'] >= 95:
            # Already a complete Y prime - keep as is
            current['yprime_class'] = size_class if size_class != 'unknown' else 'complete'
            current['merged'] = False
            current['fragment_count'] = 1
            merged_results.append(current)
            i += 1
            continue

        # Try to merge with next 1-2 hits
        merge_group = [current]
        j = i + 1

        while j < len(sorted_hits) and len(merge_group) < MAX_FRAGMENTS_TO_MERGE:
            next_hit = sorted_hits[j]

            # Must be on same chromosome
            if next_hit['chr'] != current['chr']:
                break

            # Check gap between fragments
            prev_end = merge_group[-1]['end']
            gap = calculate_gap(prev_end, next_hit['start'])

            if gap > MAX_FRAGMENT_GAP:
                break  # Too far apart

            merge_group.append(next_hit)
            j += 1

        # Try merging if we have 2+ fragments
        if len(merge_group) >= 2:
            merged_start = merge_group[0]['start']
            merged_end = merge_group[-1]['end']
            merged_length = merged_end - merged_start

            # Check if merged size is valid
            size_class = classify_yprime_size(merged_length)

            if size_class != 'unknown':
                # Valid merged Y prime!
                total_aligned = sum(f['length'] for f in merge_group)
                avg_pident = sum(f['pident'] * f['length'] for f in merge_group) / total_aligned
                max_bitscore = max(f['bitscore'] for f in merge_group)
                sources = [f['source'] for f in merge_group]

                # For merged fragments, calculate adjusted identity based on combined alignment
                # Use the maximum qlen from fragments as the reference Y prime size
                max_qlen = max(f.get('qlen', merged_length) for f in merge_group)

                # Calculate adjusted identity: total matches / max query length
                # Total matches = sum of (pident * aligned_length) for each fragment
                total_matches = sum((f['pident'] / 100.0) * f['length'] for f in merge_group)
                if USE_ADJUSTED_IDENTITY:
                    adjusted_pident = (total_matches / max_qlen) * 100.0
                else:
                    adjusted_pident = avg_pident

                merged_hit = {
                    'chr': merge_group[0]['chr'],
                    'start': merged_start,
                    'end': merged_end,
                    'length': merged_length,
                    'strand': merge_group[0]['strand'],
                    'type': 'y_prime',
                    'source': f"MERGED_{len(merge_group)}_fragments",
                    'pident': avg_pident,
                    'adjusted_pident': adjusted_pident,
                    'bitscore': max_bitscore,
                    'query_coverage': (total_aligned / merged_length) * 100,
                    'qlen': max_qlen,
                    'yprime_class': size_class,
                    'merged': True,
                    'fragment_count': len(merge_group),
                    'fragment_sources': sources,
                    'note': 'Merged {} fragments: {}'.format(len(merge_group), ", ".join([str(f["length"]) + "bp" for f in merge_group]))
                }
                merged_results.append(merged_hit)
                i = j  # Skip merged hits
                continue

        # Couldn't merge or merge didn't make valid size - keep original
        current['yprime_class'] = 'fragment'
        current['merged'] = False
        current['fragment_count'] = 1
        current['note'] = f'{current_length}bp fragment - too small or couldn\'t merge'
        merged_results.append(current)
        i += 1

    return merged_results


def calculate_overlap(yp1: Dict, yp2: Dict) -> int:
    """Calculate overlap in bp between two Y primes."""
    if yp1['chr'] != yp2['chr']:
        return 0
    overlap_start = max(yp1['start'], yp2['start'])
    overlap_end = min(yp1['end'], yp2['end'])
    return max(0, overlap_end - overlap_start)


def get_quality_score(yp: Dict) -> float:
    """
    Get a quality score for ranking Y primes.

    Uses adjusted_pident if available and USE_ADJUSTED_IDENTITY is enabled,
    otherwise falls back to standard pident.

    The adjusted_pident is ALWAYS the primary factor for selection. A higher
    adjusted_pident will always win, regardless of alignment length.

    Length is ONLY used as a tiebreaker when two Y primes have the exact same
    adjusted_pident. In that case, the longer alignment is preferred.

    The length bonus is intentionally tiny (0.00000001 per bp) so that even
    a 10,000 bp length difference (0.0001) cannot overcome a 0.001% identity
    difference. This ensures identity is strictly prioritized.
    """
    if USE_ADJUSTED_IDENTITY and 'adjusted_pident' in yp:
        base_score = yp['adjusted_pident']
    else:
        base_score = yp['pident']

    # Add an extremely small length bonus as a tiebreaker ONLY
    # This bonus is so small it can never affect identity-based ranking:
    # - Maximum realistic length: 10,000 bp
    # - Maximum bonus: 10,000 * 0.00000001 = 0.0001
    # - Minimum identity difference that matters: 0.001% (3 decimal places)
    # So the length bonus can ONLY matter when identities are exactly equal
    length = yp.get('length', 0)
    length_bonus = length * 0.00000001

    return base_score + length_bonus


def truncate_yprime(yp: Dict, truncate_start: bool, truncate_bp: int) -> Optional[Dict]:
    """
    Truncate a Y prime by removing bp from start or end.

    Returns None if truncation would make Y prime invalid (< YPRIME_SHORT_MIN).
    """
    new_yp = deepcopy(yp)

    if truncate_start:
        new_yp['start'] = yp['start'] + truncate_bp
    else:
        new_yp['end'] = yp['end'] - truncate_bp

    new_yp['length'] = new_yp['end'] - new_yp['start']

    # Check if still valid
    if new_yp['length'] < YPRIME_SHORT_MIN:
        return None

    # Update size class
    new_yp['yprime_class'] = classify_yprime_size(new_yp['length'])
    new_yp['note'] = new_yp.get('note', '') + f' [truncated {truncate_bp}bp from {"start" if truncate_start else "end"}]'

    return new_yp


def find_alternative_alignment(
    target_yp: Dict,
    all_candidates: List[Dict],
    existing_yprimes: List[Dict],
    quality_threshold: float
) -> Optional[Dict]:
    """
    Find an alternative alignment for a Y prime that doesn't overlap with existing ones.

    Args:
        target_yp: The Y prime we need an alternative for
        all_candidates: All candidate alignments (including ones not yet selected)
        existing_yprimes: Y primes already in the final set
        quality_threshold: Maximum identity difference to consider as alternative

    Returns:
        Alternative Y prime dict, or None if no suitable alternative found
    """
    target_quality = get_quality_score(target_yp)
    target_region = (target_yp['chr'], target_yp['start'], target_yp['end'])

    best_alternative = None
    best_alt_quality = 0

    for candidate in all_candidates:
        # Skip if it's the same alignment
        if (candidate['chr'] == target_yp['chr'] and
            candidate['start'] == target_yp['start'] and
            candidate['end'] == target_yp['end']):
            continue

        # Must be on same chromosome
        if candidate['chr'] != target_yp['chr']:
            continue

        # Check quality threshold
        candidate_quality = get_quality_score(candidate)
        if target_quality - candidate_quality > quality_threshold:
            continue

        # Check if candidate overlaps with target region (should cover similar area)
        # Allow alternatives that cover at least 50% of the same region
        overlap_with_target = calculate_overlap(candidate, target_yp)
        min_overlap_needed = 0.3 * min(candidate['length'], target_yp['length'])
        if overlap_with_target < min_overlap_needed:
            continue

        # Check if candidate overlaps with any existing Y primes
        has_overlap = False
        for existing in existing_yprimes:
            if existing['chr'] != candidate['chr']:
                continue
            overlap = calculate_overlap(candidate, existing)
            if overlap > 0:
                has_overlap = True
                break

        if not has_overlap and candidate_quality > best_alt_quality:
            best_alternative = candidate
            best_alt_quality = candidate_quality

    return best_alternative


def resolve_overlapping_yprimes(yprimes: List[Dict], all_candidates: List[Dict]) -> List[Dict]:
    """
    Resolve overlapping Y primes using intelligent truncation and alternative selection.

    Strategy:
    1. For overlaps <= 10bp: truncate if quality diff > 0.5%, else prefer non-overlapping alternative
    2. For overlaps 11-100bp: prefer non-overlapping alternative within 2% identity, else truncate
    3. For overlaps > 100bp: find alternative or drop lower quality alignment

    Args:
        yprimes: List of Y primes after merging (may have overlaps)
        all_candidates: All high-quality hits (for finding alternatives)

    Returns:
        List of non-overlapping Y primes
    """
    if len(yprimes) <= 1:
        return yprimes

    # Sort by quality (best first) then by position
    sorted_yp = sorted(yprimes, key=lambda x: (-get_quality_score(x), x['chr'], x['start']))

    final_yprimes = []
    processed_indices = set()

    for i, current in enumerate(sorted_yp):
        if i in processed_indices:
            continue

        # Check for overlaps with already-selected Y primes
        overlaps_with_selected = []
        for selected in final_yprimes:
            overlap = calculate_overlap(current, selected)
            if overlap > 0:
                overlaps_with_selected.append((selected, overlap))

        if not overlaps_with_selected:
            # No overlap - add directly
            final_yprimes.append(deepcopy(current))
            processed_indices.add(i)
            continue

        # Handle overlaps
        for selected, overlap_bp in overlaps_with_selected:
            current_quality = get_quality_score(current)
            selected_quality = get_quality_score(selected)
            quality_diff = abs(current_quality - selected_quality)

            # Determine which side to truncate
            # If current starts after selected starts, truncate current's start
            # If current ends before selected ends, truncate current's end
            if current['start'] >= selected['start']:
                truncate_start = True
                truncate_amount = selected['end'] - current['start']
            else:
                truncate_start = False
                truncate_amount = current['end'] - selected['start']

            truncate_amount = max(0, truncate_amount)

            # Decision logic based on overlap size
            if overlap_bp <= SMALL_OVERLAP_THRESHOLD:
                # Small overlap (≤10bp)
                if quality_diff > SMALL_OVERLAP_QUALITY_DIFF:
                    # Quality difference is significant - truncate
                    truncated = truncate_yprime(current, truncate_start, truncate_amount)
                    if truncated:
                        current = truncated
                    else:
                        # Truncation would make it invalid - look for alternative
                        alt = find_alternative_alignment(current, all_candidates, final_yprimes, ALTERNATIVE_QUALITY_THRESHOLD)
                        if alt:
                            current = alt
                        else:
                            current = None  # Drop it
                            break
                else:
                    # Similar quality - prefer non-overlapping alternative if available
                    alt = find_alternative_alignment(current, all_candidates, final_yprimes, SMALL_OVERLAP_QUALITY_DIFF)
                    if alt:
                        current = alt
                    else:
                        # No alternative - truncate anyway
                        truncated = truncate_yprime(current, truncate_start, truncate_amount)
                        if truncated:
                            current = truncated
                        else:
                            current = None
                            break

            elif overlap_bp <= MAX_TRUNCATION_BP:
                # Medium overlap (11-100bp) - prefer alternative within 2%
                alt = find_alternative_alignment(current, all_candidates, final_yprimes, ALTERNATIVE_QUALITY_THRESHOLD)
                if alt:
                    current = alt
                else:
                    # No alternative - truncate
                    truncated = truncate_yprime(current, truncate_start, truncate_amount)
                    if truncated:
                        current = truncated
                    else:
                        current = None
                        break

            else:
                # Large overlap (>100bp) - must find alternative or drop
                alt = find_alternative_alignment(current, all_candidates, final_yprimes, 100.0)  # Accept any quality alternative
                if alt:
                    current = alt
                else:
                    # No alternative - drop this Y prime
                    print(f"      Dropping Y prime at {current['chr']}:{current['start']}-{current['end']} "
                          f"due to {overlap_bp}bp overlap with no alternative")
                    current = None
                    break

        # Add current if it survived
        if current is not None:
            # Final check - make sure no overlap remains
            still_overlaps = False
            for selected in final_yprimes:
                if calculate_overlap(current, selected) > 0:
                    still_overlaps = True
                    break

            if not still_overlaps:
                final_yprimes.append(current)

        processed_indices.add(i)

    # Sort final results by position
    final_yprimes.sort(key=lambda x: (x['chr'], x['start']))

    return final_yprimes


def filter_overlapping_yprimes(yprimes: List[Dict], all_candidates: List[Dict] = None) -> List[Dict]:
    """
    Wrapper for backward compatibility - resolve overlapping Y primes.

    If all_candidates not provided, falls back to simple best-selection logic.
    """
    if all_candidates is not None:
        return resolve_overlapping_yprimes(yprimes, all_candidates)

    # Fallback: simple overlap resolution (keep best when heavily overlapping)
    if len(yprimes) <= 1:
        return yprimes

    sorted_yp = sorted(yprimes, key=lambda x: (x['chr'], x['start']))
    filtered = []
    skip_indices = set()

    for i in range(len(sorted_yp)):
        if i in skip_indices:
            continue

        current = sorted_yp[i]
        overlapping = [current]

        for j in range(i + 1, len(sorted_yp)):
            if sorted_yp[j]['chr'] != current['chr']:
                break

            other = sorted_yp[j]
            overlap = calculate_overlap(current, other)

            # If overlap > 50% of smaller, consider them overlapping
            if overlap > 0.5 * min(current['length'], other['length']):
                overlapping.append(other)
                skip_indices.add(j)

        if len(overlapping) > 1:
            best = max(overlapping, key=lambda x: (
                x.get('yprime_class', 'unknown') != 'unknown',
                x.get('yprime_class', 'unknown') != 'fragment',
                x.get('merged', False),
                x['bitscore'],
                x['pident']
            ))
            filtered.append(best)
        else:
            filtered.append(current)

    return filtered


def assign_yprimes_to_chr_ends(yprime_df: pd.DataFrame, anchor_regions: Dict) -> Tuple[Dict, List[Dict]]:
    """
    Assign Y primes to chromosome ends based on proximity to anchors.

    Y primes should be telomere-proximal (beyond the anchor toward the telomere).
    """
    # Step 1: Filter for high-quality hits only
    identity_type = "adjusted identity" if USE_ADJUSTED_IDENTITY else "identity"
    print(f"    Filtering Y prime hits for substantial matches (>2.5kb, >80% {identity_type})...")
    high_quality_hits = filter_high_quality_hits(yprime_df)
    print(f"    Found {len(high_quality_hits)} high-quality Y prime hits (filtered from {len(yprime_df)} total)")
    if USE_ADJUSTED_IDENTITY:
        print("    (Using adjusted identity: unaligned bases penalized as mismatches)")

    # Step 2: Group by chromosome
    hits_by_chr = {}
    for hit in high_quality_hits:
        chr_name = hit['chr']
        if chr_name not in hits_by_chr:
            hits_by_chr[chr_name] = []
        hits_by_chr[chr_name].append(hit)

    # Step 3: Merge adjacent fragments for each chromosome
    print("    Merging adjacent fragments into complete Y primes...")
    merged_by_chr = {}
    for chr_name, hits in hits_by_chr.items():
        merged = merge_adjacent_fragments(hits)
        # Pass all high-quality hits for this chromosome as candidates for alternative selection
        filtered = filter_overlapping_yprimes(merged, all_candidates=hits)
        merged_by_chr[chr_name] = filtered

    total_after_merge = sum(len(yps) for yps in merged_by_chr.values())
    print(f"    After merging and overlap resolution: {total_after_merge} Y prime regions detected")

    # Step 4: Assign to chromosome ends based on anchor positions
    chr_end_yprimes = {}
    quality_report = []

    for chr_end, regions in anchor_regions.items():
        anchors = regions.get('anchor', [])
        if not anchors:
            continue

        anchor = anchors[0]
        anchor_chr = anchor['chr']
        anchor_start = anchor['start']
        anchor_end = anchor['end']
        arm = chr_end[-1]  # L or R

        # Get Y primes on same chromosome
        chr_yprimes = merged_by_chr.get(anchor_chr, [])

        # Assign based on arm
        assigned_yprimes = []
        for yp in chr_yprimes:
            if arm == 'L':
                # Y primes should be telomere-proximal (before anchor)
                if yp['end'] <= anchor_start:
                    assigned_yprimes.append(yp)
            else:  # R arm
                # Y primes should be telomere-proximal (after anchor)
                if yp['start'] >= anchor_end:
                    assigned_yprimes.append(yp)

        chr_end_yprimes[chr_end] = assigned_yprimes

        # Add to quality report
        for yp in assigned_yprimes:
            quality_report.append({
                'chr_end': chr_end,
                'chr': yp['chr'],
                'start': yp['start'],
                'end': yp['end'],
                'length': yp['length'],
                'pident': yp['pident'],
                'adjusted_pident': yp.get('adjusted_pident', yp['pident']),
                'bitscore': yp['bitscore'],
                'query_coverage': yp.get('query_coverage', 0),
                'class': yp.get('yprime_class', 'unknown'),
                'merged': yp.get('merged', False),
                'fragment_count': yp.get('fragment_count', 1),
                'note': yp.get('note', '')
            })

    return chr_end_yprimes, quality_report
