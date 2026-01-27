#!/usr/bin/env python3
"""
X prime detection with adjusted identity scoring.

Strategy:
1. BLAST all X prime sequences against the reference
2. For each chromosome end, find the best matching X prime
3. Require full-length matches (penalize fragments with adjusted identity)
4. Warn if best match is from unexpected chromosome end
5. Validate X prime position relative to Y primes (X prime should be closer to chr middle)
"""

import pandas as pd
from typing import List, Dict, Tuple, Optional
from collections import defaultdict


# X prime size thresholds (based on 6991_xprimes.fasta analysis)
XPRIME_MIN_SIZE = 100    # Minimum X prime size (core-only cases)
XPRIME_MAX_SIZE = 900    # Maximum X prime size (core + variable)

# Detection parameters
MIN_XPRIME_IDENTITY = 80.0       # Minimum percent identity for X prime hits
MIN_XPRIME_ADJUSTED_IDENTITY = 75.0  # Minimum adjusted identity (penalizes fragments)
MIN_HIGH_QUALITY_IDENTITY = 95.0  # Threshold for high-quality matches
MIN_XPRIME_COVERAGE = 80.0       # Minimum query coverage for valid match


def calculate_adjusted_pident(pident: float, alignment_length: int, query_length: int) -> float:
    """
    Calculate adjusted percent identity that penalizes partial alignments.

    Standard BLAST pident = matches / alignment_length
    Adjusted pident = matches / query_length

    This treats unaligned bases as mismatches, penalizing fragments.

    Args:
        pident: BLAST percent identity (0-100)
        alignment_length: Length of the alignment in bp
        query_length: Total length of the query sequence in bp

    Returns:
        Adjusted percent identity (0-100)
    """
    if query_length <= 0:
        return pident

    matches = (pident / 100.0) * alignment_length
    adjusted = (matches / query_length) * 100.0

    return adjusted


def extract_chr_end_from_xprime_query(qseqid: str) -> Optional[str]:
    """
    Extract chromosome end from X prime query sequence ID.

    Examples:
        chr1L_xprime -> chr1L
        chr10R_xprime -> chr10R

    Args:
        qseqid: Query sequence identifier

    Returns:
        Chromosome end (e.g., "chr1L") or None
    """
    if '_xprime' in qseqid:
        return qseqid.replace('_xprime', '')
    return None


def parse_xprime_fasta_for_boundaries(xprime_fasta_path: str) -> Dict[str, Dict]:
    """
    Parse X prime FASTA file to extract core/variable boundary information.

    Header format: >chr1L_xprime chr1:213-962(-) len=749 core_len=463 var_len=287

    Args:
        xprime_fasta_path: Path to X prime FASTA file

    Returns:
        Dict mapping chr_end to {'core_len': int, 'var_len': int, 'strand': str}
    """
    boundaries = {}

    try:
        with open(xprime_fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Parse header
                    parts = line.strip().split()
                    if len(parts) >= 1:
                        # Extract chr_end from first part (e.g., >chr1L_xprime -> chr1L)
                        chr_end = parts[0][1:].replace('_xprime', '')

                        # Parse core_len and var_len from header
                        core_len = 0
                        var_len = 0
                        strand = '+'

                        for part in parts:
                            if part.startswith('core_len='):
                                try:
                                    core_len = int(part.split('=')[1])
                                except (ValueError, IndexError):
                                    pass
                            elif part.startswith('var_len='):
                                try:
                                    var_len = int(part.split('=')[1])
                                except (ValueError, IndexError):
                                    pass
                            elif '(' in part and ')' in part:
                                # Extract strand from coordinates like chr1:213-962(-)
                                strand_match = part[part.find('(')+1:part.find(')')]
                                if strand_match in ['+', '-']:
                                    strand = strand_match

                        boundaries[chr_end] = {
                            'core_len': core_len,
                            'var_len': var_len,
                            'strand': strand
                        }
    except FileNotFoundError:
        pass

    return boundaries


def split_xprime_to_core_variable(xprime_hit: Dict, boundaries: Dict[str, Dict]) -> Tuple[Optional[Dict], Optional[Dict]]:
    """
    Split a detected X prime hit into core and variable elements based on reference proportions.

    The split is based on the ratio of core_len to var_len from the source X prime.
    For L arms (minus strand): variable is at lower coordinates, core at higher
    For R arms (plus strand): core is at lower coordinates, variable at higher

    Args:
        xprime_hit: X prime hit dictionary with 'start', 'end', 'source_chr_end', 'strand'
        boundaries: Dict mapping chr_end to core/variable lengths

    Returns:
        Tuple of (core_element_dict, variable_element_dict) - either can be None
    """
    source_chr_end = xprime_hit.get('source_chr_end')
    if not source_chr_end or source_chr_end not in boundaries:
        # No boundary info available, return None
        return None, None

    boundary_info = boundaries[source_chr_end]
    core_len = boundary_info['core_len']
    var_len = boundary_info['var_len']
    ref_strand = boundary_info['strand']

    # If no variable element, return only core
    if var_len == 0:
        core_element = xprime_hit.copy()
        core_element['type'] = 'x_core_element'
        return core_element, None

    # Calculate the total reference length and proportions
    total_ref_len = core_len + var_len
    if total_ref_len == 0:
        return None, None

    core_ratio = core_len / total_ref_len
    var_ratio = var_len / total_ref_len

    # Get the detected X prime coordinates
    detected_start = xprime_hit['start']
    detected_end = xprime_hit['end']
    detected_len = detected_end - detected_start

    # Calculate split lengths based on ratios
    split_core_len = int(detected_len * core_ratio)
    split_var_len = detected_len - split_core_len  # Remainder goes to variable

    # Determine arm from target chr_end
    target_chr_end = xprime_hit.get('target_chr_end', source_chr_end)
    arm = target_chr_end[-1] if target_chr_end else 'L'

    # Determine order based on arm
    # L arm (minus strand orientation): telomere at low coords
    #   Order from low to high: Telomere - [Y'] - Variable - Core - Space - Anchor
    # R arm (plus strand orientation): telomere at high coords
    #   Order from low to high: Anchor - Space - Core - Variable - [Y'] - Telomere

    if arm == 'L':
        # L arm: variable at lower coordinates, core at higher
        var_start = detected_start
        var_end = detected_start + split_var_len
        core_start = var_end
        core_end = detected_end
        strand = '-'
    else:
        # R arm: core at lower coordinates, variable at higher
        core_start = detected_start
        core_end = detected_start + split_core_len
        var_start = core_end
        var_end = detected_end
        strand = '+'

    # Create core element
    core_element = xprime_hit.copy()
    core_element['start'] = core_start
    core_element['end'] = core_end
    core_element['type'] = 'x_core_element'
    core_element['strand'] = strand
    core_element['length'] = core_end - core_start

    # Create variable element
    var_element = xprime_hit.copy()
    var_element['start'] = var_start
    var_element['end'] = var_end
    var_element['type'] = 'x_variable_element'
    var_element['strand'] = strand
    var_element['length'] = var_end - var_start

    return core_element, var_element


def extract_chr_end_from_subject(sseqid: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract chromosome and determine arm from subject chromosome name.

    Args:
        sseqid: Subject sequence ID (e.g., "chr1_extended", "chr10")

    Returns:
        Tuple of (chromosome_base, None) - arm determined by position later
    """
    # Extract chromosome number
    sseqid_lower = sseqid.lower()

    # Handle formats like chr1_extended, chr1, chromosome_1, etc.
    if 'chr' in sseqid_lower:
        # Extract the number after 'chr'
        parts = sseqid_lower.replace('chr', '').split('_')[0]
        if parts.isdigit():
            return f"chr{parts}", None

    return None, None


def determine_arm_from_position(position: int, chr_length: int) -> str:
    """
    Determine chromosome arm (L or R) based on position.

    L arm = position < midpoint (closer to start)
    R arm = position >= midpoint (closer to end)

    Args:
        position: Position on chromosome
        chr_length: Total chromosome length

    Returns:
        'L' or 'R'
    """
    midpoint = chr_length / 2
    return 'L' if position < midpoint else 'R'


def filter_xprime_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter X prime BLAST hits for quality and calculate adjusted identity.

    Args:
        df: DataFrame with BLAST results

    Returns:
        Filtered DataFrame with adjusted_pident column
    """
    if len(df) == 0:
        return df

    # Calculate adjusted percent identity
    df = df.copy()
    df['adjusted_pident'] = df.apply(
        lambda row: calculate_adjusted_pident(
            row['pident'],
            row['length'],
            row['qlen']
        ),
        axis=1
    )

    # Filter by adjusted identity
    df = df[df['adjusted_pident'] >= MIN_XPRIME_ADJUSTED_IDENTITY]

    # Filter by standard identity
    df = df[df['pident'] >= MIN_XPRIME_IDENTITY]

    # Calculate query coverage
    df['query_coverage'] = (df['length'] / df['qlen']) * 100

    # Filter by coverage
    df = df[df['query_coverage'] >= MIN_XPRIME_COVERAGE]

    return df


def assign_xprimes_to_chr_ends(xprime_df: pd.DataFrame, chr_end_regions: Dict,
                                xprime_fasta_path: Optional[str] = None) -> Tuple[Dict, List[Dict], Dict]:
    """
    Assign X prime hits to chromosome ends, selecting the best match for each end.

    For each chromosome end in the reference:
    1. Find all X prime hits on that chromosome arm
    2. Select the best hit (highest adjusted identity)
    3. Check if the best hit's source matches the expected chr end
    4. Warn if there's a mismatch or low quality
    5. Split into core and variable elements if boundary info available

    Args:
        xprime_df: DataFrame with filtered X prime BLAST results
        chr_end_regions: Dict of chromosome end regions (from anchor detection)
        xprime_fasta_path: Optional path to X prime FASTA file for core/variable boundary info

    Returns:
        Tuple of:
        - Dict mapping chr_end to list of X prime hits (combined)
        - List of quality report entries
        - Dict mapping chr_end to {'core': core_element, 'variable': var_element}
    """
    chr_end_xprimes = defaultdict(list)
    chr_end_xprime_split = {}  # Store split core/variable elements
    quality_report = []

    # Parse boundary info from FASTA if provided
    boundaries = {}
    if xprime_fasta_path:
        boundaries = parse_xprime_fasta_for_boundaries(xprime_fasta_path)

    if len(xprime_df) == 0:
        return chr_end_xprimes, quality_report, chr_end_xprime_split

    # Get chromosome sizes from the data
    chr_sizes = xprime_df.groupby('sseqid')['slen'].first().to_dict()

    # Process hits: determine which chr end each hit belongs to
    hits_by_chr_end = defaultdict(list)

    for _, hit in xprime_df.iterrows():
        subject_chr = hit['sseqid']
        chr_base, _ = extract_chr_end_from_subject(subject_chr)

        if chr_base is None:
            continue

        # Get position (midpoint of hit)
        hit_start = min(hit['sstart'], hit['send'])
        hit_end = max(hit['sstart'], hit['send'])
        hit_midpoint = (hit_start + hit_end) / 2

        # Determine arm based on position
        chr_length = chr_sizes.get(subject_chr, 0)
        if chr_length == 0:
            continue

        arm = determine_arm_from_position(hit_midpoint, chr_length)
        target_chr_end = f"{chr_base}{arm}"

        # Get the source chr end from the query
        source_chr_end = extract_chr_end_from_xprime_query(hit['qseqid'])

        hit_info = {
            'chr': subject_chr,
            'start': hit_start,
            'end': hit_end,
            'strand': '+' if hit['sstart'] < hit['send'] else '-',
            'length': hit_end - hit_start,
            'pident': hit['pident'],
            'adjusted_pident': hit['adjusted_pident'],
            'query_coverage': hit['query_coverage'],
            'bitscore': hit['bitscore'],
            'source': hit['qseqid'],
            'source_chr_end': source_chr_end,
            'target_chr_end': target_chr_end,
            'qlen': hit['qlen']
        }

        hits_by_chr_end[target_chr_end].append(hit_info)

    # For each chr end, select the best X prime hit
    for chr_end in sorted(hits_by_chr_end.keys()):
        hits = hits_by_chr_end[chr_end]

        if not hits:
            continue

        # Sort by adjusted identity (descending), then bitscore
        hits_sorted = sorted(hits, key=lambda x: (x['adjusted_pident'], x['bitscore']), reverse=True)
        best_hit = hits_sorted[0]

        # Check for mismatches and quality issues
        source_chr_end = best_hit['source_chr_end']

        # Check 1: Does the best match come from the expected chr end?
        if source_chr_end and source_chr_end != chr_end:
            quality_report.append({
                'type': 'xprime_mismatch',
                'chr_end': chr_end,
                'expected_source': chr_end,
                'actual_source': source_chr_end,
                'adjusted_pident': best_hit['adjusted_pident'],
                'pident': best_hit['pident'],
                'message': f"X prime at {chr_end} best matches {source_chr_end}_xprime (not {chr_end}_xprime)"
            })

        # Check 2: Is the match high quality?
        if best_hit['adjusted_pident'] < MIN_HIGH_QUALITY_IDENTITY:
            quality_report.append({
                'type': 'xprime_low_quality',
                'chr_end': chr_end,
                'source': best_hit['source'],
                'adjusted_pident': best_hit['adjusted_pident'],
                'pident': best_hit['pident'],
                'query_coverage': best_hit['query_coverage'],
                'message': f"X prime at {chr_end} has low adjusted identity ({best_hit['adjusted_pident']:.1f}% < {MIN_HIGH_QUALITY_IDENTITY}%)"
            })

        # Add type field for output
        best_hit['type'] = 'x_prime'
        chr_end_xprimes[chr_end].append(best_hit)

        # Split into core and variable elements if boundary info available
        if boundaries:
            core_element, var_element = split_xprime_to_core_variable(best_hit, boundaries)
            chr_end_xprime_split[chr_end] = {
                'core': core_element,
                'variable': var_element
            }

    return chr_end_xprimes, quality_report, chr_end_xprime_split


def validate_xprime_positions(chr_end_xprimes: Dict, chr_end_yprimes: Dict, chr_sizes: Dict) -> List[Dict]:
    """
    Validate that X primes are positioned correctly relative to Y primes.

    X prime should always be closer to the chromosome middle than the first Y prime.

    Args:
        chr_end_xprimes: Dict mapping chr_end to X prime hits
        chr_end_yprimes: Dict mapping chr_end to Y prime hits
        chr_sizes: Dict mapping chromosome name to size

    Returns:
        List of position warnings
    """
    position_warnings = []

    for chr_end, xprimes in chr_end_xprimes.items():
        if not xprimes:
            continue

        # Get Y primes for this chr end
        yprimes = chr_end_yprimes.get(chr_end, [])
        if not yprimes:
            # No Y primes to compare against - that's OK
            continue

        xprime = xprimes[0]  # Best X prime
        xprime_chr = xprime['chr']
        xprime_midpoint = (xprime['start'] + xprime['end']) / 2

        # Get chromosome size
        chr_size = chr_sizes.get(xprime_chr, 0)
        if chr_size == 0:
            continue

        chr_midpoint = chr_size / 2

        # Determine which arm we're on
        arm = chr_end[-1]  # L or R

        # Find the Y prime closest to chromosome middle
        closest_yprime = None
        closest_yprime_dist = float('inf')

        for yp in yprimes:
            yp_midpoint = (yp['start'] + yp['end']) / 2
            dist_to_middle = abs(yp_midpoint - chr_midpoint)
            if dist_to_middle < closest_yprime_dist:
                closest_yprime_dist = dist_to_middle
                closest_yprime = yp

        if closest_yprime is None:
            continue

        yprime_midpoint = (closest_yprime['start'] + closest_yprime['end']) / 2

        # Calculate distances to chromosome middle
        xprime_dist_to_middle = abs(xprime_midpoint - chr_midpoint)
        yprime_dist_to_middle = abs(yprime_midpoint - chr_midpoint)

        # X prime should be closer to middle than Y prime
        if xprime_dist_to_middle >= yprime_dist_to_middle:
            position_warnings.append({
                'type': 'xprime_position_error',
                'chr_end': chr_end,
                'xprime_pos': f"{xprime['start']}-{xprime['end']}",
                'yprime_pos': f"{closest_yprime['start']}-{closest_yprime['end']}",
                'xprime_dist_to_middle': int(xprime_dist_to_middle),
                'yprime_dist_to_middle': int(yprime_dist_to_middle),
                'message': f"X prime at {chr_end} is NOT closer to chr middle than Y prime! "
                          f"X prime dist: {int(xprime_dist_to_middle):,}bp, Y prime dist: {int(yprime_dist_to_middle):,}bp"
            })

    return position_warnings


def get_chr_sizes_from_blast(df: pd.DataFrame) -> Dict[str, int]:
    """Extract chromosome sizes from BLAST results."""
    if len(df) == 0:
        return {}
    return df.groupby('sseqid')['slen'].first().to_dict()
