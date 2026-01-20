#!/usr/bin/env python3
"""
Label pre-telomeric regions (anchor, Y prime, X prime) in a reference genome.

This script uses BLAST to align known anchor and Y prime sequences from a reference strain
to a new strain's reference genome, then annotates the regions.

Includes Y prime probe verification to ensure detected Y prime counts match expected counts.
"""

import argparse
import os
import subprocess
import sys
import pandas as pd
from collections import defaultdict

# Import Y prime detection module (v2 with better fragment filtering)
from yprime_detection_v2 import (
    assign_yprimes_to_chr_ends,
    MIN_YPRIME_IDENTITY,
    MIN_YPRIME_COVERAGE
)

# Import X prime detection module
from xprime_detection import (
    assign_xprimes_to_chr_ends,
    filter_xprime_hits,
    validate_xprime_positions,
    get_chr_sizes_from_blast,
    MIN_XPRIME_IDENTITY,
    MIN_HIGH_QUALITY_IDENTITY
)


def run_probe_blast(probe_fasta, reference_fasta, output_file, min_identity=90, num_threads=4):
    """
    BLAST the Y prime probe against the reference to count expected Y primes.

    Args:
        probe_fasta: Path to probe sequence FASTA
        reference_fasta: Path to reference genome FASTA
        output_file: Output file for BLAST results
        min_identity: Minimum percent identity (default 90%)
        num_threads: Number of threads for BLAST

    Returns:
        Path to output file
    """
    blast_cmd = [
        'blastn',
        '-query', probe_fasta,
        '-subject', reference_fasta,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
        '-perc_identity', str(min_identity),
        '-num_threads', str(num_threads),
        '-task', 'blastn'
    ]

    subprocess.run(blast_cmd, check=True, capture_output=True)
    return output_file


def count_yprimes_from_probe(probe_blast_file, reference_fasta, min_length=80):
    """
    Count expected Y primes per chromosome end using probe BLAST results.

    Args:
        probe_blast_file: Path to probe BLAST output
        reference_fasta: Path to reference FASTA (for chromosome sizes)
        min_length: Minimum alignment length for probe hits

    Returns:
        Dict mapping chr_end to expected Y prime count
    """
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen'
    ]

    if not os.path.exists(probe_blast_file) or os.path.getsize(probe_blast_file) == 0:
        return {}

    df = pd.read_csv(probe_blast_file, sep='\t', names=columns)

    # Filter for high-quality matches
    df = df[df['length'] >= min_length]

    if len(df) == 0:
        return {}

    # Calculate position on chromosome (use midpoint of hit)
    df['hit_position'] = (df['sstart'] + df['send']) / 2

    # Get chromosome sizes from slen column
    chr_sizes = df.groupby('sseqid')['slen'].first().to_dict()

    # Count Y primes per chromosome end
    yprime_counts = defaultdict(int)

    for _, hit in df.iterrows():
        chr_name = hit['sseqid']
        position = hit['hit_position']
        chr_size = chr_sizes.get(chr_name, 0)

        if chr_size == 0:
            continue

        # Extract chromosome number from name (e.g., chr1_extended -> 1)
        chr_num = None
        for part in chr_name.lower().replace('chr', '').replace('_extended', '').split('_'):
            if part.isdigit():
                chr_num = part
                break

        if chr_num is None:
            continue

        # Determine which arm based on position
        midpoint = chr_size / 2
        arm = 'L' if position < midpoint else 'R'

        chr_end = f"chr{chr_num}{arm}"
        yprime_counts[chr_end] += 1

    return dict(yprime_counts)


def verify_yprime_counts(detected_counts, expected_counts):
    """
    Verify that detected Y prime counts match expected counts from probe.

    Args:
        detected_counts: Dict mapping chr_end to detected Y prime count
        expected_counts: Dict mapping chr_end to expected Y prime count (from probe)

    Returns:
        Tuple: (matches: bool, mismatches: list of dicts with details)
    """
    mismatches = []

    # Get all chromosome ends from both sources
    all_chr_ends = set(detected_counts.keys()) | set(expected_counts.keys())

    for chr_end in sorted(all_chr_ends):
        detected = detected_counts.get(chr_end, 0)
        expected = expected_counts.get(chr_end, 0)

        if detected != expected:
            mismatches.append({
                'chr_end': chr_end,
                'detected': detected,
                'expected': expected,
                'difference': detected - expected
            })

    return len(mismatches) == 0, mismatches


def run_blast(query_fasta, subject_fasta, output_file, evalue=1e-10, num_threads=4):
    """
    Run BLAST to find similar regions between query and subject sequences.

    Args:
        query_fasta: Path to query FASTA file
        subject_fasta: Path to subject FASTA file
        output_file: Path to output BLAST results
        evalue: E-value threshold
        num_threads: Number of threads to use

    Returns:
        Path to output file
    """
    blast_cmd = [
        'blastn',
        '-query', query_fasta,
        '-subject', subject_fasta,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
        '-evalue', str(evalue),
        '-num_threads', str(num_threads),
        '-task', 'blastn'
    ]

    print(f"Running BLAST: {' '.join(blast_cmd)}")
    subprocess.run(blast_cmd, check=True)

    return output_file


def parse_blast_results(blast_file, min_pident=80, min_length=100):
    """
    Parse BLAST output and filter for high-quality alignments.

    Args:
        blast_file: Path to BLAST output file
        min_pident: Minimum percent identity
        min_length: Minimum alignment length

    Returns:
        pd.DataFrame with filtered BLAST results
    """
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen'
    ]

    if os.path.getsize(blast_file) == 0:
        print(f"Warning: BLAST file {blast_file} is empty")
        return pd.DataFrame(columns=columns)

    df = pd.read_csv(blast_file, sep='\t', names=columns)

    # Filter by quality
    df = df[df['pident'] >= min_pident]
    df = df[df['length'] >= min_length]

    # Calculate coverage
    df['query_coverage'] = (df['length'] / df['qlen']) * 100
    df['subject_coverage'] = (df['length'] / df['slen']) * 100

    # Ensure proper strand orientation (start < end)
    df['strand'] = df.apply(lambda row: '+' if row['sstart'] < row['send'] else '-', axis=1)
    df['start'] = df[['sstart', 'send']].min(axis=1)
    df['end'] = df[['sstart', 'send']].max(axis=1)

    return df


def extract_chr_end_from_query(qseqid):
    """
    Extract chromosome end information from query sequence ID.

    Examples:
        chr1L_anchor -> (chr1, L, anchor)
        Y_Prime_chr2L1#Short/Solo/ID4_Green-Light -> (chr2, L, y_prime)

    Args:
        qseqid: Query sequence identifier

    Returns:
        tuple: (chromosome, arm, region_type) or (None, None, None)
    """
    if '_anchor' in qseqid:
        # Format: chr1L_anchor
        parts = qseqid.replace('_anchor', '')
        if len(parts) >= 4:
            chrom = parts[:-1]  # Everything except last character
            arm = parts[-1]     # Last character (L or R)
            return (chrom, arm, 'anchor')
    elif qseqid.startswith('Y_Prime_'):
        # Format: Y_Prime_chr2L1#...
        parts = qseqid.split('#')[0]  # Remove everything after #
        parts = parts.replace('Y_Prime_', '')

        # Extract chr and arm from patterns like chr2L1 or chr4R1,2,3,6,7
        chr_parts = parts.split(';')[0].split(',')[0]  # Take first chr if multiple

        # Find where the arm letter is (look for L or R)
        for i, char in enumerate(chr_parts):
            if char in ['L', 'R']:
                chrom = chr_parts[:i+1][:-1]  # Everything before L/R
                arm = char
                return (chrom, arm, 'y_prime')

    return (None, None, None)


def check_chromosome_match(expected_chr, subject_chr):
    """
    Check if subject chromosome matches expected chromosome.

    Handles naming variations like:
    - chr1 matches chr1_extended
    - chr10 matches chr10_extended

    Args:
        expected_chr: Expected chromosome (e.g., "chr1")
        subject_chr: Actual subject chromosome (e.g., "chr1_extended")

    Returns:
        bool: True if chromosomes match
    """
    # Normalize both to lowercase for comparison
    expected_lower = expected_chr.lower()
    subject_lower = subject_chr.lower()

    # Check if subject starts with expected chromosome
    # This handles chr1_extended, chr1_7575, etc.
    if subject_lower.startswith(expected_lower):
        # Make sure we don't match chr1 with chr10
        next_char_idx = len(expected_lower)
        if next_char_idx >= len(subject_lower):
            return True  # Exact match
        next_char = subject_lower[next_char_idx]
        # Next character should be non-digit (like '_' in chr1_extended)
        return not next_char.isdigit()

    return False


def assign_regions_to_chr_ends(anchor_df, yprime_df, enforce_synteny=True, min_high_quality_pident=95.0):
    """
    Assign detected anchor and Y prime regions to specific chromosome ends.

    With strict synteny enforcement:
    - Anchors must match their expected chromosome
    - Warnings issued for low-quality matches (<95% identity)
    - Mismatches are reported but not included

    Args:
        anchor_df: DataFrame with anchor BLAST results
        yprime_df: DataFrame with Y prime BLAST results
        enforce_synteny: If True, only accept matches to expected chromosome
        min_high_quality_pident: Minimum % identity for high-quality match (default: 95.0)

    Returns:
        tuple: (chr_end_regions dict, quality_report dict)
    """
    chr_end_regions = defaultdict(lambda: {'anchor': [], 'y_prime': [], 'x_prime': []})

    # Quality tracking
    quality_report = {
        'anchor_mismatches': [],
        'anchor_low_quality': [],
        'anchor_missing': [],
        'yprime_mismatches': [],
        'yprime_low_quality': []
    }

    # Process anchors with strict synteny checking
    for _, row in anchor_df.iterrows():
        chrom, arm, region_type = extract_chr_end_from_query(row['qseqid'])
        if chrom and arm:
            chr_end = f"{chrom}{arm}"
            subject_chr = row['sseqid']
            pident = row['pident']

            # Check if chromosome matches
            chr_matches = check_chromosome_match(chrom, subject_chr)

            if enforce_synteny and not chr_matches:
                # Mismatch - anchor mapped to wrong chromosome
                quality_report['anchor_mismatches'].append({
                    'query': row['qseqid'],
                    'expected_chr': chrom,
                    'found_chr': subject_chr,
                    'chr_end': chr_end,
                    'arm': arm,
                    'pident': pident,
                    'bitscore': row['bitscore'],
                    'query_coverage': row['query_coverage']
                })
                continue  # Skip this hit

            # Check quality thresholds (BOTH identity AND coverage)
            quality_issues = []

            # Check 1: Percent identity (should be ≥95%)
            if pident < min_high_quality_pident:
                quality_issues.append(f"low identity ({pident:.1f}% < {min_high_quality_pident}%)")

            # Check 2: Query coverage (should be ≥95% - complete anchor, not fragment)
            min_coverage = 95.0
            if row['query_coverage'] < min_coverage:
                quality_issues.append(f"incomplete anchor ({row['query_coverage']:.1f}% coverage < {min_coverage}%)")

            # If any quality issues, add to report
            if quality_issues:
                quality_report['anchor_low_quality'].append({
                    'query': row['qseqid'],
                    'chr': subject_chr,
                    'chr_end': chr_end,
                    'arm': arm,
                    'pident': pident,
                    'bitscore': row['bitscore'],
                    'query_coverage': row['query_coverage'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'issues': ', '.join(quality_issues)
                })

            # Add to results (even if low quality - will be filtered later if needed)
            chr_end_regions[chr_end]['anchor'].append({
                'chr': subject_chr,
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': row['strand'],
                'type': 'anchor',
                'source': row['qseqid'],
                'pident': row['pident'],
                'bitscore': row['bitscore'],
                'query_coverage': row['query_coverage']
            })

    # NOTE: Y prime processing is now handled separately with fragment merging
    # See assign_yprimes_to_chr_ends() in yprime_detection module
    # This allows for:
    # - No chromosome constraint (any Y prime can be on any chr end)
    # - Fragment merging (combine small fragments into complete Y primes)
    # - Tandem repeat detection (multiple Y primes close together)
    # - Proper size classification (short vs long Y primes)

    # CRITICAL FIX: Keep only the BEST anchor hit per chromosome end
    # BLAST often finds multiple partial/spurious matches
    # We want only the highest-scoring (best) match for each chr_end
    print("\n    Filtering to keep only BEST anchor per chromosome end...")
    for chr_end in list(chr_end_regions.keys()):
        anchors = chr_end_regions[chr_end]['anchor']
        if len(anchors) > 1:
            # Multiple anchors found - keep only the best (highest bitscore)
            best_anchor = max(anchors, key=lambda x: x['bitscore'])
            removed_count = len(anchors) - 1
            print(f"      {chr_end}: {len(anchors)} hits found → keeping best (bitscore: {best_anchor['bitscore']:.1f})")

            # Replace with single best anchor
            chr_end_regions[chr_end]['anchor'] = [best_anchor]

            # Add removed hits to quality report for transparency
            for anchor in anchors:
                if anchor != best_anchor:
                    quality_report['anchor_low_quality'].append({
                        'query': anchor['source'],
                        'chr': anchor['chr'],
                        'chr_end': chr_end,
                        'arm': chr_end[-1],
                        'pident': anchor['pident'],
                        'bitscore': anchor['bitscore'],
                        'query_coverage': anchor['query_coverage'],
                        'start': anchor['start'],
                        'end': anchor['end'],
                        'reason': f'Secondary hit (removed, kept best with bitscore {best_anchor["bitscore"]:.1f})'
                    })

    # Check for missing anchors (expected chr ends with no anchor found)
    # Expected: chr1L, chr1R, chr2L, chr2R, ..., chr16L, chr16R (32 total)
    expected_chr_ends = []
    for i in range(1, 17):  # chr1 to chr16
        expected_chr_ends.append(f"chr{i}L")
        expected_chr_ends.append(f"chr{i}R")

    found_chr_ends_with_anchors = set(chr_end for chr_end, regions in chr_end_regions.items() if regions['anchor'])
    missing_anchors = set(expected_chr_ends) - found_chr_ends_with_anchors
    quality_report['anchor_missing'] = sorted(missing_anchors)

    return chr_end_regions, quality_report


def infer_x_prime_regions(chr_end_regions):
    """
    Infer X prime regions based on anchor and Y prime positions.

    X prime is typically located between anchor and Y prime.
    If no Y prime is found, X prime extends from anchor to the expected telomere position.

    Args:
        chr_end_regions: Dict of chromosome end regions

    Returns:
        Updated chr_end_regions with X prime annotations
    """
    for chr_end, regions in chr_end_regions.items():
        anchors = regions['anchor']
        yprimes = regions['y_prime']

        if not anchors:
            continue

        # Get the best anchor (highest bitscore)
        best_anchor = max(anchors, key=lambda x: x['bitscore'])

        if yprimes:
            # Get the best Y prime
            best_yprime = max(yprimes, key=lambda x: x['bitscore'])

            # X prime is between anchor and Y prime
            # Determine order based on position
            if best_anchor['end'] < best_yprime['start']:
                # Anchor comes before Y prime
                x_start = best_anchor['end'] + 1
                x_end = best_yprime['start'] - 1
            else:
                # Y prime comes before anchor (unusual but possible)
                x_start = best_yprime['end'] + 1
                x_end = best_anchor['start'] - 1

            if x_start < x_end and (x_end - x_start) > 100:  # Minimum 100bp for X prime
                regions['x_prime'].append({
                    'chr': best_anchor['chr'],
                    'start': x_start,
                    'end': x_end,
                    'strand': best_anchor['strand'],
                    'type': 'x_prime',
                    'source': 'inferred_between_anchor_and_yprime',
                    'pident': None,
                    'bitscore': None,
                    'query_coverage': None
                })
        else:
            # No Y prime found - X prime extends from anchor toward telomere
            # Assume X prime is ~1-5kb from anchor end
            x_start = best_anchor['end'] + 1
            x_end = x_start + 3000  # Default 3kb X prime region

            regions['x_prime'].append({
                'chr': best_anchor['chr'],
                'start': x_start,
                'end': x_end,
                'strand': best_anchor['strand'],
                'type': 'x_prime',
                'source': 'inferred_no_yprime',
                'pident': None,
                'bitscore': None,
                'query_coverage': None
            })

    return chr_end_regions


def write_gff3(chr_end_regions, output_file, source='TeloTracker'):
    """
    Write annotations in GFF3 format.

    Args:
        chr_end_regions: Dict of chromosome end regions
        output_file: Path to output GFF3 file
        source: Source field for GFF3
    """
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")

        feature_id = 1
        for chr_end, regions in sorted(chr_end_regions.items()):
            for region_type in ['anchor', 'x_prime', 'y_prime']:
                for region in regions[region_type]:
                    attributes = [
                        f"ID={region_type}_{feature_id}",
                        f"Name={chr_end}_{region_type}",
                        f"chr_end={chr_end}"
                    ]

                    if region['source']:
                        attributes.append(f"source_seq={region['source']}")
                    if region['pident'] is not None:
                        attributes.append(f"pident={region['pident']:.2f}")
                    if region['bitscore'] is not None:
                        attributes.append(f"bitscore={region['bitscore']:.2f}")
                    if region['query_coverage'] is not None:
                        attributes.append(f"query_coverage={region['query_coverage']:.2f}")

                    attributes_str = ';'.join(attributes)

                    f.write(f"{region['chr']}\t{source}\t{region_type}\t"
                           f"{region['start']}\t{region['end']}\t.\t"
                           f"{region['strand']}\t.\t{attributes_str}\n")

                    feature_id += 1


def write_bed(chr_end_regions, output_file):
    """
    Write annotations in BED format.

    Args:
        chr_end_regions: Dict of chromosome end regions
        output_file: Path to output BED file
    """
    with open(output_file, 'w') as f:
        for chr_end, regions in sorted(chr_end_regions.items()):
            for region_type in ['anchor', 'x_prime', 'y_prime']:
                for region in regions[region_type]:
                    # BED format: chr start end name score strand
                    name = f"{chr_end}_{region_type}"
                    score = int(region['bitscore']) if region['bitscore'] is not None else 0

                    f.write(f"{region['chr']}\t{region['start']}\t{region['end']}\t"
                           f"{name}\t{score}\t{region['strand']}\n")


def write_bed_simplified(chr_end_regions, output_file):
    """
    Write annotations in simplified BED format (matching 6991_final_features.bed format).

    Format: chr start end name strand length
    (No score column, length added at end)

    Includes space_between_anchor entries showing the gap between anchor and X prime.

    Args:
        chr_end_regions: Dict of chromosome end regions
        output_file: Path to output BED file
    """
    def clean_chr_name(chr_name):
        """Remove _extended suffix from chromosome names."""
        if chr_name.endswith('_extended'):
            return chr_name[:-9]  # Remove '_extended' (9 characters)
        return chr_name

    with open(output_file, 'w') as f:
        for chr_end, regions in sorted(chr_end_regions.items()):
            # Get anchor and x_prime positions for calculating space_between_anchor
            anchors = regions.get('anchor', [])
            xprimes = regions.get('x_prime', [])

            anchor_info = None
            xprime_info = None

            if anchors:
                anchor_info = anchors[0]
            if xprimes:
                xprime_info = xprimes[0]

            for region_type in ['anchor', 'x_prime', 'y_prime']:
                for i, region in enumerate(regions[region_type]):
                    # Simplified BED format: chr start end name strand length
                    # Name format matches 6991_final_features.bed style
                    if region_type == 'y_prime':
                        # Y primes are numbered: chr1L_Y_Prime_1, chr1L_Y_Prime_2, etc.
                        name = f"{chr_end}_Y_Prime_{i+1}"
                    elif region_type == 'x_prime':
                        name = f"{chr_end}_x_prime"
                    else:
                        name = f"{chr_end}_{region_type}"

                    length = region['end'] - region['start']
                    chr_clean = clean_chr_name(region['chr'])

                    f.write(f"{chr_clean}\t{region['start']}\t{region['end']}\t"
                           f"{name}\t{region['strand']}\t{length}\n")

            # Calculate and write space_between_anchor if both anchor and x_prime exist
            if anchor_info and xprime_info:
                chr_name = clean_chr_name(anchor_info['chr'])
                anchor_start = anchor_info['start']
                anchor_end = anchor_info['end']
                xprime_start = xprime_info['start']
                xprime_end = xprime_info['end']

                # Determine arm (L or R) from chr_end name
                arm = chr_end[-1]  # Last character is L or R

                if arm == 'L':
                    # L arm: telomere at low coordinates, anchor closer to middle
                    # Structure: TELO---[Y']---[X']---[space]---[ANCHOR]---middle
                    # Space is between X prime end and anchor start
                    space_start = xprime_end + 1
                    space_end = anchor_start - 1
                    strand = '+'
                else:
                    # R arm: telomere at high coordinates, anchor closer to middle
                    # Structure: middle---[ANCHOR]---[space]---[X']---[Y']---TELO
                    # Space is between anchor end and X prime start
                    space_start = anchor_end + 1
                    space_end = xprime_start - 1
                    strand = '+'

                # Only write if there's actually a gap (space_end > space_start)
                if space_end >= space_start:
                    space_length = space_end - space_start + 1
                    space_name = f"{chr_end}_space_between_anchor"
                    f.write(f"{chr_name}\t{space_start}\t{space_end}\t"
                           f"{space_name}\t{strand}\t{space_length}\n")


def write_tsv(chr_end_regions, output_file):
    """
    Write annotations in TSV format for downstream analysis.

    Args:
        chr_end_regions: Dict of chromosome end regions
        output_file: Path to output TSV file
    """
    rows = []

    for chr_end, regions in sorted(chr_end_regions.items()):
        for region_type in ['anchor', 'x_prime', 'y_prime']:
            for region in regions[region_type]:
                rows.append({
                    'chr_end': chr_end,
                    'chr': region['chr'],
                    'start': region['start'],
                    'end': region['end'],
                    'length': region['end'] - region['start'] + 1,
                    'strand': region['strand'],
                    'region_type': region_type,
                    'source': region['source'],
                    'pident': region['pident'],
                    'bitscore': region['bitscore'],
                    'query_coverage': region['query_coverage']
                })

    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)


def print_quality_report(chr_end_regions, quality_report, yprime_quality=None, min_high_quality_pident=95.0):
    """
    Print a concise quality report focusing only on regions that are USED in final output.

    Args:
        chr_end_regions: Dict of chromosome end regions (the actual output)
        quality_report: Dict with quality tracking information
        yprime_quality: List with Y prime quality information (optional)
        min_high_quality_pident: Minimum % identity for high-quality anchors
    """
    print("\n" + "=" * 80)
    print("QUALITY REPORT - USED REGIONS ONLY")
    print("=" * 80)

    has_issues = False

    # Check anchors that are actually USED for quality issues
    low_quality_used_anchors = []
    for chr_end, regions in chr_end_regions.items():
        for anchor in regions.get('anchor', []):
            pident = anchor.get('pident', 0)
            coverage = anchor.get('query_coverage', 0)

            issues = []
            if pident < min_high_quality_pident:
                issues.append(f"low identity ({pident:.1f}% < {min_high_quality_pident}%)")
            if coverage < 95.0:
                issues.append(f"incomplete ({coverage:.1f}% coverage)")

            if issues:
                low_quality_used_anchors.append({
                    'chr_end': chr_end,
                    'chr': anchor['chr'],
                    'start': anchor['start'],
                    'end': anchor['end'],
                    'pident': pident,
                    'coverage': coverage,
                    'bitscore': anchor.get('bitscore', 0),
                    'issues': ', '.join(issues)
                })

    # Report low-quality anchors that ARE USED
    if low_quality_used_anchors:
        has_issues = True
        print(f"\n⚠️  WARNING: {len(low_quality_used_anchors)} LOW QUALITY ANCHOR(S) IN USE")
        print("-" * 80)
        for a in low_quality_used_anchors:
            print(f"  {a['chr_end']:8s} | {a['chr']:15s} | {a['start']:>10,}-{a['end']:>10,}")
            print(f"           Identity: {a['pident']:.1f}% | Coverage: {a['coverage']:.1f}%")
            print(f"           Issues: {a['issues']}")
        print()

    # Report missing anchors (WARNING)
    if quality_report['anchor_missing']:
        has_issues = True
        print(f"\n⚠️  WARNING: {len(quality_report['anchor_missing'])} MISSING ANCHOR(S)")
        print("-" * 80)
        missing = quality_report['anchor_missing']
        # Print in rows of 8
        for i in range(0, len(missing), 8):
            row = missing[i:i+8]
            print("  " + "  ".join(f"{x:6s}" for x in row))
        print()

    # Report Y prime summary (concise)
    if yprime_quality:
        # Count by type
        merged_count = sum(1 for yp in yprime_quality if yp.get('merged', False))
        fragment_count = sum(1 for yp in yprime_quality if yp.get('class') == 'fragment')
        short_count = sum(1 for yp in yprime_quality if yp.get('class') == 'short')
        long_count = sum(1 for yp in yprime_quality if yp.get('class') == 'long')

        print(f"\nℹ️  Y PRIME SUMMARY: {len(yprime_quality)} total")
        print(f"   Short: {short_count} | Long: {long_count} | Fragments: {fragment_count} | Merged: {merged_count}")

        # Only show warning if there are fragments
        if fragment_count > 0:
            has_issues = True
            print(f"\n   ⚠️  {fragment_count} fragment(s) could not form complete Y primes")

    if not has_issues:
        print("\n✅ All used anchors and Y primes have high quality!")

    print("\n" + "=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description='Label pre-telomeric regions (anchor, Y prime, X prime) in a reference genome'
    )

    parser.add_argument(
        '--reference',
        required=True,
        help='Reference FASTA file (output from create_ref.sh)'
    )
    parser.add_argument(
        '--anchors',
        required=True,
        help='Anchor sequences FASTA file (e.g., test_anchors.fasta)'
    )
    parser.add_argument(
        '--yprimes',
        required=True,
        help='Y prime sequences FASTA file (e.g., repeatmasker_6991_all_y_primes.fasta)'
    )
    parser.add_argument(
        '--xprimes',
        default=None,
        help='X prime sequences FASTA file (e.g., 6991_xprimes.fasta). '
             'If provided, X primes will be detected using BLAST alignment.'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--prefix',
        default='pretelomeric_labels',
        help='Prefix for output files (default: pretelomeric_labels)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for BLAST (default: 4)'
    )
    parser.add_argument(
        '--min-pident',
        type=float,
        default=80.0,
        help='Minimum percent identity for BLAST hits (default: 80.0)'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=100,
        help='Minimum alignment length for BLAST hits (default: 100)'
    )
    parser.add_argument(
        '--evalue',
        type=float,
        default=1e-10,
        help='E-value threshold for BLAST (default: 1e-10)'
    )
    parser.add_argument(
        '--probe',
        default=None,
        help='Y prime probe FASTA file for verification (e.g., references/probe.fasta). '
             'If provided, verifies detected Y prime counts match probe-based expected counts.'
    )
    parser.add_argument(
        '--skip-probe-verification',
        action='store_true',
        help='Skip Y prime probe verification even if probe file is provided'
    )

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Define output files
    anchor_blast_file = os.path.join(args.output_dir, f"{args.prefix}_anchor_blast.txt")
    yprime_blast_file = os.path.join(args.output_dir, f"{args.prefix}_yprime_blast.txt")
    xprime_blast_file = os.path.join(args.output_dir, f"{args.prefix}_xprime_blast.txt")
    gff_file = os.path.join(args.output_dir, f"{args.prefix}.gff3")
    bed_file = os.path.join(args.output_dir, f"{args.prefix}.bed")
    bed_simp_file = os.path.join(args.output_dir, f"{args.prefix}_simp.bed")
    tsv_file = os.path.join(args.output_dir, f"{args.prefix}.tsv")
    visualization_file = os.path.join(args.output_dir, f"{args.prefix}_structure.txt")

    print("=" * 80)
    print("Pre-telomeric Region Labeling Pipeline")
    print("=" * 80)
    print(f"Reference: {args.reference}")
    print(f"Anchors: {args.anchors}")
    print(f"Y primes: {args.yprimes}")
    print(f"X primes: {args.xprimes if args.xprimes else 'Not provided (X prime detection skipped)'}")
    print(f"Output directory: {args.output_dir}")
    print()

    # Step 1: Run BLAST for anchors
    print("Step 1: Running BLAST for anchor regions...")
    run_blast(args.anchors, args.reference, anchor_blast_file,
              evalue=args.evalue, num_threads=args.threads)
    print(f"Anchor BLAST results: {anchor_blast_file}")
    print()

    # Step 2: Run BLAST for Y primes
    print("Step 2: Running BLAST for Y prime regions...")
    run_blast(args.yprimes, args.reference, yprime_blast_file,
              evalue=args.evalue, num_threads=args.threads)
    print(f"Y prime BLAST results: {yprime_blast_file}")
    print()

    # Step 2b: Run BLAST for X primes (if provided)
    xprime_df = pd.DataFrame()
    if args.xprimes:
        print("Step 2b: Running BLAST for X prime regions...")
        run_blast(args.xprimes, args.reference, xprime_blast_file,
                  evalue=args.evalue, num_threads=args.threads)
        print(f"X prime BLAST results: {xprime_blast_file}")
        print()

    # Step 3: Parse BLAST results
    print("Step 3: Parsing and filtering BLAST results...")
    anchor_df = parse_blast_results(anchor_blast_file,
                                    min_pident=args.min_pident,
                                    min_length=args.min_length)
    yprime_df = parse_blast_results(yprime_blast_file,
                                    min_pident=MIN_YPRIME_IDENTITY,
                                    min_length=args.min_length)

    # Parse X prime BLAST results if provided
    if args.xprimes and os.path.exists(xprime_blast_file):
        xprime_df = parse_blast_results(xprime_blast_file,
                                        min_pident=MIN_XPRIME_IDENTITY,
                                        min_length=50)  # X primes can be smaller
        xprime_df = filter_xprime_hits(xprime_df)
        print(f"Found {len(xprime_df)} X prime hits (after adjusted identity filtering)")

    print(f"Found {len(anchor_df)} anchor hits")
    print(f"Found {len(yprime_df)} Y prime hits (before filtering and merging)")
    print()

    # Step 4: Assign ANCHORS to chromosome ends with synteny enforcement
    print("Step 4: Assigning ANCHOR regions to chromosome ends (with synteny enforcement)...")
    print("         - Enforcing strict chromosome matching")
    print("         - Requiring ≥95% identity for high-quality anchors")
    # Pass empty yprime_df since Y primes are handled separately now
    chr_end_regions, quality_report = assign_regions_to_chr_ends(
        anchor_df, pd.DataFrame(),
        enforce_synteny=True,
        min_high_quality_pident=95.0
    )
    print(f"Detected anchors for {len(chr_end_regions)} chromosome ends")
    print()

    # Step 5: Assign Y PRIMES to chromosome ends with fragment merging
    print("Step 5: Assigning Y PRIME regions with fragment merging...")
    print(f"         - Min identity: {MIN_YPRIME_IDENTITY}%")
    print(f"         - Min coverage: {MIN_YPRIME_COVERAGE}%")
    print("         - Merging adjacent fragments into complete Y primes")
    chr_end_yprimes, yprime_quality = assign_yprimes_to_chr_ends(yprime_df, chr_end_regions)

    # Add Y primes to chr_end_regions
    for chr_end, yprimes in chr_end_yprimes.items():
        if chr_end in chr_end_regions:
            chr_end_regions[chr_end]['y_prime'] = yprimes
        else:
            # Y prime found without anchor (unusual but possible)
            chr_end_regions[chr_end] = {'anchor': [], 'y_prime': yprimes, 'x_prime': []}

    total_yprimes = sum(len(yprimes) for yprimes in chr_end_yprimes.values())
    print(f"Detected {total_yprimes} Y prime regions across {len(chr_end_yprimes)} chromosome ends")
    print()

    # Step 5b: Assign X PRIMES to chromosome ends (if X prime file provided)
    xprime_quality = []
    xprime_position_warnings = []

    if args.xprimes and len(xprime_df) > 0:
        print("Step 5b: Assigning X PRIME regions...")
        print(f"         - Min identity: {MIN_XPRIME_IDENTITY}%")
        print(f"         - Min adjusted identity: 75% (penalizes fragments)")
        print(f"         - High quality threshold: {MIN_HIGH_QUALITY_IDENTITY}%")
        print("         - Selecting best match per chromosome end")

        chr_end_xprimes, xprime_quality = assign_xprimes_to_chr_ends(xprime_df, chr_end_regions)

        # Get chromosome sizes for position validation
        chr_sizes = get_chr_sizes_from_blast(xprime_df)

        # Validate X prime positions relative to Y primes
        xprime_position_warnings = validate_xprime_positions(chr_end_xprimes, chr_end_yprimes, chr_sizes)

        # Add X primes to chr_end_regions
        for chr_end, xprimes in chr_end_xprimes.items():
            if chr_end in chr_end_regions:
                chr_end_regions[chr_end]['x_prime'] = xprimes
            else:
                chr_end_regions[chr_end] = {'anchor': [], 'y_prime': [], 'x_prime': xprimes}

        total_xprimes = sum(len(xprimes) for xprimes in chr_end_xprimes.values())
        print(f"Detected {total_xprimes} X prime regions across {len(chr_end_xprimes)} chromosome ends")

        # Report quality issues
        mismatches = [q for q in xprime_quality if q['type'] == 'xprime_mismatch']
        low_quality = [q for q in xprime_quality if q['type'] == 'xprime_low_quality']

        if mismatches:
            print()
            print("=" * 80)
            print(f"⚠️  X PRIME SOURCE MISMATCHES: {len(mismatches)}")
            print("=" * 80)
            print("The following X primes matched a different chr end than expected:")
            print("-" * 80)
            print(f"{'Target Chr End':<15} {'Best Match From':<20} {'Adj. Identity':<15}")
            print("-" * 80)
            for m in mismatches:
                print(f"{m['chr_end']:<15} {m['actual_source']}_xprime{'':<5} {m['adjusted_pident']:.1f}%")
            print("-" * 80)

        if low_quality:
            print()
            print(f"⚠️  X PRIME LOW QUALITY MATCHES: {len(low_quality)}")
            for lq in low_quality:
                print(f"   - {lq['message']}")

        if xprime_position_warnings:
            print()
            print("=" * 80)
            print(f"⚠️  X PRIME POSITION ERRORS: {len(xprime_position_warnings)}")
            print("=" * 80)
            print("X prime should be closer to chr middle than Y prime:")
            print("-" * 80)
            for w in xprime_position_warnings:
                print(f"   - {w['message']}")
            print("-" * 80)

        print()
    elif args.xprimes:
        print("Step 5b: X PRIME detection SKIPPED (no valid X prime hits found)")
        print()
    else:
        print("Step 5b: X PRIME detection SKIPPED (no X prime file provided)")
        print()

    # Step 6: Y prime probe verification (if probe file provided)
    probe_verification_passed = True
    probe_mismatches = []

    if args.probe and not args.skip_probe_verification:
        print("Step 6: Verifying Y prime counts with probe...")
        print(f"         - Probe file: {args.probe}")

        # Run probe BLAST
        probe_blast_file = os.path.join(args.output_dir, f"{args.prefix}_probe_blast.txt")
        try:
            run_probe_blast(args.probe, args.reference, probe_blast_file,
                           min_identity=90, num_threads=args.threads)

            # Count expected Y primes from probe
            expected_counts = count_yprimes_from_probe(probe_blast_file, args.reference)

            # Get detected counts
            detected_counts = {chr_end: len(yprimes) for chr_end, yprimes in chr_end_yprimes.items()}

            # Verify counts match
            probe_verification_passed, probe_mismatches = verify_yprime_counts(detected_counts, expected_counts)

            # Print verification results
            total_expected = sum(expected_counts.values())
            print(f"         - Expected Y primes (from probe): {total_expected}")
            print(f"         - Detected Y primes: {total_yprimes}")

            if probe_verification_passed:
                print("         ✅ VERIFICATION PASSED: Detected Y prime counts match expected counts!")
            else:
                print("         ❌ VERIFICATION FAILED: Y prime count mismatch detected!")
                print()
                print("=" * 80)
                print("⚠️  Y PRIME COUNT MISMATCH ERROR ⚠️")
                print("=" * 80)
                print()
                print("The following chromosome ends have mismatched Y prime counts:")
                print("-" * 80)
                print(f"{'Chr End':<10} {'Detected':<10} {'Expected':<10} {'Difference':<12}")
                print("-" * 80)
                for m in probe_mismatches:
                    diff_str = f"+{m['difference']}" if m['difference'] > 0 else str(m['difference'])
                    print(f"{m['chr_end']:<10} {m['detected']:<10} {m['expected']:<10} {diff_str:<12}")
                print("-" * 80)
                print()
                print("Possible causes:")
                print("  - Detected > Expected: Over-detection (fragments counted as separate Y primes)")
                print("  - Detected < Expected: Under-detection (Y primes missed or merged incorrectly)")
                print("  - Assembly issues: Reference may have structural variants")
                print()
                print("Recommendations:")
                print("  - Check the Y prime BLAST results manually")
                print("  - Verify assembly quality in affected regions")
                print("  - Adjust detection parameters if needed")
                print("=" * 80)

        except Exception as e:
            print(f"         ⚠️  WARNING: Probe verification failed with error: {e}")
            print("         Continuing without probe verification...")

        print()
    elif args.probe and args.skip_probe_verification:
        print("Step 6: Y prime probe verification SKIPPED (--skip-probe-verification flag)")
        print()
    else:
        print("Step 6: Y prime probe verification SKIPPED (no probe file provided)")
        print("         To enable verification, use --probe references/probe.fasta")
        print()

    ## COMMENTED OUT FOR TESTING - X PRIME INFERENCE
    # # Step 7: Infer X prime regions
    # print("Step 7: Inferring X prime regions...")
    # chr_end_regions = infer_x_prime_regions(chr_end_regions)

    # Print summary
    print("Summary by chromosome end:")
    for chr_end in sorted(chr_end_regions.keys()):
        regions = chr_end_regions[chr_end]
        print(f"  {chr_end}: anchor={len(regions['anchor'])}, "
              f"y_prime={len(regions['y_prime'])}, "
              f"x_prime={len(regions['x_prime'])}")
    print()

    # Print quality report
    print_quality_report(chr_end_regions, quality_report, yprime_quality)

    # Step 7: Write output files
    print("\nStep 7: Writing output files...")
    write_gff3(chr_end_regions, gff_file)
    print(f"GFF3 file: {gff_file}")

    write_bed(chr_end_regions, bed_file)
    print(f"BED file: {bed_file}")

    write_bed_simplified(chr_end_regions, bed_simp_file)
    print(f"Simplified BED file: {bed_simp_file}")

    write_tsv(chr_end_regions, tsv_file)
    print(f"TSV file: {tsv_file}")

    # Step 8: Generate structure visualization (with comparison to 6991 reference)
    print("\nStep 8: Generating chromosome end structure visualization...")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    visualize_script = os.path.join(script_dir, "visualize_chr_end_structure.py")

    # Reference BED file for comparison (6991 strain)
    ref_bed_file = os.path.join(os.path.dirname(script_dir), "references", "6991_final_features.bed")

    if os.path.exists(visualize_script):
        try:
            # Build command with optional comparison
            viz_cmd = [sys.executable, visualize_script, "-i", tsv_file, "-o", visualization_file]

            if os.path.exists(ref_bed_file):
                viz_cmd.extend(["--compare", ref_bed_file])
                print(f"         - Comparing against reference: {ref_bed_file}")

            result = subprocess.run(
                viz_cmd,
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                print(f"Structure visualization: {visualization_file}")
            else:
                print(f"Warning: Visualization failed: {result.stderr}")
        except Exception as e:
            print(f"Warning: Could not run visualization script: {e}")
    else:
        print(f"Warning: Visualization script not found: {visualize_script}")

    # Write quality report to file
    quality_report_file = os.path.join(args.output_dir, f"{args.prefix}_quality_report.txt")
    with open(quality_report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PRE-TELOMERIC REGION LABELING - QUALITY REPORT\n")
        f.write("=" * 80 + "\n\n")

        if quality_report['anchor_mismatches']:
            f.write("ANCHOR CHROMOSOME MISMATCHES (EXCLUDED FROM RESULTS):\n")
            f.write("-" * 80 + "\n")
            for m in quality_report['anchor_mismatches']:
                f.write(f"Query: {m['query']}\n")
                f.write(f"  Expected chromosome: {m['expected_chr']}\n")
                f.write(f"  Found in: {m['found_chr']}\n")
                f.write(f"  Identity: {m['pident']:.2f}% | Coverage: {m['query_coverage']:.2f}% | Bitscore: {m['bitscore']:.2f}\n\n")

        if quality_report['anchor_low_quality']:
            # Separate primary low quality from secondary hits
            primary_lq = [m for m in quality_report['anchor_low_quality'] if 'reason' not in m]
            secondary = [m for m in quality_report['anchor_low_quality'] if 'reason' in m]

            if primary_lq:
                f.write("\nANCHOR LOW QUALITY MATCHES (Expected: ≥95% identity AND ≥95% coverage):\n")
                f.write("-" * 80 + "\n")
                for m in primary_lq:
                    f.write(f"Query: {m['query']}\n")
                    f.write(f"  Chromosome: {m['chr']}\n")
                    f.write(f"  Position: {m['start']:,}-{m['end']:,}\n")
                    f.write(f"  Identity: {m['pident']:.2f}% | Coverage: {m['query_coverage']:.2f}% | Bitscore: {m['bitscore']:.2f}\n")
                    if 'issues' in m:
                        f.write(f"  Issues: {m['issues']}\n")
                    f.write("\n")

            if secondary:
                f.write(f"\nSECONDARY/DUPLICATE HITS FILTERED OUT: {len(secondary)} total\n")
                f.write("(Only best hit per chromosome end kept)\n\n")

        if quality_report['anchor_missing']:
            f.write("\nMISSING ANCHORS:\n")
            f.write("-" * 80 + "\n")
            f.write(f"Chromosome ends without anchors: {', '.join(quality_report['anchor_missing'])}\n\n")

        if quality_report['yprime_mismatches']:
            f.write("\nY PRIME CHROMOSOME MISMATCHES (EXCLUDED FROM RESULTS):\n")
            f.write("-" * 80 + "\n")
            for m in quality_report['yprime_mismatches']:
                f.write(f"Query: {m['query']}\n")
                f.write(f"  Expected: {m['expected_chr']} | Found: {m['found_chr']}\n")
                f.write(f"  Identity: {m['pident']:.2f}%\n\n")

        if quality_report['yprime_low_quality']:
            f.write(f"\nY PRIME LOW QUALITY MATCHES (<80% identity): {len(quality_report['yprime_low_quality'])} total\n")

        if not any([quality_report['anchor_mismatches'],
                   quality_report['anchor_low_quality'],
                   quality_report['anchor_missing'],
                   quality_report['yprime_mismatches']]):
            f.write("\nNo quality issues detected!\n")
            f.write("All anchors found with high identity (≥95%) on expected chromosomes.\n")

        # Add probe verification results
        if args.probe and not args.skip_probe_verification:
            f.write("\n" + "=" * 80 + "\n")
            f.write("Y PRIME PROBE VERIFICATION\n")
            f.write("=" * 80 + "\n\n")

            if probe_verification_passed:
                f.write("✅ VERIFICATION PASSED\n")
                f.write("Detected Y prime counts match expected counts from probe.\n")
            else:
                f.write("❌ VERIFICATION FAILED\n")
                f.write("Y prime count mismatches detected:\n\n")
                f.write(f"{'Chr End':<10} {'Detected':<10} {'Expected':<10} {'Difference':<12}\n")
                f.write("-" * 42 + "\n")
                for m in probe_mismatches:
                    diff_str = f"+{m['difference']}" if m['difference'] > 0 else str(m['difference'])
                    f.write(f"{m['chr_end']:<10} {m['detected']:<10} {m['expected']:<10} {diff_str:<12}\n")

    print(f"Quality report: {quality_report_file}")
    print()

    # Final status
    if not probe_verification_passed:
        print("=" * 80)
        print("⚠️  Pipeline completed with Y PRIME COUNT MISMATCH WARNING")
        print("    Review the mismatches above and verify assembly quality.")
        print("=" * 80)
    else:
        print("=" * 80)
        print("Pipeline completed successfully!")
        print("=" * 80)


if __name__ == '__main__':
    main()
