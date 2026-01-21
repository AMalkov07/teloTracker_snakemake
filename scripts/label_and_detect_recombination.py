#!/usr/bin/env python3
"""
Label pre-telomeric regions AND detect recombination events in a reference genome.

This script extends label_pretelomeric_regions.py by:
1. Tracking the source chromosome end (from 6991) for each detected feature
2. Identifying recombination breakpoints where feature sources change
3. Generating a recombination report showing where the new strain differs from 6991

Recombination detection logic:
- For each chromosome end in the new reference:
  - Detect anchor → which 6991 anchor matches best? (should match expected chr end)
  - Detect X prime → which 6991 X prime matches best?
  - Detect Y prime(s) → which 6991 Y prime(s) match best?
- Compare sources:
  - If all features match same 6991 chr end → No recombination
  - If sources differ → Recombination breakpoint between features with different sources
"""

import argparse
import os
import subprocess
import sys
import pandas as pd
from collections import defaultdict

# Import from existing modules
from yprime_detection_v2 import (
    assign_yprimes_to_chr_ends,
    MIN_YPRIME_IDENTITY,
    MIN_YPRIME_COVERAGE
)

from xprime_detection import (
    assign_xprimes_to_chr_ends,
    filter_xprime_hits,
    validate_xprime_positions,
    get_chr_sizes_from_blast,
    MIN_XPRIME_IDENTITY,
    MIN_HIGH_QUALITY_IDENTITY
)

# Import functions from original labeling script
from label_pretelomeric_regions import (
    run_blast,
    parse_blast_results,
    extract_chr_end_from_query,
    check_chromosome_match,
    assign_regions_to_chr_ends,
    write_gff3,
    write_bed,
    write_bed_simplified,
    write_tsv,
    print_quality_report,
    run_probe_blast,
    count_yprimes_from_probe,
    verify_yprime_counts
)


def extract_source_chr_end(source_id):
    """
    Extract the source chromosome end from a feature's source ID.

    Examples:
        chr1L_anchor -> chr1L
        Y_Prime_chr2L1#Short/Solo/ID4_Green -> chr2L
        chr3R_xprime -> chr3R

    Args:
        source_id: Source identifier string

    Returns:
        str: Source chromosome end (e.g., 'chr1L') or None
    """
    if source_id is None:
        return None

    source_str = str(source_id)

    # Anchor format: chr1L_anchor
    if '_anchor' in source_str:
        return source_str.replace('_anchor', '')

    # X prime format: chr1L_xprime or chr1L_x_prime
    if '_xprime' in source_str.lower() or '_x_prime' in source_str.lower():
        parts = source_str.lower().replace('_xprime', '').replace('_x_prime', '')
        # Normalize capitalization
        if len(parts) >= 4:
            return parts[:3] + parts[3].upper() + parts[4:] if len(parts) > 4 else parts[:3] + parts[3].upper()

    # Y prime format: Y_Prime_chr2L1#... or Y_Prime_chr2L1,2,3#...
    if source_str.startswith('Y_Prime_'):
        parts = source_str.split('#')[0].replace('Y_Prime_', '')
        # Extract first chr_end from patterns like chr2L1 or chr4R1,2,3
        chr_parts = parts.split(';')[0].split(',')[0]

        # Find where the arm letter is (look for L or R)
        for i, char in enumerate(chr_parts):
            if char in ['L', 'R']:
                return chr_parts[:i+1]

    return None


def detect_recombination_events(chr_end_regions):
    """
    Detect recombination events by comparing feature sources within each chromosome end.

    For each chromosome end, compares the source (6991 chr end) of:
    - Anchor
    - X prime
    - Y prime(s)

    If sources differ, identifies recombination breakpoints.

    Args:
        chr_end_regions: Dict of chromosome end regions with features

    Returns:
        list: List of recombination event dicts
    """
    recombination_events = []

    for chr_end in sorted(chr_end_regions.keys()):
        regions = chr_end_regions[chr_end]

        # Get sources for each feature type
        sources = {
            'anchor': None,
            'x_prime': None,
            'y_primes': []
        }

        # Get anchor source
        if regions['anchor']:
            anchor = regions['anchor'][0]  # Best anchor
            sources['anchor'] = extract_source_chr_end(anchor.get('source'))

        # Get X prime source
        if regions['x_prime']:
            xprime = regions['x_prime'][0]  # Best X prime
            sources['x_prime'] = extract_source_chr_end(xprime.get('source'))

        # Get Y prime sources (can have multiple)
        for yprime in regions['y_prime']:
            yp_source = extract_source_chr_end(yprime.get('source'))
            if yp_source:
                sources['y_primes'].append({
                    'source': yp_source,
                    'position': yprime.get('start', 0),
                    'yprime_id': yprime.get('source')
                })

        # Expected source is the chromosome end itself (e.g., chr1L should have chr1L features)
        expected_source = chr_end

        # Check for recombination
        event = {
            'chr_end': chr_end,
            'expected_source': expected_source,
            'anchor_source': sources['anchor'],
            'xprime_source': sources['x_prime'],
            'yprime_sources': [yp['source'] for yp in sources['y_primes']],
            'has_recombination': False,
            'recombination_type': None,
            'breakpoint_location': None,
            'details': []
        }

        # Check anchor match
        if sources['anchor'] and sources['anchor'] != expected_source:
            event['has_recombination'] = True
            event['details'].append(f"Anchor from {sources['anchor']} (expected {expected_source})")

        # Check X prime match
        if sources['x_prime'] and sources['x_prime'] != expected_source:
            event['has_recombination'] = True
            event['details'].append(f"X prime from {sources['x_prime']} (expected {expected_source})")

            # If anchor matches but X prime doesn't, breakpoint is between anchor and X prime
            if sources['anchor'] == expected_source:
                event['breakpoint_location'] = 'between_anchor_and_xprime'
                event['recombination_type'] = 'xprime_swap'

        # Check Y prime matches
        mismatched_yprimes = [yp for yp in sources['y_primes'] if yp['source'] != expected_source]
        if mismatched_yprimes:
            event['has_recombination'] = True
            for yp in mismatched_yprimes:
                event['details'].append(f"Y prime from {yp['source']} (expected {expected_source})")

            # Determine breakpoint location
            if sources['anchor'] == expected_source and sources['x_prime'] == expected_source:
                event['breakpoint_location'] = 'between_xprime_and_yprime'
                event['recombination_type'] = 'yprime_swap'
            elif sources['anchor'] == expected_source:
                event['breakpoint_location'] = 'between_anchor_and_yprime'
                event['recombination_type'] = 'xprime_and_yprime_swap'

        # Check for complete chromosome end swap
        all_sources = [s for s in [sources['anchor'], sources['x_prime']] + [yp['source'] for yp in sources['y_primes']] if s]
        unique_sources = set(all_sources)

        if len(unique_sources) == 1 and expected_source not in unique_sources:
            # All features from same source, but not expected
            event['recombination_type'] = 'complete_chr_end_swap'
            event['breakpoint_location'] = 'before_anchor'
            event['swap_source'] = list(unique_sources)[0]
        elif len(unique_sources) > 1:
            # Multiple different sources
            event['recombination_type'] = 'complex_recombination'

        recombination_events.append(event)

    return recombination_events


def write_recombination_report(recombination_events, output_file):
    """
    Write recombination events to a TSV file.

    Args:
        recombination_events: List of recombination event dicts
        output_file: Path to output TSV file
    """
    rows = []

    for event in recombination_events:
        rows.append({
            'chr_end': event['chr_end'],
            'expected_source': event['expected_source'],
            'anchor_source': event['anchor_source'],
            'xprime_source': event['xprime_source'],
            'yprime_sources': ';'.join(event['yprime_sources']) if event['yprime_sources'] else '',
            'has_recombination': event['has_recombination'],
            'recombination_type': event['recombination_type'] or '',
            'breakpoint_location': event['breakpoint_location'] or '',
            'details': ' | '.join(event['details']) if event['details'] else ''
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)


def print_recombination_summary(recombination_events):
    """
    Print a summary of detected recombination events.

    Args:
        recombination_events: List of recombination event dicts
    """
    print("\n" + "=" * 80)
    print("RECOMBINATION DETECTION SUMMARY")
    print("=" * 80)

    # Count events by type
    events_with_recomb = [e for e in recombination_events if e['has_recombination']]
    events_no_recomb = [e for e in recombination_events if not e['has_recombination']]

    print(f"\nChromosome ends analyzed: {len(recombination_events)}")
    print(f"  - No recombination detected: {len(events_no_recomb)}")
    print(f"  - Recombination detected: {len(events_with_recomb)}")

    if events_with_recomb:
        # Group by type
        by_type = defaultdict(list)
        for e in events_with_recomb:
            by_type[e['recombination_type'] or 'unknown'].append(e)

        print("\nRecombination events by type:")
        for rtype, events in sorted(by_type.items()):
            print(f"  - {rtype}: {len(events)}")
            for e in events:
                print(f"      {e['chr_end']}: {' | '.join(e['details'])}")

    print("\n" + "=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description='Label pre-telomeric regions AND detect recombination events'
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
        help='X prime sequences FASTA file (e.g., 6991_xprimes.fasta)'
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
        help='Y prime probe FASTA file for verification'
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
    recomb_file = os.path.join(args.output_dir, f"{args.prefix}_recombination.tsv")
    visualization_file = os.path.join(args.output_dir, f"{args.prefix}_structure.txt")

    print("=" * 80)
    print("Pre-telomeric Region Labeling with Recombination Detection")
    print("=" * 80)
    print(f"Reference: {args.reference}")
    print(f"Anchors: {args.anchors}")
    print(f"Y primes: {args.yprimes}")
    print(f"X primes: {args.xprimes if args.xprimes else 'Not provided'}")
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

    if args.xprimes and os.path.exists(xprime_blast_file):
        xprime_df = parse_blast_results(xprime_blast_file,
                                        min_pident=MIN_XPRIME_IDENTITY,
                                        min_length=50)
        xprime_df = filter_xprime_hits(xprime_df)
        print(f"Found {len(xprime_df)} X prime hits")

    print(f"Found {len(anchor_df)} anchor hits")
    print(f"Found {len(yprime_df)} Y prime hits")
    print()

    # Step 4: Assign anchors to chromosome ends
    print("Step 4: Assigning ANCHOR regions to chromosome ends...")
    chr_end_regions, quality_report = assign_regions_to_chr_ends(
        anchor_df, pd.DataFrame(),
        enforce_synteny=True,
        min_high_quality_pident=95.0
    )
    print(f"Detected anchors for {len(chr_end_regions)} chromosome ends")
    print()

    # Step 5: Assign Y primes to chromosome ends
    print("Step 5: Assigning Y PRIME regions with fragment merging...")
    chr_end_yprimes, yprime_quality = assign_yprimes_to_chr_ends(yprime_df, chr_end_regions)

    for chr_end, yprimes in chr_end_yprimes.items():
        if chr_end in chr_end_regions:
            chr_end_regions[chr_end]['y_prime'] = yprimes
        else:
            chr_end_regions[chr_end] = {'anchor': [], 'y_prime': yprimes, 'x_prime': []}

    total_yprimes = sum(len(yprimes) for yprimes in chr_end_yprimes.values())
    print(f"Detected {total_yprimes} Y prime regions across {len(chr_end_yprimes)} chromosome ends")
    print()

    # Step 5b: Assign X primes (if provided)
    if args.xprimes and len(xprime_df) > 0:
        print("Step 5b: Assigning X PRIME regions...")
        chr_end_xprimes, xprime_quality = assign_xprimes_to_chr_ends(xprime_df, chr_end_regions)

        for chr_end, xprimes in chr_end_xprimes.items():
            if chr_end in chr_end_regions:
                chr_end_regions[chr_end]['x_prime'] = xprimes
            else:
                chr_end_regions[chr_end] = {'anchor': [], 'y_prime': [], 'x_prime': xprimes}

        total_xprimes = sum(len(xprimes) for xprimes in chr_end_xprimes.values())
        print(f"Detected {total_xprimes} X prime regions across {len(chr_end_xprimes)} chromosome ends")
        print()

    # Step 6: Y prime probe verification (if provided)
    if args.probe and not args.skip_probe_verification:
        print("Step 6: Verifying Y prime counts with probe...")
        probe_blast_file = os.path.join(args.output_dir, f"{args.prefix}_probe_blast.txt")
        try:
            run_probe_blast(args.probe, args.reference, probe_blast_file,
                           min_identity=90, num_threads=args.threads)
            expected_counts = count_yprimes_from_probe(probe_blast_file, args.reference)
            detected_counts = {chr_end: len(yprimes) for chr_end, yprimes in chr_end_yprimes.items()}
            probe_verification_passed, probe_mismatches = verify_yprime_counts(detected_counts, expected_counts)

            if probe_verification_passed:
                print("         ✅ VERIFICATION PASSED")
            else:
                print("         ❌ VERIFICATION FAILED - see quality report")
        except Exception as e:
            print(f"         ⚠️  WARNING: Probe verification failed: {e}")
        print()
    else:
        print("Step 6: Y prime probe verification SKIPPED")
        print()

    # Step 7: DETECT RECOMBINATION EVENTS
    print("Step 7: Detecting recombination events...")
    recombination_events = detect_recombination_events(chr_end_regions)

    # Print recombination summary
    print_recombination_summary(recombination_events)

    # Write recombination report
    write_recombination_report(recombination_events, recomb_file)
    print(f"Recombination report: {recomb_file}")
    print()

    # Step 8: Write output files
    print("Step 8: Writing output files...")

    write_gff3(chr_end_regions, gff_file)
    print(f"GFF3 file: {gff_file}")

    write_bed(chr_end_regions, bed_file)
    print(f"BED file: {bed_file}")

    write_bed_simplified(chr_end_regions, bed_simp_file)
    print(f"Simplified BED file: {bed_simp_file}")

    write_tsv(chr_end_regions, tsv_file)
    print(f"TSV file: {tsv_file}")

    # Step 9: Generate structure visualization
    print("\nStep 9: Generating chromosome end structure visualization...")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    visualize_script = os.path.join(script_dir, "visualize_chr_end_structure.py")
    ref_bed_file = os.path.join(os.path.dirname(script_dir), "references", "6991_final_features.bed")

    if os.path.exists(visualize_script):
        try:
            viz_cmd = [sys.executable, visualize_script, "-i", tsv_file, "-o", visualization_file]
            if os.path.exists(ref_bed_file):
                viz_cmd.extend(["--compare", ref_bed_file])
            result = subprocess.run(viz_cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print(f"Structure visualization: {visualization_file}")
            else:
                print(f"Warning: Visualization failed: {result.stderr}")
        except Exception as e:
            print(f"Warning: Could not run visualization script: {e}")

    # Print quality report
    print_quality_report(chr_end_regions, quality_report, yprime_quality)

    # Final summary
    events_with_recomb = [e for e in recombination_events if e['has_recombination']]

    print("\n" + "=" * 80)
    if events_with_recomb:
        print(f"Pipeline completed with {len(events_with_recomb)} RECOMBINATION EVENT(S) detected!")
        print("Review the recombination report for details.")
    else:
        print("Pipeline completed successfully!")
        print("No recombination events detected - all features match expected sources.")
    print("=" * 80)


if __name__ == '__main__':
    main()
