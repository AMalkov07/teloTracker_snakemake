#!/usr/bin/env python3
"""
Visualize and summarize pre-telomeric region labels.

This script reads the TSV output from label_pretelomeric_regions.py
and generates summary statistics and simple visualizations.
"""

import argparse
import pandas as pd
import sys
from collections import defaultdict


def load_labels(tsv_file):
    """Load labeled regions from TSV file."""
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading TSV file: {e}")
        sys.exit(1)


def print_summary_stats(df):
    """Print summary statistics about labeled regions."""
    print("=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    print()

    # Overall counts
    print("Region counts by type:")
    region_counts = df['region_type'].value_counts()
    for region_type, count in region_counts.items():
        print(f"  {region_type:15s}: {count:3d}")
    print()

    # Chromosome end coverage
    chr_ends = df['chr_end'].nunique()
    print(f"Total chromosome ends with labels: {chr_ends}")
    print()

    # Region size statistics
    print("Region size statistics (bp):")
    for region_type in ['anchor', 'x_prime', 'y_prime']:
        subset = df[df['region_type'] == region_type]
        if len(subset) > 0:
            print(f"\n  {region_type}:")
            print(f"    Count:   {len(subset)}")
            print(f"    Mean:    {subset['length'].mean():.0f}")
            print(f"    Median:  {subset['length'].median():.0f}")
            print(f"    Min:     {subset['length'].min()}")
            print(f"    Max:     {subset['length'].max()}")
    print()

    # Quality metrics for BLAST-derived regions
    blast_regions = df[df['pident'].notna()]
    if len(blast_regions) > 0:
        print("BLAST alignment quality (for anchor and Y prime regions):")
        print(f"  Percent identity:")
        print(f"    Mean:   {blast_regions['pident'].mean():.2f}%")
        print(f"    Median: {blast_regions['pident'].median():.2f}%")
        print(f"    Min:    {blast_regions['pident'].min():.2f}%")
        print(f"    Max:    {blast_regions['pident'].max():.2f}%")
        print()
        print(f"  Query coverage:")
        print(f"    Mean:   {blast_regions['query_coverage'].mean():.2f}%")
        print(f"    Median: {blast_regions['query_coverage'].median():.2f}%")
        print(f"    Min:    {blast_regions['query_coverage'].min():.2f}%")
        print(f"    Max:    {blast_regions['query_coverage'].max():.2f}%")
        print()


def print_chr_end_details(df):
    """Print detailed information per chromosome end."""
    print("=" * 80)
    print("DETAILED BREAKDOWN BY CHROMOSOME END")
    print("=" * 80)
    print()

    chr_ends = sorted(df['chr_end'].unique())

    for chr_end in chr_ends:
        subset = df[df['chr_end'] == chr_end]

        print(f"{chr_end}:")

        for region_type in ['anchor', 'x_prime', 'y_prime']:
            regions = subset[subset['region_type'] == region_type]

            if len(regions) == 0:
                print(f"  {region_type:10s}: NOT FOUND")
            else:
                for _, region in regions.iterrows():
                    if pd.notna(region['pident']):
                        print(f"  {region_type:10s}: {region['start']:>10,} - {region['end']:>10,} "
                              f"({region['length']:>6,} bp) | "
                              f"identity: {region['pident']:5.1f}%, "
                              f"coverage: {region['query_coverage']:5.1f}%")
                    else:
                        print(f"  {region_type:10s}: {region['start']:>10,} - {region['end']:>10,} "
                              f"({region['length']:>6,} bp) | "
                              f"inferred ({region['source']})")

        print()


def check_quality(df):
    """Check for potential quality issues."""
    print("=" * 80)
    print("QUALITY CHECKS")
    print("=" * 80)
    print()

    issues = []

    # Check for chromosome ends missing anchors
    chr_ends_all = set(df['chr_end'].unique())
    chr_ends_with_anchor = set(df[df['region_type'] == 'anchor']['chr_end'].unique())
    missing_anchors = chr_ends_all - chr_ends_with_anchor

    if missing_anchors:
        issues.append(f"Chromosome ends without anchor: {', '.join(sorted(missing_anchors))}")

    # Check for low percent identity anchors
    anchors = df[df['region_type'] == 'anchor']
    low_pident = anchors[anchors['pident'] < 80]
    if len(low_pident) > 0:
        issues.append(f"{len(low_pident)} anchor(s) with <80% identity (may need manual review)")

    # Check for low coverage anchors
    low_coverage = anchors[anchors['query_coverage'] < 70]
    if len(low_coverage) > 0:
        issues.append(f"{len(low_coverage)} anchor(s) with <70% query coverage (partial alignments)")

    # Check for very short X primes
    xprimes = df[df['region_type'] == 'x_prime']
    short_xprimes = xprimes[xprimes['length'] < 500]
    if len(short_xprimes) > 0:
        issues.append(f"{len(short_xprimes)} X prime region(s) shorter than 500bp")

    # Check for overlapping regions on same chr_end
    for chr_end in df['chr_end'].unique():
        subset = df[df['chr_end'] == chr_end].sort_values('start')
        for i in range(len(subset) - 1):
            if subset.iloc[i]['end'] > subset.iloc[i+1]['start']:
                issues.append(f"Overlapping regions on {chr_end}: "
                            f"{subset.iloc[i]['region_type']} and {subset.iloc[i+1]['region_type']}")

    if issues:
        print("Potential issues found:")
        for i, issue in enumerate(issues, 1):
            print(f"  {i}. {issue}")
        print()
        print("Note: Some issues may be expected (e.g., missing Y primes on certain chr ends)")
    else:
        print("No major quality issues detected!")

    print()


def export_summary(df, output_file):
    """Export a summary table to TSV."""
    summary_rows = []

    for chr_end in sorted(df['chr_end'].unique()):
        subset = df[df['chr_end'] == chr_end]

        row = {'chr_end': chr_end}

        for region_type in ['anchor', 'x_prime', 'y_prime']:
            regions = subset[subset['region_type'] == region_type]

            if len(regions) > 0:
                region = regions.iloc[0]  # Take first if multiple
                row[f'{region_type}_start'] = region['start']
                row[f'{region_type}_end'] = region['end']
                row[f'{region_type}_length'] = region['length']

                if pd.notna(region['pident']):
                    row[f'{region_type}_pident'] = region['pident']
                    row[f'{region_type}_coverage'] = region['query_coverage']
            else:
                row[f'{region_type}_start'] = None
                row[f'{region_type}_end'] = None
                row[f'{region_type}_length'] = None
                row[f'{region_type}_pident'] = None
                row[f'{region_type}_coverage'] = None

        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(output_file, sep='\t', index=False)
    print(f"Summary table exported to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize and summarize pre-telomeric region labels'
    )

    parser.add_argument(
        '--input-tsv',
        required=True,
        help='TSV file from label_pretelomeric_regions.py'
    )
    parser.add_argument(
        '--export-summary',
        help='Export summary table to this TSV file (optional)'
    )

    args = parser.parse_args()

    # Load data
    df = load_labels(args.input_tsv)

    # Print summaries
    print_summary_stats(df)
    print_chr_end_details(df)
    check_quality(df)

    # Export summary if requested
    if args.export_summary:
        export_summary(df, args.export_summary)
        print()

    print("=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
