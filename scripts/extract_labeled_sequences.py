#!/usr/bin/env python3
"""
Extract sequences for labeled pre-telomeric regions.

This script reads labeled regions from TSV and extracts corresponding
sequences from the reference FASTA file.
"""

import argparse
import pandas as pd
import sys


def load_labels(tsv_file):
    """Load labeled regions from TSV file."""
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading TSV file: {e}")
        sys.exit(1)


def parse_fasta(fasta_file):
    """Parse FASTA file into a dictionary."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def extract_sequence(sequences, chrom, start, end, strand):
    """Extract a subsequence from the sequences dict."""
    if chrom not in sequences:
        return None

    seq = sequences[chrom][start-1:end]  # Convert to 0-based

    # Reverse complement if on minus strand
    if strand == '-':
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        seq = ''.join(complement.get(base, base) for base in reversed(seq))

    return seq


def write_fasta(output_file, sequences_dict):
    """Write sequences to FASTA file."""
    with open(output_file, 'w') as f:
        for seq_id, seq in sequences_dict.items():
            f.write(f">{seq_id}\n")
            # Write in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Extract sequences for labeled pre-telomeric regions'
    )

    parser.add_argument(
        '--input-tsv',
        required=True,
        help='TSV file from label_pretelomeric_regions.py'
    )
    parser.add_argument(
        '--reference',
        required=True,
        help='Reference FASTA file'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for extracted sequences'
    )
    parser.add_argument(
        '--prefix',
        default='extracted',
        help='Prefix for output files (default: extracted)'
    )
    parser.add_argument(
        '--region-types',
        nargs='+',
        choices=['anchor', 'x_prime', 'y_prime'],
        default=['anchor', 'x_prime', 'y_prime'],
        help='Which region types to extract (default: all)'
    )
    parser.add_argument(
        '--chr-ends',
        nargs='+',
        help='Specific chromosome ends to extract (e.g., chr1L chr1R). If not specified, extracts all.'
    )

    args = parser.parse_args()

    import os
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 80)
    print("Extracting labeled region sequences")
    print("=" * 80)
    print(f"Input TSV: {args.input_tsv}")
    print(f"Reference: {args.reference}")
    print(f"Output directory: {args.output_dir}")
    print(f"Region types: {', '.join(args.region_types)}")
    print()

    # Load labels
    print("Loading labels...")
    df = load_labels(args.input_tsv)

    # Filter by region type
    df = df[df['region_type'].isin(args.region_types)]

    # Filter by chromosome ends if specified
    if args.chr_ends:
        df = df[df['chr_end'].isin(args.chr_ends)]

    print(f"Found {len(df)} regions to extract")
    print()

    # Load reference sequences
    print("Loading reference sequences...")
    sequences = parse_fasta(args.reference)
    print(f"Loaded {len(sequences)} sequences from reference")
    print()

    # Extract sequences for each region type
    for region_type in args.region_types:
        subset = df[df['region_type'] == region_type]

        if len(subset) == 0:
            print(f"No {region_type} regions found, skipping...")
            continue

        extracted = {}
        skipped = 0

        for _, row in subset.iterrows():
            seq_id = f"{row['chr_end']}_{region_type}_{row['chr']}:{row['start']}-{row['end']}({row['strand']})"

            seq = extract_sequence(sequences, row['chr'], row['start'], row['end'], row['strand'])

            if seq:
                extracted[seq_id] = seq
            else:
                skipped += 1
                print(f"Warning: Could not extract sequence for {seq_id}")

        if extracted:
            output_file = os.path.join(args.output_dir, f"{args.prefix}_{region_type}.fasta")
            write_fasta(output_file, extracted)
            print(f"{region_type:10s}: Extracted {len(extracted)} sequences -> {output_file}")
            if skipped > 0:
                print(f"             (Skipped {skipped} due to missing chromosomes)")

    print()
    print("=" * 80)
    print("Extraction complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
