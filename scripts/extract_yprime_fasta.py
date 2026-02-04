#!/usr/bin/env python3
"""
Extract Y prime sequences from a reference genome based on labeled regions.

Creates a FASTA file with unique Y prime sequences, merging identical sequences
and tracking all locations where each Y prime was found.

Header format matches 6991 reference style for pipeline compatibility:
>Y_Prime_chr4R1;chr12R3#Short/Solo/ID1_Blue

The format after # is: Size/Type/ID_Color
- Size: Short (5-6.4kb) or Long (6.4-6.9kb)
- Type: Solo (single location) or Tandem (multiple locations)
- ID_Color: Y prime group ID and associated color
"""

import argparse
import os
import sys
import pandas as pd
from collections import defaultdict


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def load_reference(fasta_path):
    """Load reference genome into memory."""
    print(f"Loading reference genome: {fasta_path}")
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]  # Get first word after >
                current_seq = []
            else:
                current_seq.append(line)

        if current_name:
            sequences[current_name] = ''.join(current_seq)

    print(f"  Loaded {len(sequences)} chromosomes")
    return sequences


def extract_sequence(sequences, chrom, start, end, strand):
    """Extract sequence from reference, handling strand orientation."""
    if chrom not in sequences:
        print(f"  Warning: Chromosome {chrom} not found in reference")
        return None

    # Convert to 0-based coordinates
    seq = sequences[chrom][start-1:end]

    if strand == '-':
        seq = reverse_complement(seq)

    return seq.upper()


def load_yprime_regions(tsv_path):
    """Load Y prime regions from labeled TSV file."""
    print(f"Loading Y prime regions from: {tsv_path}")
    df = pd.read_csv(tsv_path, sep='\t')

    # Filter for Y primes only
    type_col = 'region_type' if 'region_type' in df.columns else 'type'
    yprimes = df[df[type_col] == 'y_prime'].copy()

    print(f"  Found {len(yprimes)} Y prime regions")
    return yprimes


def parse_yprime_id_from_source(source):
    """
    Parse Y' prime ID and color from the source column.

    Source format: Y_Prime_chr10L1#Long/Solo/ID6_Purple-Dark
    Returns: (id_str, color) e.g., ('ID6', 'Purple-Dark')
    """
    try:
        # Split on '#' to get the metadata part
        if '#' in source:
            metadata = source.split('#')[1]
            # Split on '/' to get Size/Type/ID_Color
            parts = metadata.split('/')
            if len(parts) >= 3:
                # ID_Color is the third part, e.g., "ID6_Purple-Dark"
                id_color = parts[2]
                # Split on first '_' to separate ID from color
                if '_' in id_color:
                    underscore_idx = id_color.index('_')
                    y_prime_id = id_color[:underscore_idx]
                    color = id_color[underscore_idx + 1:]
                    return (y_prime_id, color)
    except (IndexError, ValueError) as e:
        print(f"  Warning: Could not parse ID from source: {source}")

    return (None, None)


def find_unique_yprimes(yprimes_df, sequences):
    """
    Find unique Y prime sequences and track all their locations.

    Returns dict mapping sequence -> {
        'locations': [(chr_end, position_num, identity, length), ...],
        'avg_identity': float,
        'length': int,
        'y_prime_id': str,  # Inherited from BLAST match
        'color': str        # Inherited from BLAST match
    }
    """
    print("Extracting and deduplicating Y prime sequences...")

    # Track Y primes per chr_end for numbering
    chr_end_counts = defaultdict(int)

    # Store sequences with their metadata
    sequence_info = defaultdict(lambda: {'locations': [], 'identities': [], 'length': 0,
                                          'y_prime_id': None, 'color': None})

    for _, row in yprimes_df.iterrows():
        chr_end = row['chr_end']
        chrom = row['chr']
        start = int(row['start'])
        end = int(row['end'])
        strand = row['strand']
        length = int(row['length'])
        pident = float(row['pident'])

        # Get source column to extract Y' prime ID
        source = row.get('source', '')

        # Increment position counter for this chr_end
        chr_end_counts[chr_end] += 1
        position_num = chr_end_counts[chr_end]

        # Extract sequence
        seq = extract_sequence(sequences, chrom, start, end, strand)
        if seq is None:
            continue

        # Store with location info
        sequence_info[seq]['locations'].append((chr_end, position_num, pident, length))
        sequence_info[seq]['identities'].append(pident)
        sequence_info[seq]['length'] = length

        # Parse and store Y' prime ID from source (use first match's ID)
        if sequence_info[seq]['y_prime_id'] is None and source:
            y_prime_id, color = parse_yprime_id_from_source(source)
            if y_prime_id:
                sequence_info[seq]['y_prime_id'] = y_prime_id
                sequence_info[seq]['color'] = color

    # Calculate average identity for each unique sequence
    for seq in sequence_info:
        identities = sequence_info[seq]['identities']
        sequence_info[seq]['avg_identity'] = sum(identities) / len(identities)

    print(f"  Found {len(sequence_info)} unique Y prime sequences")
    print(f"  Total Y prime instances: {sum(len(v['locations']) for v in sequence_info.values())}")

    return dict(sequence_info)


def classify_yprime_size(length):
    """Classify Y prime as Short or Long based on size."""
    if 5000 <= length < 6400:
        return "Short"
    elif 6400 <= length <= 6900:
        return "Long"
    else:
        return "Unknown"


def format_header(locations, avg_identity, length, y_prime_id=None, color=None):
    """
    Format FASTA header with all locations and metadata.

    Format matches 6991 reference style:
    >Y_Prime_chr4R1;chr12R3,5;chr14L2#Short/Solo/ID1_Blue

    This format is required for compatibility with get_stats_of_recombination.py
    which parses the y_prime_group field by splitting on '/'.
    """
    # Group locations by chr_end
    chr_end_positions = defaultdict(list)
    for chr_end, pos_num, pident, seq_len in locations:
        chr_end_positions[chr_end].append(pos_num)

    # Format location string
    location_parts = []
    for chr_end in sorted(chr_end_positions.keys()):
        positions = sorted(chr_end_positions[chr_end])
        if len(positions) == 1:
            location_parts.append(f"{chr_end}{positions[0]}")
        else:
            # Multiple positions on same chr_end: chr4R1,2,3
            pos_str = ','.join(str(p) for p in positions)
            location_parts.append(f"{chr_end}{pos_str}")

    location_str = ';'.join(location_parts)

    # Classify size
    size_class = classify_yprime_size(length)

    # Determine if Solo or Tandem based on total number of positions
    total_positions = sum(len(positions) for positions in chr_end_positions.values())
    tandem_type = "Solo" if total_positions == 1 else "Tandem"

    # Use inherited ID and color from BLAST match, or fall back to Unknown
    if y_prime_id is None:
        y_prime_id = "Unknown"
        color = "Gray"
    elif color is None:
        color = "Gray"

    # Format header to match 6991 style: Y_Prime_location#Size/Type/ID_Color
    header = f">Y_Prime_{location_str}#{size_class}/{tandem_type}/{y_prime_id}_{color}"

    return header


def write_fasta(unique_yprimes, output_path):
    """Write unique Y primes to FASTA file."""
    print(f"Writing FASTA to: {output_path}")

    # Sort by number of locations (most common first), then by first chr_end
    sorted_yprimes = sorted(
        unique_yprimes.items(),
        key=lambda x: (-len(x[1]['locations']), x[1]['locations'][0][0])
    )

    with open(output_path, 'w') as f:
        for idx, (seq, info) in enumerate(sorted_yprimes):
            header = format_header(
                info['locations'],
                info['avg_identity'],
                info['length'],
                info.get('y_prime_id'),  # Pass inherited ID from BLAST match
                info.get('color')        # Pass inherited color from BLAST match
            )
            f.write(f"{header}\n")

            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")

    print(f"  Wrote {len(sorted_yprimes)} unique Y prime sequences")


def print_summary(unique_yprimes):
    """Print summary of unique Y primes found."""
    print("\n" + "=" * 80)
    print("Y PRIME SUMMARY")
    print("=" * 80)

    # Count by size class
    short_count = 0
    long_count = 0
    unknown_count = 0

    # Count by number of locations
    single_location = 0
    multi_location = 0

    for seq, info in unique_yprimes.items():
        length = info['length']
        size_class = classify_yprime_size(length)

        if size_class == "Short":
            short_count += 1
        elif size_class == "Long":
            long_count += 1
        else:
            unknown_count += 1

        if len(info['locations']) == 1:
            single_location += 1
        else:
            multi_location += 1

    print(f"\nUnique Y prime sequences: {len(unique_yprimes)}")
    print(f"  - Short (5-6.4kb): {short_count}")
    print(f"  - Long (6.4-6.9kb): {long_count}")
    print(f"  - Unknown size: {unknown_count}")
    print(f"\nLocation distribution:")
    print(f"  - Found at single location: {single_location}")
    print(f"  - Found at multiple locations: {multi_location}")

    # Show multi-location Y primes
    if multi_location > 0:
        print(f"\nY primes found at multiple locations:")
        for seq, info in sorted(unique_yprimes.items(),
                                key=lambda x: -len(x[1]['locations'])):
            if len(info['locations']) > 1:
                locs = [f"{loc[0]}{loc[1]}" for loc in info['locations']]
                print(f"  {info['length']}bp ({info['avg_identity']:.1f}%): {', '.join(locs)}")


def main():
    parser = argparse.ArgumentParser(
        description='Extract unique Y prime sequences to FASTA file'
    )

    parser.add_argument(
        '--labeled-tsv',
        required=True,
        help='Labeled regions TSV from label_pretelomeric_regions.py'
    )
    parser.add_argument(
        '--reference',
        required=True,
        help='Reference FASTA file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output FASTA file for Y prime sequences'
    )
    parser.add_argument(
        '--strain',
        default='unknown',
        help='Strain identifier (for output naming)'
    )

    args = parser.parse_args()

    print("=" * 80)
    print("Y Prime Sequence Extractor")
    print("=" * 80)
    print(f"Labeled TSV: {args.labeled_tsv}")
    print(f"Reference: {args.reference}")
    print(f"Output: {args.output}")
    print()

    # Load reference genome
    sequences = load_reference(args.reference)

    # Load Y prime regions
    yprimes_df = load_yprime_regions(args.labeled_tsv)

    if len(yprimes_df) == 0:
        print("No Y primes found in labeled TSV!")
        sys.exit(1)

    # Find unique Y primes
    unique_yprimes = find_unique_yprimes(yprimes_df, sequences)

    # Write FASTA
    write_fasta(unique_yprimes, args.output)

    # Print summary
    print_summary(unique_yprimes)

    print("\n" + "=" * 80)
    print("Done!")
    print("=" * 80)


if __name__ == '__main__':
    main()
