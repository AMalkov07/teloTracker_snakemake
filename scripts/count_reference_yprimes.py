#!/usr/bin/env python3
"""
Count expected Y primes per chromosome end in a reference genome.

This script uses the Y prime probe sequence to count how many Y primes exist
at each chromosome end in the reference assembly. This provides the expected
Y prime counts that can be compared against read-level Y prime detection.

Methods:
1. BLAST the probe against the reference to find all Y prime locations
2. Group Y prime hits by chromosome and position (telomere-proximal vs centromere-proximal)
3. Assign Y primes to L or R arm based on position relative to chromosome center
4. Output expected Y prime counts per chr_end
"""

import argparse
import os
import subprocess
import sys
import pandas as pd
from collections import defaultdict


def run_blast_probe(probe_fasta, reference_fasta, output_file, min_identity=90, num_threads=4):
    """
    BLAST the Y prime probe against the reference genome.

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

    print(f"Running BLAST: {' '.join(blast_cmd)}")
    subprocess.run(blast_cmd, check=True)

    return output_file


def parse_blast_results(blast_file, min_length=80):
    """
    Parse BLAST results and filter for high-quality probe matches.

    Args:
        blast_file: Path to BLAST output
        min_length: Minimum alignment length

    Returns:
        DataFrame with filtered BLAST results
    """
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen'
    ]

    if os.path.getsize(blast_file) == 0:
        print("Warning: No BLAST hits found")
        return pd.DataFrame(columns=columns)

    df = pd.read_csv(blast_file, sep='\t', names=columns)

    # Filter for high-quality matches
    df = df[df['length'] >= min_length]

    # Calculate position on chromosome (use midpoint of hit)
    df['hit_position'] = (df['sstart'] + df['send']) / 2

    # Determine strand
    df['strand'] = df.apply(lambda row: '+' if row['sstart'] < row['send'] else '-', axis=1)

    return df


def get_chromosome_sizes(reference_fasta):
    """
    Get chromosome sizes from FASTA index or by parsing FASTA.

    Args:
        reference_fasta: Path to reference FASTA

    Returns:
        Dict mapping chromosome name to size
    """
    fai_file = f"{reference_fasta}.fai"

    # Try to use .fai index if exists
    if os.path.exists(fai_file):
        sizes = {}
        with open(fai_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                sizes[parts[0]] = int(parts[1])
        return sizes

    # Otherwise parse FASTA headers
    print("No .fai index found, parsing FASTA...")
    sizes = {}
    try:
        from Bio import SeqIO
        for record in SeqIO.parse(reference_fasta, 'fasta'):
            sizes[record.id] = len(record.seq)
    except ImportError:
        # Fallback: use seqkit
        result = subprocess.run(
            ['seqkit', 'fx2tab', '-n', '-l', reference_fasta],
            capture_output=True, text=True
        )
        for line in result.stdout.strip().split('\n'):
            parts = line.split('\t')
            sizes[parts[0]] = int(parts[1])

    return sizes


def assign_yprimes_to_chr_ends(probe_hits, chr_sizes, anchor_positions=None):
    """
    Assign Y prime probe hits to chromosome ends (L or R arm).

    Strategy:
    - If anchor positions are provided, use them to determine L vs R
    - Otherwise, use chromosome midpoint as boundary
    - Hits in first half → L arm (telomere-proximal for L)
    - Hits in second half → R arm (telomere-proximal for R)

    Args:
        probe_hits: DataFrame with probe BLAST hits
        chr_sizes: Dict of chromosome sizes
        anchor_positions: Optional dict of anchor positions per chr_end

    Returns:
        Dict mapping chr_end to Y prime count
    """
    yprime_counts = defaultdict(int)
    yprime_positions = defaultdict(list)

    for _, hit in probe_hits.iterrows():
        chr_name = hit['sseqid']
        position = hit['hit_position']

        # Extract chromosome number from name (e.g., chr1_extended -> 1)
        chr_num = None
        for part in chr_name.lower().replace('chr', '').replace('_extended', '').split('_'):
            if part.isdigit():
                chr_num = part
                break

        if chr_num is None:
            print(f"Warning: Could not extract chromosome number from {chr_name}")
            continue

        # Get chromosome size
        chr_size = chr_sizes.get(chr_name, 0)
        if chr_size == 0:
            print(f"Warning: Unknown chromosome size for {chr_name}")
            continue

        # Determine which arm based on position
        # First half of chromosome = L arm (positions closer to 0)
        # Second half of chromosome = R arm (positions closer to chr_size)
        midpoint = chr_size / 2

        if position < midpoint:
            arm = 'L'
        else:
            arm = 'R'

        chr_end = f"chr{chr_num}{arm}"
        yprime_counts[chr_end] += 1
        yprime_positions[chr_end].append({
            'position': int(position),
            'identity': hit['pident'],
            'strand': hit['strand']
        })

    return dict(yprime_counts), dict(yprime_positions)


def count_yprimes_from_labeled_regions(labeled_tsv):
    """
    Alternative: Count Y primes from already-labeled TSV file.

    This uses the output of label_pretelomeric_regions.py to get counts.

    Args:
        labeled_tsv: Path to labeled regions TSV

    Returns:
        Dict mapping chr_end to Y prime count
    """
    df = pd.read_csv(labeled_tsv, sep='\t')

    # Filter for Y primes only (column may be 'type' or 'region_type')
    type_col = 'region_type' if 'region_type' in df.columns else 'type'
    yprimes = df[df[type_col] == 'y_prime']

    # Count per chr_end
    counts = yprimes.groupby('chr_end').size().to_dict()

    return counts


def generate_output(yprime_counts, output_file, strain_id="unknown"):
    """
    Generate output files with expected Y prime counts.

    Args:
        yprime_counts: Dict mapping chr_end to count
        output_file: Path to output TSV
        strain_id: Strain identifier for output
    """
    # Generate all expected chr ends
    all_chr_ends = [f"chr{i}{arm}" for i in range(1, 17) for arm in ['L', 'R']]

    rows = []
    total = 0
    for chr_end in all_chr_ends:
        count = yprime_counts.get(chr_end, 0)
        rows.append({
            'chr_end': chr_end,
            'anchor_name': f"{chr_end}_anchor",
            'expected_y_primes': count
        })
        total += count

    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)

    print(f"\nExpected Y prime counts written to: {output_file}")
    print(f"Total Y primes across all chr ends: {total}")

    return df


def generate_python_dict(yprime_counts, strain_id):
    """
    Generate Python dict code for use in y_prime_analysis.py.

    Args:
        yprime_counts: Dict mapping chr_end to count
        strain_id: Strain identifier
    """
    print(f"\n# Python dict for strain {strain_id} (paste into y_prime_analysis.py):")
    print(f"elif '{strain_id}' in base_name:")
    print("    y_prime_numbers = {")

    for i in range(1, 17):
        l_count = yprime_counts.get(f"chr{i}L", 0)
        r_count = yprime_counts.get(f"chr{i}R", 0)
        comma = "," if i < 16 else ""
        print(f"        'chr{i}L_anchor': {l_count}, 'chr{i}R_anchor': {r_count}{comma}")

    print("    }")


def main():
    parser = argparse.ArgumentParser(
        description='Count expected Y primes per chromosome end in a reference genome'
    )

    parser.add_argument(
        '--reference',
        required=True,
        help='Reference FASTA file'
    )
    parser.add_argument(
        '--probe',
        default='references/probe.fasta',
        help='Y prime probe FASTA file (default: references/probe.fasta)'
    )
    parser.add_argument(
        '--labeled-tsv',
        help='Alternative: Use labeled regions TSV from label_pretelomeric_regions.py'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output TSV file for Y prime counts'
    )
    parser.add_argument(
        '--strain',
        default='unknown',
        help='Strain identifier (for generating Python dict code)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for BLAST (default: 4)'
    )
    parser.add_argument(
        '--min-identity',
        type=float,
        default=90.0,
        help='Minimum percent identity for probe matches (default: 90)'
    )

    args = parser.parse_args()

    print("=" * 80)
    print("Reference Y Prime Counter")
    print("=" * 80)

    if args.labeled_tsv:
        # Use pre-labeled TSV file
        print(f"\nUsing labeled regions from: {args.labeled_tsv}")
        yprime_counts = count_yprimes_from_labeled_regions(args.labeled_tsv)
        yprime_positions = {}  # Not available from TSV

    else:
        # Run BLAST with probe
        print(f"\nReference: {args.reference}")
        print(f"Probe: {args.probe}")
        print(f"Min identity: {args.min_identity}%")

        # Create temp file for BLAST output
        blast_output = args.output.replace('.tsv', '_probe_blast.txt')

        # Run BLAST
        print("\nStep 1: Running BLAST for Y prime probe...")
        run_blast_probe(args.probe, args.reference, blast_output,
                       min_identity=args.min_identity, num_threads=args.threads)

        # Parse results
        print("\nStep 2: Parsing BLAST results...")
        probe_hits = parse_blast_results(blast_output, min_length=80)
        print(f"Found {len(probe_hits)} probe hits")

        # Get chromosome sizes
        print("\nStep 3: Getting chromosome sizes...")
        chr_sizes = get_chromosome_sizes(args.reference)
        print(f"Found {len(chr_sizes)} chromosomes")

        # Assign to chr ends
        print("\nStep 4: Assigning Y primes to chromosome ends...")
        yprime_counts, yprime_positions = assign_yprimes_to_chr_ends(probe_hits, chr_sizes)

    # Generate output
    print("\nStep 5: Generating output...")
    df = generate_output(yprime_counts, args.output, args.strain)

    # Print summary
    print("\n" + "=" * 80)
    print("EXPECTED Y PRIMES PER CHROMOSOME END")
    print("=" * 80)

    for i in range(1, 17):
        l_end = f"chr{i}L"
        r_end = f"chr{i}R"
        l_count = yprime_counts.get(l_end, 0)
        r_count = yprime_counts.get(r_end, 0)
        print(f"  chr{i:2d}:  L={l_count}  R={r_count}")

    # Generate Python dict code
    generate_python_dict(yprime_counts, args.strain)

    # Show positions if available
    if yprime_positions:
        print("\n" + "=" * 80)
        print("Y PRIME POSITIONS (for verification)")
        print("=" * 80)
        for chr_end in sorted(yprime_positions.keys()):
            positions = yprime_positions[chr_end]
            print(f"\n{chr_end} ({len(positions)} Y primes):")
            for pos in sorted(positions, key=lambda x: x['position']):
                print(f"  Position: {pos['position']:>12,} | Identity: {pos['identity']:.1f}% | Strand: {pos['strand']}")

    print("\n" + "=" * 80)
    print("Done!")
    print("=" * 80)


if __name__ == '__main__':
    main()
