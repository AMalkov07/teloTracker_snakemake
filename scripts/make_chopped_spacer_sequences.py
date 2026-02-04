#!/usr/bin/env python3
"""
Generate chopped spacer sequences from a BED file and reference genome.

This script extracts "space_between_anchor" regions from a BED file,
retrieves the sequences from a reference genome, and chops them into
250bp chunks for use with RepeatMasker pairing analysis.

Two methods are available:
1. Variable-size (default): Uses actual spacer region sizes from BED file
2. Fixed 50kb (--fixed-50kb): Uses a fixed 50kb window from each telomere end

Usage:
    # Variable-size method (default)
    python make_chopped_spacer_sequences.py <strain_id> <bed_file> <reference_fasta> <output_dir>

    # Fixed 50kb method
    python make_chopped_spacer_sequences.py <strain_id> <bed_file> <reference_fasta> <output_dir> --fixed-50kb

Example:
    python make_chopped_spacer_sequences.py 7302 \
        labeling_test/new_ref_7302_bed/pretelomeric_regions_7302_simp.bed \
        references/7302_reference.fasta \
        references/

    # With fixed 50kb window
    python make_chopped_spacer_sequences.py 7302 \
        labeling_test/new_ref_7302_bed/pretelomeric_regions_7302_simp.bed \
        references/7302_reference.fasta \
        references/ --fixed-50kb
"""

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq


# Fixed window size for the old method (50kb)
FIXED_WINDOW_SIZE = 50000


def parse_bed_for_spacers(bed_file):
    """
    Parse BED file and extract space_between_anchor regions.

    Returns a list of dicts with: chr, start, end, name, strand, chr_end
    """
    spacer_regions = []

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue

            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            strand = fields[4] if len(fields) > 4 else '+'

            # Only process space_between_anchor regions
            if 'space_between_anchor' in name:
                # Extract chr end name (e.g., "chr1L" from "chr1L_space_between_anchor")
                chr_end = name.replace('_space_between_anchor', '')

                spacer_regions.append({
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'strand': strand,
                    'chr_end': chr_end
                })

    return spacer_regions


def parse_bed_for_telomere_ends(bed_file):
    """
    Parse BED file to get telomere and spacer positions for each chromosome end.

    This is used for the fixed 50kb method to determine where to start
    extracting sequences from. The extraction starts from near the spacer/X element
    boundary and goes OUTWARD toward the telomere (capturing Y' primes).

    Returns a list of dicts with: chr, spacer_start, spacer_end, chr_end, strand
    Only returns entries for actual chromosome ends (chr#L or chr#R), not ITS regions.
    """
    # First pass: collect all features by chromosome end
    features_by_chr_end = {}

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue

            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            strand = fields[4] if len(fields) > 4 else '+'

            # Skip ITS regions - they're not chromosome ends
            if name.startswith('ITS_'):
                continue

            # Extract chr end name from feature name
            # Only accept valid chromosome end patterns (chr#L or chr#R)
            chr_end = None
            for pattern in ['_Telomere', '_telomere', '_space_between', '_x_core', '_x_variable', '_anchor', '_Y_Prime']:
                if pattern in name:
                    potential_chr_end = name.split(pattern)[0]
                    # Validate it's a real chromosome end (chr#L or chr#R format)
                    if potential_chr_end.startswith('chr') and (potential_chr_end.endswith('L') or potential_chr_end.endswith('R')):
                        chr_end = potential_chr_end
                    break

            if chr_end is None:
                continue

            if chr_end not in features_by_chr_end:
                features_by_chr_end[chr_end] = {
                    'chr': chrom,
                    'chr_end': chr_end,
                    'strand': strand,
                    'spacer_start': None,
                    'spacer_end': None,
                    'telomere_start': None,
                    'telomere_end': None
                }

            # Track spacer positions
            if 'space_between_anchor' in name:
                features_by_chr_end[chr_end]['spacer_start'] = start
                features_by_chr_end[chr_end]['spacer_end'] = end

            # Track telomere positions
            if 'Telomere_Repeat' in name or 'telomere' in name.lower():
                features_by_chr_end[chr_end]['telomere_start'] = start
                features_by_chr_end[chr_end]['telomere_end'] = end

    return list(features_by_chr_end.values())


def load_reference_genome(fasta_file):
    """Load reference genome into a dictionary."""
    genome = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        # Handle different chromosome naming conventions
        chrom_name = record.id
        genome[chrom_name] = str(record.seq)

        # Also store without 'chr' prefix if it has one
        if chrom_name.startswith('chr'):
            genome[chrom_name[3:]] = str(record.seq)
        else:
            genome['chr' + chrom_name] = str(record.seq)

        # Handle '_extended' suffix (e.g., chr1_extended -> chr1)
        if '_extended' in chrom_name:
            base_name = chrom_name.replace('_extended', '')
            genome[base_name] = str(record.seq)
            # Also without 'chr' prefix
            if base_name.startswith('chr'):
                genome[base_name[3:]] = str(record.seq)

    return genome


def chop_sequence(sequence, chunk_size=250):
    """Chop a sequence into chunks of specified size."""
    chunks = []
    for i in range(0, len(sequence), chunk_size):
        chunk = sequence[i:i + chunk_size]
        if len(chunk) >= chunk_size // 2:  # Only include chunks that are at least half size
            chunks.append(chunk)
    return chunks


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(sequence).reverse_complement())


def generate_variable_size_spacers(strain_id, bed_file, reference_fasta, output_dir, chunk_size=250):
    """
    Generate spacer sequences using actual spacer region sizes from BED file.

    This is the NEW method that uses precise spacer boundaries.
    """
    print(f"Processing strain {strain_id} (variable-size method)...")
    print(f"BED file: {bed_file}")
    print(f"Reference: {reference_fasta}")
    print(f"Output directory: {output_dir}")

    # Parse BED file for spacer regions
    print("\nParsing BED file for space_between_anchor regions...")
    spacer_regions = parse_bed_for_spacers(bed_file)
    print(f"Found {len(spacer_regions)} spacer regions")

    # Load reference genome
    print("\nLoading reference genome...")
    genome = load_reference_genome(reference_fasta)
    print(f"Loaded {len(genome) // 2} chromosomes")

    # Output file
    output_file = os.path.join(output_dir, f"{strain_id}_50kb_chopped_up_spacer_sequences.fasta")

    print(f"\nGenerating chopped spacer sequences (variable-size)...")
    total_chunks = 0

    with open(output_file, 'w') as out_f:
        for region in spacer_regions:
            chrom = region['chr']
            start = region['start']
            end = region['end']
            chr_end = region['chr_end']
            strand = region['strand']

            # Get chromosome sequence
            if chrom not in genome:
                print(f"  Warning: Chromosome {chrom} not found in reference, skipping {chr_end}")
                continue

            # Extract sequence (BED is 0-based, half-open)
            sequence = genome[chrom][start:end]

            # Reverse complement if on negative strand
            if strand == '-':
                sequence = reverse_complement(sequence)

            # Chop into chunks
            chunks = chop_sequence(sequence, chunk_size)

            print(f"  {chr_end}: {len(sequence)} bp -> {len(chunks)} chunks")

            # Write chunks to output
            for i, chunk in enumerate(chunks, 1):
                header = f">{chr_end}_50kb_spacer_section_{i}_from_repeat_to_plus_50kb#{i}"
                out_f.write(f"{header}\n{chunk}\n")
                total_chunks += 1

    print(f"\nDone! Wrote {total_chunks} chunks to {output_file}")
    return output_file


def generate_fixed_50kb_spacers(strain_id, bed_file, reference_fasta, output_dir, chunk_size=250):
    """
    Generate spacer sequences using a fixed 50kb window from each chromosome end.

    This is the OLD method that extracts exactly 50kb (200 chunks) per chromosome end.

    The extraction starts from near the spacer/X element boundary and goes OUTWARD
    toward the telomere, capturing the Y' primes, ITS regions, and telomere sequences.
    This is useful for detecting recombination by comparing these regions between
    different chromosome ends.
    """
    print(f"Processing strain {strain_id} (fixed 50kb method)...")
    print(f"BED file: {bed_file}")
    print(f"Reference: {reference_fasta}")
    print(f"Output directory: {output_dir}")
    print(f"Window size: {FIXED_WINDOW_SIZE} bp")

    # Parse BED file for feature positions
    print("\nParsing BED file for chromosome end positions...")
    chr_end_features = parse_bed_for_telomere_ends(bed_file)
    print(f"Found {len(chr_end_features)} chromosome ends")

    # Load reference genome
    print("\nLoading reference genome...")
    genome = load_reference_genome(reference_fasta)
    print(f"Loaded {len(genome) // 2} chromosomes")

    # Output file
    output_file = os.path.join(output_dir, f"{strain_id}_50kb_chopped_up_spacer_sequences.fasta")

    print(f"\nGenerating chopped spacer sequences (fixed 50kb)...")
    print("Extraction direction: from spacer/X-element boundary OUTWARD toward telomere")
    total_chunks = 0

    with open(output_file, 'w') as out_f:
        for features in chr_end_features:
            chrom = features['chr']
            chr_end = features['chr_end']
            strand = features['strand']
            spacer_start = features['spacer_start']
            spacer_end = features['spacer_end']

            # Get chromosome sequence
            if chrom not in genome:
                print(f"  Warning: Chromosome {chrom} not found in reference, skipping {chr_end}")
                continue

            chrom_seq = genome[chrom]
            chrom_len = len(chrom_seq)

            # Determine extraction region based on chromosome arm
            # The goal is to extract 50kb going from the spacer boundary TOWARD the telomere
            # This captures the Y' primes, ITS regions, and telomere sequences

            if chr_end.endswith('L') or strand == '-':
                # L arm: spacer is after telomere (telomere at start of chromosome)
                # Extract 50kb starting from position 0 (telomere end)
                # This captures: telomere, Y' primes (if any), X elements, spacer, etc.
                start = 0
                end = min(FIXED_WINDOW_SIZE, chrom_len)

                sequence = chrom_seq[start:end]
                # Reverse complement for L arms to get sequence in consistent orientation
                # (telomere -> centromere direction)
                sequence = reverse_complement(sequence)
            else:
                # R arm: spacer is before telomere (telomere at end of chromosome)
                # Extract exactly 50kb ending at chromosome end
                # This captures the telomeric region including Y' primes, X elements, etc.
                end = chrom_len
                start = max(0, end - FIXED_WINDOW_SIZE)

                sequence = chrom_seq[start:end]

            # Chop into chunks
            chunks = chop_sequence(sequence, chunk_size)

            print(f"  {chr_end}: {len(sequence)} bp -> {len(chunks)} chunks (pos {start}-{end})")

            # Write chunks to output
            for i, chunk in enumerate(chunks, 1):
                header = f">{chr_end}_50kb_spacer_section_{i}_from_repeat_to_plus_50kb#{i}"
                out_f.write(f"{header}\n{chunk}\n")
                total_chunks += 1

    print(f"\nDone! Wrote {total_chunks} chunks to {output_file}")
    return output_file


def main():
    if len(sys.argv) < 5:
        print("Usage: python make_chopped_spacer_sequences.py <strain_id> <bed_file> <reference_fasta> <output_dir> [--fixed-50kb]")
        print("\nMethods:")
        print("  Default (no flag): Variable-size method using actual spacer regions from BED")
        print("  --fixed-50kb: Fixed 50kb window from each telomere end (old method)")
        print("\nExample:")
        print("  # Variable-size method (default)")
        print("  python make_chopped_spacer_sequences.py 7302 \\")
        print("      pretelomeric_regions_7302_simp.bed \\")
        print("      references/7302_reference.fasta \\")
        print("      references/")
        print("")
        print("  # Fixed 50kb method")
        print("  python make_chopped_spacer_sequences.py 7302 \\")
        print("      pretelomeric_regions_7302_simp.bed \\")
        print("      references/7302_reference.fasta \\")
        print("      references/ --fixed-50kb")
        sys.exit(1)

    strain_id = sys.argv[1]
    bed_file = sys.argv[2]
    reference_fasta = sys.argv[3]
    output_dir = sys.argv[4]

    # Check for --fixed-50kb flag
    use_fixed_50kb = '--fixed-50kb' in sys.argv

    # Chunk size in bp
    chunk_size = 250

    if use_fixed_50kb:
        print("=" * 60)
        print("Using FIXED 50kb method (old method)")
        print("=" * 60)
        output_file = generate_fixed_50kb_spacers(strain_id, bed_file, reference_fasta, output_dir, chunk_size)
    else:
        print("=" * 60)
        print("Using VARIABLE-SIZE method (new method)")
        print("=" * 60)
        output_file = generate_variable_size_spacers(strain_id, bed_file, reference_fasta, output_dir, chunk_size)

    return output_file


if __name__ == "__main__":
    main()
