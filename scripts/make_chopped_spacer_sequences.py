#!/usr/bin/env python3
"""
Generate chopped spacer sequences from a BED file and reference genome.

This script extracts "space_between_anchor" regions from a BED file,
retrieves the sequences from a reference genome, and chops them into
250bp chunks for use with RepeatMasker pairing analysis.

Usage:
    python make_chopped_spacer_sequences.py <strain_id> <bed_file> <reference_fasta> <output_dir>

Example:
    python make_chopped_spacer_sequences.py 7302 \
        labeling_test/new_ref_7302_bed/pretelomeric_regions_7302_simp.bed \
        references/7302_reference.fasta \
        references/
"""

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq


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


def main():
    if len(sys.argv) < 5:
        print("Usage: python make_chopped_spacer_sequences.py <strain_id> <bed_file> <reference_fasta> <output_dir>")
        print("\nExample:")
        print("  python make_chopped_spacer_sequences.py 7302 \\")
        print("      labeling_test/new_ref_7302_bed/pretelomeric_regions_7302_simp.bed \\")
        print("      references/7302_reference.fasta \\")
        print("      references/")
        sys.exit(1)

    strain_id = sys.argv[1]
    bed_file = sys.argv[2]
    reference_fasta = sys.argv[3]
    output_dir = sys.argv[4]

    # Chunk size in bp
    chunk_size = 250

    print(f"Processing strain {strain_id}...")
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
    print(f"Loaded {len(genome) // 2} chromosomes")  # Divide by 2 because we store with/without 'chr' prefix

    # Output file
    output_file = os.path.join(output_dir, f"{strain_id}_50kb_chopped_up_spacer_sequences.fasta")

    print(f"\nGenerating chopped spacer sequences...")
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


if __name__ == "__main__":
    main()
