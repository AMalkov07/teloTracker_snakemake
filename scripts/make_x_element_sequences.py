#!/usr/bin/env python3
"""
Generate x element sequences from a BED file and reference genome.

This script extracts x_core_element and x_variable_element regions from a BED file,
combines them for each chromosome end, and outputs a FASTA file for use with
RepeatMasker pairing analysis.

Usage:
    python make_x_element_sequences.py <strain_id> <bed_file> <reference_fasta> <output_dir>

Example:
    python make_x_element_sequences.py 7302 \
        labeling_test/new_ref_7302_bed/pretelomeric_regions_7302_simp.bed \
        labeling_test/new_ref_7302_bed/assembly_7302_dorado_reference.fasta \
        references/
"""

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict


def parse_bed_for_x_elements(bed_file):
    """
    Parse BED file and extract x_core_element and x_variable_element regions.

    Returns a dict of chr_end -> list of regions
    """
    x_regions = defaultdict(list)

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

            # Only process x element regions
            if 'x_core_element' in name or 'x_variable_element' in name:
                # Extract chr end name (e.g., "chr1L" from "chr1L_x_core_element")
                chr_end = name.split('_x_')[0]

                x_regions[chr_end].append({
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'strand': strand,
                    'type': 'core' if 'core' in name else 'variable'
                })

    return x_regions


def load_reference_genome(fasta_file):
    """Load reference genome into a dictionary."""
    genome = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
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
            if base_name.startswith('chr'):
                genome[base_name[3:]] = str(record.seq)

    return genome


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(sequence).reverse_complement())


def main():
    if len(sys.argv) < 5:
        print("Usage: python make_x_element_sequences.py <strain_id> <bed_file> <reference_fasta> <output_dir>")
        print("\nExample:")
        print("  python make_x_element_sequences.py 7302 \\")
        print("      labeling_test/new_ref_7302_bed/pretelomeric_regions_7302_simp.bed \\")
        print("      labeling_test/new_ref_7302_bed/assembly_7302_dorado_reference.fasta \\")
        print("      references/")
        sys.exit(1)

    strain_id = sys.argv[1]
    bed_file = sys.argv[2]
    reference_fasta = sys.argv[3]
    output_dir = sys.argv[4]

    print(f"Processing strain {strain_id}...")
    print(f"BED file: {bed_file}")
    print(f"Reference: {reference_fasta}")
    print(f"Output directory: {output_dir}")

    # Parse BED file for x element regions
    print("\nParsing BED file for x element regions...")
    x_regions = parse_bed_for_x_elements(bed_file)
    print(f"Found x elements for {len(x_regions)} chromosome ends")

    # Load reference genome
    print("\nLoading reference genome...")
    genome = load_reference_genome(reference_fasta)

    # Output file
    output_file = os.path.join(output_dir, f"{strain_id}_whole_x_regions_sequences.fasta")

    print(f"\nGenerating x element sequences...")
    total_sequences = 0

    with open(output_file, 'w') as out_f:
        for chr_end in sorted(x_regions.keys()):
            regions = x_regions[chr_end]

            # Sort regions by start position
            regions.sort(key=lambda r: r['start'])

            # Get chromosome
            chrom = regions[0]['chr']
            if chrom not in genome:
                print(f"  Warning: Chromosome {chrom} not found in reference, skipping {chr_end}")
                continue

            # Get the combined region (min start to max end)
            min_start = min(r['start'] for r in regions)
            max_end = max(r['end'] for r in regions)
            strand = regions[0]['strand']

            # Extract sequence
            sequence = genome[chrom][min_start:max_end]

            # Reverse complement if on negative strand
            if strand == '-':
                sequence = reverse_complement(sequence)

            # Write to output with x_ends suffix for RepeatMasker matching
            header = f">{chr_end}_x_ends"
            out_f.write(f"{header}\n{sequence}\n")
            total_sequences += 1

            print(f"  {chr_end}: {len(sequence)} bp (from {len(regions)} regions)")

    print(f"\nDone! Wrote {total_sequences} x element sequences to {output_file}")

    return output_file


if __name__ == "__main__":
    main()
