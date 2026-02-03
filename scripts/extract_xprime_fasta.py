#!/usr/bin/env python3
"""
Extract X prime sequences from a reference FASTA based on BED annotations.

X prime = x_core_element + x_variable_element (combined, must be adjacent)
"""

import argparse
import sys
from collections import defaultdict


def parse_fasta(fasta_path):
    """Parse a FASTA file and return dict of {seq_id: sequence}."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def parse_bed(bed_path):
    """
    Parse BED file and extract x_core_element and x_variable_element entries.

    Returns dict of {chr_end: {'x_core': (chr, start, end, strand), 'x_variable': (chr, start, end, strand)}}
    """
    x_elements = defaultdict(dict)

    with open(bed_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            strand = parts[4]

            # Extract chr end from name (e.g., chr1L_x_core_element -> chr1L)
            if '_x_core_element' in name:
                chr_end = name.replace('_x_core_element', '')
                x_elements[chr_end]['x_core'] = (chrom, start, end, strand)
            elif '_x_variable_element' in name:
                chr_end = name.replace('_x_variable_element', '')
                x_elements[chr_end]['x_variable'] = (chrom, start, end, strand)

    return x_elements


def extract_xprimes(fasta_sequences, x_elements):
    """
    Extract X prime sequences by combining x_core and x_variable elements.
    If x_variable is missing, extracts only the x_core with a warning.

    Returns list of (chr_end, sequence, info_dict), list of warnings, and list of errors.
    """
    xprime_sequences = []
    warnings = []
    errors = []

    # Expected chr ends (chr1L, chr1R, chr2L, chr2R, ..., chr16L, chr16R)
    expected_chr_ends = []
    for i in range(1, 17):
        expected_chr_ends.append(f"chr{i}L")
        expected_chr_ends.append(f"chr{i}R")

    for chr_end in expected_chr_ends:
        if chr_end not in x_elements:
            errors.append(f"ERROR: {chr_end} - No x_core_element or x_variable_element found")
            continue

        elements = x_elements[chr_end]

        # Check if x_core_element exists
        if 'x_core' not in elements:
            errors.append(f"ERROR: {chr_end} - Missing x_core_element")
            continue

        core_chrom, core_start, core_end, core_strand = elements['x_core']

        # Check if x_variable_element exists
        if 'x_variable' not in elements:
            # Extract just the core element with a warning
            warnings.append(f"WARNING: {chr_end} - Missing x_variable_element, extracting only x_core_element ({core_chrom}:{core_start}-{core_end})")

            # Check chromosome exists in FASTA
            if core_chrom not in fasta_sequences:
                errors.append(f"ERROR: {chr_end} - Chromosome {core_chrom} not found in FASTA")
                continue

            chrom_seq = fasta_sequences[core_chrom]

            if core_end > len(chrom_seq):
                errors.append(f"ERROR: {chr_end} - Coordinates {core_start}-{core_end} exceed chromosome length {len(chrom_seq)}")
                continue

            xprime_seq = chrom_seq[core_start:core_end]

            # Reverse complement if on minus strand
            if core_strand == '-':
                xprime_seq = reverse_complement(xprime_seq)

            info = {
                'chr': core_chrom,
                'start': core_start,
                'end': core_end,
                'strand': core_strand,
                'length': core_end - core_start,
                'core_coords': f"{core_start}-{core_end}",
                'var_coords': 'MISSING',
                'core_only': True,
                'core_len': core_end - core_start,
                'var_len': 0
            }

            xprime_sequences.append((chr_end, xprime_seq, info))
            continue

        var_chrom, var_start, var_end, var_strand = elements['x_variable']

        # Check chromosomes match
        if core_chrom != var_chrom:
            errors.append(f"ERROR: {chr_end} - x_core_element on {core_chrom} but x_variable_element on {var_chrom}")
            continue

        # Check strands match
        if core_strand != var_strand:
            errors.append(f"ERROR: {chr_end} - Strand mismatch: x_core={core_strand}, x_variable={var_strand}")
            continue

        # Check if they are adjacent (one ends where other begins, allowing for BED 0-based coordinates)
        # BED format is 0-based, half-open: [start, end)
        # So adjacent regions would have one's end == other's start

        # Determine order based on strand and position
        # For minus strand: variable comes before core (lower coordinates)
        # For plus strand: core comes before variable (lower coordinates)

        is_adjacent = False
        if core_end == var_start:
            # core is before variable
            combined_start = core_start
            combined_end = var_end
            is_adjacent = True
        elif var_end == core_start:
            # variable is before core
            combined_start = var_start
            combined_end = core_end
            is_adjacent = True
        elif core_end + 1 == var_start:
            # core is before variable with 1bp gap (BED coordinates can be off by 1)
            combined_start = core_start
            combined_end = var_end
            is_adjacent = True
        elif var_end + 1 == core_start:
            # variable is before core with 1bp gap
            combined_start = var_start
            combined_end = core_end
            is_adjacent = True

        if not is_adjacent:
            gap = min(abs(core_start - var_end), abs(var_start - core_end))
            errors.append(f"ERROR: {chr_end} - x_core_element and x_variable_element are NOT adjacent! "
                         f"x_core: {core_chrom}:{core_start}-{core_end}, "
                         f"x_variable: {var_chrom}:{var_start}-{var_end}, "
                         f"gap/overlap: {gap}bp")
            continue

        # Check chromosome exists in FASTA
        if core_chrom not in fasta_sequences:
            errors.append(f"ERROR: {chr_end} - Chromosome {core_chrom} not found in FASTA")
            continue

        # Extract sequence (BED is 0-based, Python slicing is 0-based)
        chrom_seq = fasta_sequences[core_chrom]

        if combined_end > len(chrom_seq):
            errors.append(f"ERROR: {chr_end} - Coordinates {combined_start}-{combined_end} exceed chromosome length {len(chrom_seq)}")
            continue

        xprime_seq = chrom_seq[combined_start:combined_end]

        # Reverse complement if on minus strand
        if core_strand == '-':
            xprime_seq = reverse_complement(xprime_seq)

        info = {
            'chr': core_chrom,
            'start': combined_start,
            'end': combined_end,
            'strand': core_strand,
            'length': combined_end - combined_start,
            'core_coords': f"{core_start}-{core_end}",
            'var_coords': f"{var_start}-{var_end}",
            'core_only': False,
            'core_len': core_end - core_start,
            'var_len': var_end - var_start
        }

        xprime_sequences.append((chr_end, xprime_seq, info))

    return xprime_sequences, warnings, errors


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def write_fasta(xprime_sequences, output_path, line_width=80):
    """Write X prime sequences to FASTA file."""
    with open(output_path, 'w') as f:
        for chr_end, seq, info in xprime_sequences:
            # Create header with useful info
            # Include core and variable lengths for downstream splitting
            core_only_tag = " [core_only]" if info.get('core_only', False) else ""

            # Add core_len and var_len to header for downstream processing
            core_len = info.get('core_len', info['length'])  # Default to full length if core_only
            var_len = info.get('var_len', 0)

            header = f">{chr_end}_xprime {info['chr']}:{info['start']}-{info['end']}({info['strand']}) len={info['length']} core_len={core_len} var_len={var_len}{core_only_tag}"
            f.write(header + '\n')

            # Write sequence with line breaks
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Extract X prime sequences from reference FASTA based on BED annotations'
    )
    parser.add_argument(
        '--reference', '-r',
        required=True,
        help='Reference FASTA file (e.g., 6991.fasta)'
    )
    parser.add_argument(
        '--bed', '-b',
        required=True,
        help='BED file with x_core_element and x_variable_element annotations'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output FASTA file for X prime sequences'
    )

    args = parser.parse_args()

    print("=" * 70)
    print("X Prime Extraction")
    print("=" * 70)
    print(f"Reference FASTA: {args.reference}")
    print(f"BED file: {args.bed}")
    print(f"Output: {args.output}")
    print()

    # Parse input files
    print("Parsing reference FASTA...")
    fasta_sequences = parse_fasta(args.reference)
    print(f"  Found {len(fasta_sequences)} chromosomes")

    print("Parsing BED file...")
    x_elements = parse_bed(args.bed)
    print(f"  Found x_core/x_variable elements for {len(x_elements)} chromosome ends")
    print()

    # Extract X primes
    print("Extracting X prime sequences...")
    xprime_sequences, warnings, errors = extract_xprimes(fasta_sequences, x_elements)

    # Report warnings
    if warnings:
        print()
        print("=" * 70)
        print(f"WARNINGS: {len(warnings)}")
        print("=" * 70)
        for warning in warnings:
            print(warning)
        print()

    # Report errors
    if errors:
        print()
        print("=" * 70)
        print(f"ERRORS FOUND: {len(errors)}")
        print("=" * 70)
        for error in errors:
            print(error)
        print()

    # Write output
    if xprime_sequences:
        write_fasta(xprime_sequences, args.output)
        print(f"Successfully extracted {len(xprime_sequences)} X prime sequences")
        print(f"Output written to: {args.output}")

        # Count core-only extractions
        core_only_count = sum(1 for _, _, info in xprime_sequences if info.get('core_only', False))
        complete_count = len(xprime_sequences) - core_only_count

        print(f"  - Complete (core + variable): {complete_count}")
        print(f"  - Core only (variable missing): {core_only_count}")

        # Print summary
        print()
        print("Summary:")
        print("-" * 85)
        print(f"{'Chr End':<10} {'Chromosome':<8} {'Start':>12} {'End':>12} {'Length':>8} {'Strand':<6} {'Note':<15}")
        print("-" * 85)
        for chr_end, seq, info in xprime_sequences:
            note = "core only" if info.get('core_only', False) else ""
            print(f"{chr_end:<10} {info['chr']:<8} {info['start']:>12,} {info['end']:>12,} {info['length']:>8} {info['strand']:<6} {note:<15}")
    else:
        print("No X prime sequences could be extracted!")
        sys.exit(1)

    # Exit with error code if there were errors (but not for warnings)
    if errors:
        print()
        print(f"ERROR: {len(errors)} chromosome ends had errors and were skipped")
        sys.exit(1)


if __name__ == '__main__':
    main()
