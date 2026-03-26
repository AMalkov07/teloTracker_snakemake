#!/usr/bin/env python3
"""
Make databases for each y_prime group and each pair of y_prime groups for x element sequences.

Usage:
    python make_databases_of_pairings_for_x_elements.py <strain_id> <y_prime_fasta> <x_element_fasta> <output_base_dir>

Example:
    python make_databases_of_pairings_for_x_elements.py 7302 \
        references/extracted_yprimes_7302.fasta \
        references/7302_whole_x_regions_sequences.fasta \
        references/pairings_for_x_element_ends/
"""

import sys
import os
from Bio import SeqIO
from itertools import combinations


def parse_y_prime_ids(y_prime_fasta):
    """Parse Y' prime FASTA to get Y' IDs and their associated chromosome ends."""
    y_prime_ids = {}

    with open(y_prime_fasta, 'r') as db_file:
        for line in db_file:
            if line.startswith(">"):
                try:
                    similar_y_prime_chr_ends = line.split("Y_Prime_")[1]
                    similar_y_prime_chr_ends = similar_y_prime_chr_ends.split("#")[0]
                    similar_y_prime_chr_ends = similar_y_prime_chr_ends.split(";")

                    y_prime_id = line.split("/")[2]
                    y_prime_id = y_prime_id.split("_")[0]

                    chr_end_group = [y_prime_chr_end.strip('_,0123456789') for y_prime_chr_end in similar_y_prime_chr_ends]
                    chr_end_group = list(set(chr_end_group))

                    if y_prime_id not in y_prime_ids:
                        y_prime_ids[y_prime_id] = chr_end_group
                    else:
                        for chr_end in chr_end_group:
                            if chr_end not in y_prime_ids[y_prime_id]:
                                y_prime_ids[y_prime_id].append(chr_end)
                except (IndexError, ValueError) as e:
                    print(f"  Warning: Could not parse header: {line.strip()}")
                    continue

    return dict(sorted(y_prime_ids.items()))


def get_y_prime_0_ends(strain_id):
    """Get the chromosome ends that have 0 Y' primes for a given strain."""
    # These are the ends without Y' primes - same for most strains
    y_prime_0_ends = ['chr1L', 'chr1R', 'chr2R', 'chr3L', 'chr3R', 'chr4L', 'chr6R',
                      'chr7L', 'chr9R', 'chr10R', 'chr11L', 'chr11R', 'chr13R', 'chr15L']
    return y_prime_0_ends


def main():
    if len(sys.argv) < 5:
        print("Usage: python make_databases_of_pairings_for_x_elements.py <strain_id> <y_prime_fasta> <x_element_fasta> <output_base_dir>")
        print("\nExample:")
        print("  python make_databases_of_pairings_for_x_elements.py 7302 \\")
        print("      references/extracted_yprimes_7302.fasta \\")
        print("      references/7302_whole_x_regions_sequences.fasta \\")
        print("      references/pairings_for_x_element_ends/")
        sys.exit(1)

    strain_id = sys.argv[1]
    y_prime_fasta = sys.argv[2]
    x_element_fasta = sys.argv[3]
    output_base_dir = sys.argv[4]

    print(f'Running {strain_id}...')
    print(f'Y\' prime FASTA: {y_prime_fasta}')
    print(f'X element FASTA: {x_element_fasta}')
    print(f'Output directory: {output_base_dir}')

    # Create output directories
    y_prime_group_solo_dir = os.path.join(output_base_dir, f'{strain_id}_individual_y_prime_groups')
    output_pairing_dir = os.path.join(output_base_dir, f'{strain_id}_pairings')

    os.makedirs(y_prime_group_solo_dir, exist_ok=True)
    os.makedirs(output_pairing_dir, exist_ok=True)

    # Parse Y' prime IDs
    print("\nParsing Y' prime IDs...")
    y_prime_ids = parse_y_prime_ids(y_prime_fasta)
    print(f"Found {len(y_prime_ids)} Y' prime groups:")
    for k, v in y_prime_ids.items():
        print(f"  {k}: {v}")

    # Create individual Y' prime group databases
    print("\nCreating individual Y' prime group databases...")
    for chr_group_id in y_prime_ids.keys():
        print(f'  Creating {strain_id}: {chr_group_id}')
        output_file = os.path.join(y_prime_group_solo_dir, f'{strain_id}_y_prime_group_{chr_group_id}.fasta')

        with open(x_element_fasta, 'r') as all_x_elements_file, open(output_file, 'w') as out_file:
            for record in SeqIO.parse(all_x_elements_file, 'fasta'):
                header = record.description
                chr_end_in_header = header.split("_")[0]
                if chr_end_in_header in y_prime_ids[chr_group_id]:
                    SeqIO.write(record, out_file, 'fasta')

    # Create paired databases for Y' prime group combinations
    print("\nCreating paired Y' prime group databases...")
    unique_pairs = sorted(combinations(y_prime_ids.keys(), 2))
    paired_y_prime_ids = {f'{pair[0]}_and_{pair[1]}': list(set(y_prime_ids[pair[0]] + y_prime_ids[pair[1]])) for pair in unique_pairs}
    paired_y_prime_ids = dict(sorted(paired_y_prime_ids.items()))

    for chr_group_id_pair in paired_y_prime_ids.keys():
        print(f'  Creating {strain_id}: {chr_group_id_pair}')
        output_file = os.path.join(output_pairing_dir, f'{strain_id}_paired_{chr_group_id_pair}.fasta')

        with open(x_element_fasta, 'r') as all_x_elements_file, open(output_file, 'w') as out_file:
            for record in SeqIO.parse(all_x_elements_file, 'fasta'):
                header = record.description
                chr_end_in_header = header.split("_")[0]
                if chr_end_in_header in paired_y_prime_ids[chr_group_id_pair]:
                    SeqIO.write(record, out_file, 'fasta')

    # Create paired databases for 0 Y' ends with each Y' prime group
    print("\nCreating 0 Y' ends paired databases...")
    y_prime_0_ends = get_y_prime_0_ends(strain_id)

    all_paired_y_prime_ids_with_0 = {}
    for chr_end in y_prime_0_ends:
        for y_prime_id in y_prime_ids:
            all_paired_y_prime_ids_with_0[f'{chr_end}_and_{y_prime_id}'] = y_prime_ids[y_prime_id] + [chr_end]
    all_paired_y_prime_ids_with_0 = dict(sorted(all_paired_y_prime_ids_with_0.items()))

    for chr_end_and_group_id_pair in all_paired_y_prime_ids_with_0.keys():
        print(f'  Creating {strain_id}: {chr_end_and_group_id_pair}')
        output_file = os.path.join(output_pairing_dir, f'{strain_id}_paired_{chr_end_and_group_id_pair}.fasta')

        with open(x_element_fasta, 'r') as all_x_elements_file, open(output_file, 'w') as out_file:
            for record in SeqIO.parse(all_x_elements_file, 'fasta'):
                header = record.description
                chr_end_in_header = header.split("_")[0]
                if chr_end_in_header in all_paired_y_prime_ids_with_0[chr_end_and_group_id_pair]:
                    SeqIO.write(record, out_file, 'fasta')

    print(f"\nDone! Created pairing databases in {output_pairing_dir}")


if __name__ == "__main__":
    main()
