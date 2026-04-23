#!/usr/bin/env python3
"""
Command-line wrapper for telomere pipeline functions.
"""

import argparse
import sys
from subtelomere_reference_pipeline_utils import (
    get_75th_percentile_reads,
    extract_selected_reads,
    extend_reference_multi
)


def main():
    parser = argparse.ArgumentParser(description='Telomere pipeline utilities')
    subparsers = parser.add_subparsers(dest='command', help='Command to run')

    # Select reads command
    select_parser = subparsers.add_parser('select_reads',
                                          help='Select 75th percentile reads')
    select_parser.add_argument('--input-tsv', required=True,
                              help='Input TSV file')
    select_parser.add_argument('--output-dir', required=True,
                              help='Output directory')
    select_parser.add_argument('--output-file', required=True,
                              help='Output file for read IDs')

    # Extract reads command
    extract_parser = subparsers.add_parser('extract_reads',
                                          help='Extract selected reads from FASTQ')
    extract_parser.add_argument('--reads-fastq', required=True,
                               help='Input FASTQ file')
    extract_parser.add_argument('--read-ids-file', required=True,
                               help='File with read IDs')
    extract_parser.add_argument('--output-fastq', required=True,
                               help='Output FASTQ file')

    # Extend reference command
    extend_parser = subparsers.add_parser('extend_reference',
                                         help='Extend reference using soft-clipped bases')
    extend_parser.add_argument('--bamfile', required=True,
                              help='Input BAM file')
    extend_parser.add_argument('--reference', required=True,
                              help='Reference FASTA file')
    extend_parser.add_argument('--read-ids-file', required=True,
                              help='File with read IDs')
    extend_parser.add_argument('--output-fasta', required=True,
                              help='Output extended FASTA file')
    extend_parser.add_argument('--trim', type=int, default=0,
                              help='Number of bp to trim from extensions (default: 0)')

    args = parser.parse_args()

    if args.command == 'select_reads':
        get_75th_percentile_reads(
            args.input_tsv,
            args.output_dir,
            args.output_file
        )
    elif args.command == 'extract_reads':
        extract_selected_reads(
            args.reads_fastq,
            args.read_ids_file,
            args.output_fastq
        )
    elif args.command == 'extend_reference':
        extend_reference_multi(
            args.bamfile,
            args.reference,
            args.read_ids_file,
            args.output_fasta,
            trim=args.trim
        )
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()

