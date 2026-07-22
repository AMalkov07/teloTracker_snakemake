#!/usr/bin/env python3
"""
Command-line wrapper for telomere pipeline functions.
"""

import argparse
import sys
from subtelomere_reference_pipeline_utils import (
    get_75th_percentile_reads,
    select_consensus_scaffold_reads,
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
    select_parser.add_argument('--all-matches-tsv', default=None,
                              help='all_matches_<base>_blasted_<anchor>.tsv. Supplying this '
                                   '(with --reads-fasta) switches on consensus selection: the '
                                   'scaffold must be corroborated by other near-percentile '
                                   'reads instead of being trusted on its own.')
    select_parser.add_argument('--reads-fasta', default=None,
                              help='<base>.fasta, with a .fai alongside it')
    select_parser.add_argument('--n-candidates', type=int, default=5,
                              help='reads to compare around the 75th percentile (default 5)')
    select_parser.add_argument('--min-agree', type=int, default=3,
                              help='how many must mutually agree (default 3)')
    select_parser.add_argument('--min-identity', type=float, default=90.0,
                              help='pident floor for two candidates to agree (default 90)')
    select_parser.add_argument('--min-coverage', type=float, default=0.90,
                              help='best-HSP coverage floor to agree (default 0.90)')
    select_parser.add_argument('--report-tsv', default=None,
                              help='optional per-end record of what was chosen and why')
    select_parser.add_argument('--threads', type=int, default=4)

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
    extend_parser.add_argument('--chr-arm-pairs-file', default=None,
                              help='<chr_end>\\t<read_id> pairs from select_reads. Binds each '
                                   'scaffold to its own arm: L reads may only supply a prefix, '
                                   'R reads only a suffix. Strongly recommended -- without it a '
                                   'read clipped on both sides can win the other arm\'s slot.')

    args = parser.parse_args()

    if args.command == 'select_reads':
        if args.all_matches_tsv and args.reads_fasta:
            select_consensus_scaffold_reads(
                args.input_tsv,
                args.all_matches_tsv,
                args.reads_fasta,
                args.output_dir,
                args.output_file,
                n_candidates=args.n_candidates,
                min_agree=args.min_agree,
                min_identity=args.min_identity,
                min_coverage=args.min_coverage,
                threads=args.threads,
                report_tsv=args.report_tsv
            )
        else:
            print('WARNING: --all-matches-tsv/--reads-fasta not supplied; falling back to '
                  'the single 75th-percentile read. That read is scaffolded with no check '
                  'that it is typical of its chromosome end.')
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
            trim=args.trim,
            chr_arm_pairs_file=args.chr_arm_pairs_file
        )
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()

