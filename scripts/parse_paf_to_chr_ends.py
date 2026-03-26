#!/usr/bin/env python3
"""
parse_paf_to_chr_ends.py

Parses a PAF file of assembly contigs mapped to S288C R64 reference,
identifies which contigs correspond to chromosome L/R ends, and outputs:
  1. A summary of how many chromosome ends are covered
  2. A FASTA with sequences labeled as chr1L, chr1R, chr2L, etc.

For contigs that span both arms of a small chromosome, the contig is split
at the position corresponding to the chromosome midpoint.

Usage:
    python scripts/parse_paf_to_chr_ends.py <paf> <assembly.fasta> <output.fasta>
    python scripts/parse_paf_to_chr_ends.py <paf> <assembly.fasta> <output.fasta> --threshold 100000 --min-mapq 30
"""

import sys
import argparse
from collections import defaultdict

# S288C R64 chromosome info: accession -> (chr_number, chr_length)
CHR_INFO = {
    'NC_001133.9':  (1,  230218),
    'NC_001134.8':  (2,  813184),
    'NC_001135.5':  (3,  316620),
    'NC_001136.10': (4,  1531933),
    'NC_001137.3':  (5,  576874),
    'NC_001138.5':  (6,  270161),
    'NC_001139.9':  (7,  1090940),
    'NC_001140.6':  (8,  562643),
    'NC_001141.2':  (9,  439888),
    'NC_001142.9':  (10, 745751),
    'NC_001143.9':  (11, 666816),
    'NC_001144.5':  (12, 1078177),
    'NC_001145.3':  (13, 924431),
    'NC_001146.8':  (14, 784333),
    'NC_001147.6':  (15, 1091291),
    'NC_001148.4':  (16, 948066),
}


def parse_paf(paf_file, min_mapq=30):
    """
    Parse PAF file and return the best alignment per contig.
    'Best' is defined as the longest alignment block among high-quality hits.
    Only nuclear chromosomes (in CHR_INFO) are considered.
    """
    best_hits = {}

    with open(paf_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            query_name   = parts[0]
            query_len    = int(parts[1])
            query_start  = int(parts[2])
            query_end    = int(parts[3])
            strand       = parts[4]
            target_name  = parts[5]
            target_start = int(parts[7])
            target_end   = int(parts[8])
            align_len    = int(parts[10])
            mapq         = int(parts[11])

            if mapq < min_mapq:
                continue
            if target_name not in CHR_INFO:
                continue

            record = {
                'query_name':   query_name,
                'query_len':    query_len,
                'query_start':  query_start,
                'query_end':    query_end,
                'strand':       strand,
                'target_name':  target_name,
                'target_start': target_start,
                'target_end':   target_end,
                'align_len':    align_len,
                'mapq':         mapq,
            }

            if query_name not in best_hits or align_len > best_hits[query_name]['align_len']:
                best_hits[query_name] = record

    return best_hits


def classify_contig(record, threshold):
    """
    Classify a contig's alignment as L arm, R arm, BOTH, or internal
    based on how close the alignment gets to the chromosome ends.
    """
    chr_num, chr_len = CHR_INFO[record['target_name']]
    near_L = record['target_start'] < threshold
    near_R = record['target_end'] > chr_len - threshold

    if near_L and near_R:
        return 'BOTH'
    elif near_L:
        return 'L'
    elif near_R:
        return 'R'
    else:
        return 'internal'


def get_split_position(record):
    """
    For a contig spanning both arms, compute the query (contig) position
    that corresponds to the chromosome midpoint using linear interpolation.
    Handles both + and - strand alignments.
    """
    chr_num, chr_len = CHR_INFO[record['target_name']]
    chr_midpoint = chr_len / 2

    t_start = record['target_start']
    t_end   = record['target_end']
    q_start = record['query_start']
    q_end   = record['query_end']

    fraction = (chr_midpoint - t_start) / (t_end - t_start)
    fraction = max(0.0, min(1.0, fraction))

    if record['strand'] == '+':
        split = int(q_start + fraction * (q_end - q_start))
    else:
        # On - strand the contig runs in reverse relative to the chromosome
        split = int(q_end - fraction * (q_end - q_start))

    return split


def reverse_complement(seq):
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]


def read_fasta(fasta_file):
    """Read a FASTA file into a dict {sequence_id: sequence}."""
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_name is not None:
        sequences[current_name] = ''.join(current_seq)

    return sequences


def write_fasta_seq(handle, label, sequence, line_width=80):
    handle.write(f'>{label}\n')
    for i in range(0, len(sequence), line_width):
        handle.write(sequence[i:i + line_width] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Parse PAF to extract and label chromosome end sequences'
    )
    parser.add_argument('paf',      help='PAF file from minimap2 (contigs vs S288C R64)')
    parser.add_argument('assembly', help='Assembly FASTA file (contigs)')
    parser.add_argument('output',   help='Output FASTA file with labeled chr end sequences')
    parser.add_argument('--threshold', type=int, default=100000,
                        help='Max distance from chr end to be considered subtelomeric (default: 100000 bp)')
    parser.add_argument('--min-mapq', type=int, default=30,
                        help='Minimum mapping quality to consider (default: 30)')
    args = parser.parse_args()

    print(f'Parsing PAF file: {args.paf}')
    best_hits = parse_paf(args.paf, min_mapq=args.min_mapq)
    print(f'  {len(best_hits)} contigs with high-quality alignments')

    print(f'Reading assembly FASTA: {args.assembly}')
    sequences = read_fasta(args.assembly)
    print(f'  {len(sequences)} contigs loaded')

    # --- Assign contigs to chromosome arms ---
    # arm_assignments: (chr_num, 'L'|'R') -> record dict
    # When multiple contigs compete for the same arm, keep the longest alignment.
    arm_assignments = {}

    for contig, record in best_hits.items():
        classification = classify_contig(record, threshold=args.threshold)
        chr_num = CHR_INFO[record['target_name']][0]

        if classification == 'internal':
            continue

        if classification == 'BOTH':
            split_pos = get_split_position(record)
            for arm in ('L', 'R'):
                key = (chr_num, arm)
                rec = dict(record, split_pos=split_pos, half=arm)
                if key not in arm_assignments or record['align_len'] > arm_assignments[key]['align_len']:
                    arm_assignments[key] = rec
        else:
            key = (chr_num, classification)
            if key not in arm_assignments or record['align_len'] > arm_assignments[key]['align_len']:
                arm_assignments[key] = record

    # --- Print summary ---
    all_chrs = set(range(1, 17))
    found_L = {c for (c, arm) in arm_assignments if arm == 'L'}
    found_R = {c for (c, arm) in arm_assignments if arm == 'R'}
    missing_L = sorted(all_chrs - found_L)
    missing_R = sorted(all_chrs - found_R)
    total_found = len(arm_assignments)

    print()
    print('=' * 65)
    print('CHROMOSOME END COVERAGE SUMMARY')
    print('=' * 65)
    print(f'L arms found : {len(found_L)}/16   chromosomes: {sorted(found_L)}')
    print(f'R arms found : {len(found_R)}/16   chromosomes: {sorted(found_R)}')
    print(f'Total        : {total_found}/32 chromosome ends found')
    if missing_L:
        print(f'Missing L    : chr{missing_L}')
    if missing_R:
        print(f'Missing R    : chr{missing_R}')

    print()
    print(f'{"Label":<10} {"Contig":<15} {"Chr region mapped":<30} {"AlignLen":>10}  Note')
    print('-' * 75)
    for chr_num in range(1, 17):
        for arm in ('L', 'R'):
            key = (chr_num, arm)
            label = f'chr{chr_num}{arm}'
            if key in arm_assignments:
                rec = arm_assignments[key]
                region = f'{rec["target_start"]}-{rec["target_end"]}'
                note = f'split at q:{rec["split_pos"]}' if 'split_pos' in rec else ''
                print(f'{label:<10} {rec["query_name"]:<15} {region:<30} {rec["align_len"]:>10}  {note}')
            else:
                print(f'{label:<10} {"NOT FOUND":<15}')

    # --- Write output FASTA ---
    print()
    print(f'Writing labeled sequences to: {args.output}')
    written = 0
    missing_seq = []

    with open(args.output, 'w') as out:
        for chr_num in range(1, 17):
            for arm in ('L', 'R'):
                key = (chr_num, arm)
                if key not in arm_assignments:
                    continue

                rec = arm_assignments[key]
                contig_name = rec['query_name']
                full_seq = sequences.get(contig_name)

                if full_seq is None:
                    missing_seq.append(contig_name)
                    continue

                label = f'chr{chr_num}{arm}'

                if 'half' in rec:
                    # Contig spans both arms — split at midpoint
                    split = rec['split_pos']
                    if rec['strand'] == '+':
                        seq = full_seq[:split] if arm == 'L' else full_seq[split:]
                    else:
                        # On - strand the contig runs reverse relative to the chromosome,
                        # so the L arm is at the 3' end of the contig and R arm at the 5' end
                        seq = full_seq[split:] if arm == 'L' else full_seq[:split]
                else:
                    seq = full_seq

                # For - strand alignments the chromosome end is at the wrong end of the
                # extracted sequence. RC so that:
                #   L-arm sequences always have the telomere at position 0 (left)
                #   R-arm sequences always have the telomere at the rightmost position
                if rec['strand'] == '-':
                    seq = reverse_complement(seq)

                write_fasta_seq(out, label, seq)
                written += 1

    if missing_seq:
        print(f'WARNING: sequences not found in assembly for: {missing_seq}')
    print(f'Done. Wrote {written} sequences.')


if __name__ == '__main__':
    main()
