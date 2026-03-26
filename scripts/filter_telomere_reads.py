"""
Filter per-chr-end FASTAs to only include reads with confirmed telomeric repeats.

Reads the post_telo_trimming.tsv to identify reads with non-empty repeat_length
(confirmed telomere detected), then writes filtered FASTAs containing only those reads.

Usage:
  python filter_telomere_reads.py \
      --telo-tsv       results/{base}/{base}_post_telo_trimming.tsv \
      --input-dir      results/{base}/blast/chr_anchor_reads/ \
      --output-dir     results/{base}/telomere_filtered_reads/ \
      --base-name      {base} \
      --anchor-set     {anchor}
"""

import argparse
import os
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description='Filter reads to telomere-confirmed only')
    p.add_argument('--telo-tsv',    required=True, help='post_telo_trimming.tsv')
    p.add_argument('--input-dir',   required=True, help='Directory with per-chr-end FASTAs')
    p.add_argument('--output-dir',  required=True, help='Output directory for filtered FASTAs')
    p.add_argument('--base-name',   required=True)
    p.add_argument('--anchor-set',  required=True)
    return p.parse_args()


def load_telomere_read_ids(telo_tsv):
    """Return set of read_ids that have confirmed telomeric repeats."""
    df = pd.read_csv(telo_tsv, sep='\t')
    # Reads with non-empty repeat_length have confirmed telomere
    has_telo = df[df['repeat_length'].notna() & (df['repeat_length'] > 0)]
    return set(has_telo['read_id'].values)


def filter_fasta(input_path, output_path, keep_ids):
    """Write only reads whose ID is in keep_ids to output_path."""
    kept = 0
    total = 0
    with open(input_path) as fin, open(output_path, 'w') as fout:
        writing = False
        for line in fin:
            if line.startswith('>'):
                total += 1
                read_id = line[1:].strip().split()[0]
                writing = read_id in keep_ids
                if writing:
                    kept += 1
                    fout.write(line)
            elif writing:
                fout.write(line)
    return total, kept


def main():
    args = parse_args()
    print(f'filter_telomere_reads.py')

    if not os.path.exists(args.telo_tsv):
        print(f'ERROR: missing {args.telo_tsv}', file=sys.stderr)
        sys.exit(1)

    # Load telomere-confirmed read IDs
    telo_ids = load_telomere_read_ids(args.telo_tsv)
    print(f'  {len(telo_ids)} reads with confirmed telomere')

    os.makedirs(args.output_dir, exist_ok=True)

    # Filter each per-chr-end FASTA
    total_reads = 0
    total_kept = 0
    chr_ends_processed = 0

    for chr_num in range(1, 17):
        for side in ['L', 'R']:
            chr_end = f'chr{chr_num}{side}'
            input_name = f'{args.base_name}_blasted_{args.anchor_set}_{chr_end}_anchor_reads.fasta'
            input_path = os.path.join(args.input_dir, input_name)
            output_path = os.path.join(args.output_dir, f'{args.base_name}_{chr_end}_telomere_reads.fasta')

            if not os.path.exists(input_path):
                # Create empty output so Snakemake doesn't complain
                open(output_path, 'w').close()
                continue

            n_total, n_kept = filter_fasta(input_path, output_path, telo_ids)
            total_reads += n_total
            total_kept += n_kept
            chr_ends_processed += 1

            if n_total > 0:
                pct = 100 * n_kept / n_total
                print(f'  {chr_end}: {n_kept}/{n_total} reads kept ({pct:.1f}%)')

    print(f'  Summary: {total_kept}/{total_reads} reads kept across {chr_ends_processed} chr ends '
          f'({100 * total_kept / max(total_reads, 1):.1f}%)')


if __name__ == '__main__':
    main()
