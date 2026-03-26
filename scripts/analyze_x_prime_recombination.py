"""
Step 11b: X prime (X element) recombination analysis.

For reads where the minimap2 breakpoint falls in an x_variable or x_core feature,
BLAST the diverged tail against the X element library to identify the source chr end.
X elements are highly similar across chr ends so results will often be ambiguous —
this is expected and reflected in low confidence scores.

Usage:
  python analyze_x_prime_recombination.py \
      --breakpoints-tsv  results/{base}/recombination/{base}_{chr_end}_breakpoints.tsv \
      --diverged-fasta   results/{base}/recombination/{base}_{chr_end}_diverged_tails.fasta \
      --x-element-lib-dir  references/pairings_for_x_element_ends/{strain}_pairings/ \
      --chr-end  chr4R  --strain 6991 \
      --output-tsv  results/{base}/recombination/{base}_{chr_end}_x_prime_recomb.tsv \
      --threads 4
"""

import argparse
import glob as glob_module
import os
import sys
import tempfile

import pandas as pd

from recombination_utils import (
    Hypothesis,
    hypotheses_to_row_dict,
    get_first_breakpoint_feature_type,
    get_first_breakpoint_is_mid_element,
    normalize_hypotheses,
    read_fasta,
    run_blast,
    write_fasta_list,
    write_results_tsv,
    MIN_CONFIDENCE_TO_REPORT,
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description='X prime recombination analysis')
    p.add_argument('--breakpoints-tsv',    required=True)
    p.add_argument('--diverged-fasta',     required=True)
    p.add_argument('--x-element-lib-dir',  required=True)
    p.add_argument('--chr-end',            required=True)
    p.add_argument('--strain',             required=True)
    p.add_argument('--output-tsv',         required=True)
    p.add_argument('--threads',            type=int, default=4)
    return p.parse_args()

# ---------------------------------------------------------------------------
# Build combined X element library
# ---------------------------------------------------------------------------

def build_combined_x_element_fasta(x_element_lib_dir, strain, tmp_dir):
    """
    Concatenate all {strain}_paired_*.fasta files into one combined FASTA.
    Returns path to the combined file.
    """
    pattern = os.path.join(x_element_lib_dir, f'{strain}_paired_*.fasta')
    fasta_files = sorted(glob_module.glob(pattern))
    if not fasta_files:
        return None

    combined_path = os.path.join(tmp_dir, f'{strain}_x_element_combined.fasta')
    with open(combined_path, 'w') as out_fh:
        for fasta_path in fasta_files:
            with open(fasta_path) as in_fh:
                for line in in_fh:
                    out_fh.write(line)
    return combined_path

# ---------------------------------------------------------------------------
# Source chr end extraction from BLAST subject name
# ---------------------------------------------------------------------------

def subject_to_chr_end(sseqid):
    """
    Extract the chr end from an X element pairing sequence ID.
    Pairing FASTA sequence names typically contain the chr end, e.g.:
      chr4R_x_ends_section1
    Returns e.g. 'chr4R' or the full sseqid if no match.
    """
    import re
    m = re.search(r'(chr\d+[LR])', sseqid)
    return m.group(1) if m else sseqid

# ---------------------------------------------------------------------------
# Hypothesis generation
# ---------------------------------------------------------------------------

def generate_x_prime_hypotheses(read_id, bp_row, blast_hits, diverged_seq_len):
    """
    Generate Hypothesis objects from X element BLAST hits for one read.
    """
    hypotheses = []
    is_mid = get_first_breakpoint_is_mid_element(bp_row)

    if blast_hits.empty:
        hypotheses.append(Hypothesis(
            h_type='ambiguous',
            description='X element breakpoint but no match in library',
            confidence=0.2,
            ambiguous=True,
        ))
        return normalize_hypotheses(hypotheses)

    # Group by source chr end; take max bitscore per chr end
    blast_hits = blast_hits.copy()
    blast_hits['source_chr_end'] = blast_hits['sseqid'].apply(subject_to_chr_end)
    best_by_source = blast_hits.groupby('source_chr_end')['bitscore'].max().reset_index()
    best_by_source = best_by_source.sort_values('bitscore', ascending=False)

    top_bitscore = best_by_source.iloc[0]['bitscore']
    second_bitscore = best_by_source.iloc[1]['bitscore'] if len(best_by_source) > 1 else 0

    # Base confidence from best hit
    best_hit = blast_hits[blast_hits['bitscore'] == top_bitscore].iloc[0]
    c_base = (best_hit['pident'] / 100.0) * min(
        1.0, best_hit['length'] / max(diverged_seq_len, 1)
    ) ** 0.5

    if is_mid:
        c_base *= 0.7

    # Determine if ambiguous (all within 10% of top bitscore)
    threshold = top_bitscore * 0.9
    ambiguous_sources = best_by_source[best_by_source['bitscore'] >= threshold]

    if len(ambiguous_sources) > 1 or top_bitscore < second_bitscore * 2:
        # Ambiguous: distribute confidence equally
        n_sources = len(ambiguous_sources)
        c_each = c_base / n_sources
        for _, row in ambiguous_sources.iterrows():
            if c_each < MIN_CONFIDENCE_TO_REPORT:
                break
            hypotheses.append(Hypothesis(
                h_type='x_prime',
                description=f"X element switch to {row['source_chr_end']} (ambiguous — {n_sources} equally-matched sources)",
                confidence=c_each,
                source_chr_ends=[row['source_chr_end']],
                ambiguous=True,
            ))
    else:
        # Clear winner
        top_source = best_by_source.iloc[0]['source_chr_end']
        hypotheses.append(Hypothesis(
            h_type='x_prime',
            description=f"X element switch to {top_source}",
            confidence=c_base,
            source_chr_ends=[top_source],
            ambiguous=False,
        ))

    return normalize_hypotheses(hypotheses)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    print(f'analyze_x_prime_recombination.py — chr_end={args.chr_end}')

    if not os.path.exists(args.breakpoints_tsv):
        print(f'ERROR: missing breakpoints TSV: {args.breakpoints_tsv}', file=sys.stderr)
        sys.exit(1)

    df_bp = pd.read_csv(args.breakpoints_tsv, sep='\t')
    if df_bp.empty:
        write_results_tsv([], args.output_tsv)
        return

    # Only process reads with x-element breakpoints
    df_bp['_first_bp_feat_type'] = df_bp['breakpoint_feature_types'].astype(str).apply(lambda v: v.split(';')[0])
    x_reads = df_bp[df_bp['_first_bp_feat_type'].isin(['x_variable', 'x_core'])]
    print(f'  {len(x_reads)} reads with X element breakpoints (of {len(df_bp)} total)')

    # Load diverged tails
    diverged_seqs = {}
    if os.path.exists(args.diverged_fasta):
        diverged_seqs = read_fasta(args.diverged_fasta)

    rows = []

    if x_reads.empty or not os.path.exists(args.x_element_lib_dir):
        # No X element reads or no library — write a pass-through row for all reads
        for _, bp_row in df_bp.iterrows():
            row = {'read_id': bp_row['read_id'], 'chr_end': args.chr_end,
                   'x_prime_analyzed': False, 'n_hypotheses': 0}
            row.update(hypotheses_to_row_dict([]))
            rows.append(row)
        write_results_tsv(rows, args.output_tsv)
        return

    with tempfile.TemporaryDirectory() as tmp_dir:
        combined_lib = build_combined_x_element_fasta(args.x_element_lib_dir, args.strain, tmp_dir)
        if combined_lib is None:
            print(f'  WARNING: no X element FASTA files found in {args.x_element_lib_dir}')
            for _, bp_row in df_bp.iterrows():
                row = {'read_id': bp_row['read_id'], 'chr_end': args.chr_end,
                       'x_prime_analyzed': False, 'n_hypotheses': 0}
                row.update(hypotheses_to_row_dict([]))
                rows.append(row)
            write_results_tsv(rows, args.output_tsv)
            return

        for _, bp_row in df_bp.iterrows():
            read_id = bp_row['read_id']
            feat_type = get_first_breakpoint_feature_type(bp_row)

            if feat_type not in ('x_variable', 'x_core'):
                row = {'read_id': read_id, 'chr_end': args.chr_end,
                       'x_prime_analyzed': False, 'n_hypotheses': 0}
                row.update(hypotheses_to_row_dict([]))
                rows.append(row)
                continue

            # Find diverged sequence for this read
            diverged_seq = ''
            for header in diverged_seqs:
                if header.startswith(read_id) or read_id in header:
                    diverged_seq = diverged_seqs[header]
                    break

            if not diverged_seq:
                hypotheses = [Hypothesis(
                    h_type='ambiguous',
                    description='X element breakpoint but no diverged tail available',
                    confidence=0.2,
                    ambiguous=True,
                )]
            else:
                # Write diverged tail to temp FASTA and BLAST
                tail_fasta = os.path.join(tmp_dir, f'{read_id}_tail.fasta')
                write_fasta_list([(read_id, diverged_seq)], tail_fasta)

                blast_hits = run_blast(tail_fasta, combined_lib, tmp_dir, label=f'{read_id}_x')
                hypotheses = generate_x_prime_hypotheses(
                    read_id, bp_row, blast_hits, len(diverged_seq)
                )

            row = {
                'read_id': read_id,
                'chr_end': args.chr_end,
                'x_prime_analyzed': True,
                'n_hypotheses': len(hypotheses),
            }
            row.update(hypotheses_to_row_dict(hypotheses))
            rows.append(row)

    os.makedirs(os.path.dirname(args.output_tsv) or '.', exist_ok=True)
    write_results_tsv(rows, args.output_tsv)
    print(f'  Output: {args.output_tsv} ({len(rows)} reads)')


if __name__ == '__main__':
    main()
