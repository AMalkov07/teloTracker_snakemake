"""
Step 12: Summarize recombination results across all chromosome ends.

Reads all per-chr-end *_features.tsv files and produces a single summary TSV
with one row per chr end.

Usage:
  python aggregate_recombination.py --summarize \
      --recombination-dir  results/{base}/recombination/ \
      --base-name          {base} \
      --output-summary     results/{base}/recombination/{base}_recombination_summary.tsv
"""

import argparse
import glob
import os
import sys

import sys

# Ensure scripts/ is on the import path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from recombination_utils import write_results_tsv


def parse_args():
    p = argparse.ArgumentParser(description='Summarize recombination results')
    p.add_argument('--recombination-dir', required=True)
    p.add_argument('--base-name',         required=True)
    p.add_argument('--output-summary',    required=True)
    return p.parse_args()


def summarize(recombination_dir, base_name, output_summary):
    """Read all *_features.tsv files, produce per-chr-end summary."""
    pattern = os.path.join(recombination_dir, f'{base_name}_*_features.tsv')
    files = sorted(glob.glob(pattern))

    if not files:
        print(f'  No features TSV files found matching: {pattern}')
        write_results_tsv([], output_summary)
        return

    summary_rows = []
    skipped_rows = []
    for fpath in files:
        # Extract chr_end from filename
        basename = os.path.basename(fpath)
        # Pattern: {base_name}_{chr_end}_features.tsv
        suffix = basename.replace(f'{base_name}_', '').replace('_features.tsv', '')
        chr_end = suffix

        # Was this chr_end skipped by analyze_features.py due to low coverage?
        sidecar_path = fpath + '.skipped'
        if os.path.exists(sidecar_path):
            try:
                with open(sidecar_path) as fh:
                    reason = fh.read().strip()
            except Exception:
                reason = 'skipped (reason file unreadable)'
            row = {
                'chr_end': chr_end,
                'total_reads': 0,
                'status': 'skipped',
                'skip_reason': reason,
            }
            summary_rows.append(row)
            skipped_rows.append(row)
            continue

        try:
            df = pd.read_csv(fpath, sep='\t')
        except Exception as e:
            print(f'  Warning: could not read {fpath}: {e}')
            continue

        if df.empty:
            summary_rows.append({
                'chr_end': chr_end,
                'total_reads': 0,
                'status': 'no_events',
                'skip_reason': '',
            })
            continue

        total = len(df)
        n_recomb = df['recombination_detected'].sum() if 'recombination_detected' in df.columns else 0
        n_no_recomb = total - n_recomb

        # Y prime status counts
        yp_col = 'y_prime_recombination_status'
        yp_no_change = len(df[df[yp_col] == 'No Change']) if yp_col in df.columns else 0
        yp_change = total - yp_no_change

        # Spacer recombination counts
        sp_col = 'spacer_recombination'
        sp_switch = len(df[df[sp_col].isin(['switch_detected', 'full_switch'])]) if sp_col in df.columns else 0

        # X element recombination counts
        xe_col = 'x_element_recombination'
        xe_switch = len(df[df[xe_col].isin(['switch_detected', 'full_switch'])]) if xe_col in df.columns else 0

        # Confidence stats
        conf_col = 'overall_confidence'
        mean_conf = df[conf_col].mean() if conf_col in df.columns else 0.0

        # Cross-feature consistency
        cf_col = 'cross_feature_consistent'
        n_consistent = df[cf_col].sum() if cf_col in df.columns else 0

        # Complex events
        cx_col = 'is_complex_event'
        n_complex = df[cx_col].sum() if cx_col in df.columns else 0

        # Most common recombination source
        src_col = 'recombination_source'
        if src_col in df.columns:
            sources = df[df[src_col] != ''][src_col]
            most_common_source = sources.mode().iloc[0] if not sources.empty else ''
        else:
            most_common_source = ''

        summary_rows.append({
            'chr_end': chr_end,
            'total_reads': total,
            'status': 'analyzed',
            'skip_reason': '',
            'n_recombination': int(n_recomb),
            'n_no_recombination': int(n_no_recomb),
            'pct_recombination': round(100 * n_recomb / max(total, 1), 1),
            'n_spacer_switch': sp_switch,
            'n_x_element_switch': xe_switch,
            'n_y_prime_change': yp_change,
            'n_y_prime_no_change': yp_no_change,
            'mean_confidence': round(mean_conf, 4),
            'n_cross_feature_consistent': int(n_consistent),
            'n_complex_events': int(n_complex),
            'most_common_source': most_common_source,
        })

    write_results_tsv(summary_rows, output_summary)
    print(f'  Summary: {len(summary_rows)} chr ends written to {output_summary}')

    # Print overview
    total_reads = sum(r['total_reads'] for r in summary_rows)
    total_recomb = sum(r.get('n_recombination', 0) for r in summary_rows)
    print(f'  Total: {total_reads} reads, {total_recomb} recombination events '
          f'({100 * total_recomb / max(total_reads, 1):.1f}%)')

    # Call out any chr_ends that were skipped due to low coverage
    if skipped_rows:
        print(f'  SKIPPED: {len(skipped_rows)} of {len(summary_rows)} chr_ends had insufficient coverage:')
        for r in skipped_rows:
            print(f'    {r["chr_end"]:<8} {r["skip_reason"]}')


def main():
    args = parse_args()
    print(f'aggregate_recombination.py -- summarize mode')
    summarize(args.recombination_dir, args.base_name, args.output_summary)


if __name__ == '__main__':
    main()
