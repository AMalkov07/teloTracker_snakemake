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

        # Most common recombination source.
        #   most_common_source        : mode over NAMED donor ends (excludes '' and 'ambiguous')
        #   most_common_source_legacy : the pre-v2 value (mode including the literal 'ambiguous')
        src_col = 'recombination_source'
        most_common_source = most_common_source_legacy = ''
        n_ambiguous = 0
        if src_col in df.columns:
            sources = df[src_col].dropna().astype(str)
            sources = sources[sources != '']
            _mode = sources.mode()
            most_common_source_legacy = _mode.iloc[0] if not _mode.empty else ''
            n_ambiguous = int((sources == 'ambiguous').sum())
            named = sources[sources != 'ambiguous']
            _mode = named.mode()
            most_common_source = _mode.iloc[0] if not _mode.empty else ''

        # Y' status breakdown (n_y_prime_change kept for compatibility = total - No Change)
        def _n(status):
            return int((df[yp_col] == status).sum()) if yp_col in df.columns else 0
        n_loss = _n("Y' Loss")
        n_loss_confirmed = 0
        if 'telomere_end_confirmed' in df.columns and yp_col in df.columns:
            n_loss_confirmed = int(((df[yp_col] == "Y' Loss") & (df['telomere_end_confirmed'].astype(str) == 'True')).sum())
        mech_col = 'recombination_mechanism'
        most_common_mechanism = ''
        if mech_col in df.columns:
            mechs = df[mech_col].dropna().astype(str)
            mechs = mechs[mechs != '']
            _mode = mechs.mode()
            most_common_mechanism = _mode.iloc[0] if not _mode.empty else ''

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
            'n_y_prime_gain': _n("Y' Gain"),
            'n_y_prime_loss': n_loss,
            'n_y_prime_loss_confirmed_end': n_loss_confirmed,
            'n_first_y_prime_change': _n("1st Y' Change"),
            'n_y_prime_recombination': _n("Y' Recombination"),
            'n_ambiguous': n_ambiguous,
            'mean_confidence': round(mean_conf, 4),
            'n_cross_feature_consistent': int(n_consistent),
            'n_complex_events': int(n_complex),
            'most_common_source': most_common_source,
            'most_common_source_legacy': most_common_source_legacy,
            'most_common_mechanism': most_common_mechanism,
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
