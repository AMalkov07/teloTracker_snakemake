#!/usr/bin/env python3
"""
Extract per-event TSVs from per-chr_end features.tsv files.

For each `<base>_<chr_end>_features.tsv` in `_pipeline/recombination/`:
  - Filter rows where `recombination_detected == True`.
  - If 0 events, skip this chr_end.
  - Else write two TSVs into the events_dir:
      <base>_<chr_end>_events_full.tsv     — all 46 columns, only event rows.
      <base>_<chr_end>_events_summary.tsv  — slim 13-column view.

Plus one cross-chr_end aggregate:
  <base>_all_events_summary.tsv  — every chr_end's events_summary concatenated.

Plus a README.md describing the columns.

Usage:
    python scripts/extract_recombination_events.py <recomb_dir> <base_name> <events_dir>
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd


SUMMARY_COLUMNS = [
    'read_id',
    'chr_end',
    'read_length',
    'recombination_source',
    'overall_confidence',
    'switch_pos',
    'evidence',
    'spacer_recombination',
    'x_element_recombination',
    'y_prime_recombination_status',
    'y_prime_observed_array',
    'is_complex_event',
    'cross_feature_consistent',
]


README_TEXT = """\
# recombination_events/

Filtered, slim view of the per-read features.tsv data in
`_pipeline/recombination/`. Only reads where `recombination_detected==True`
appear here. chr_ends with zero such events have no files in this folder.

## Files

- `<base>_<chr_end>_events_full.tsv`    — complete 46-column features.tsv,
  restricted to event rows. Use when you need every column.
- `<base>_<chr_end>_events_summary.tsv` — slim 13-column view of the same
  rows. Easier to scan when you just want the recombination call.
- `<base>_all_events_summary.tsv`       — every chr_end's events_summary
  concatenated, sorted by chr_end then read_id.

## Slim-summary columns

| Column | Meaning |
|---|---|
| read_id | Read identifier (UUID, possibly compound for merged reads) |
| chr_end | Which chr_end's analysis this event came from (e.g. chr14L) |
| read_length | Total read length in bp |
| recombination_source | Most likely donor chr_end ('ambiguous' if no clear donor) |
| overall_confidence | Confidence score 0-1 (higher = more reliable call) |
| switch_pos | Earliest non-(-1) switch position from spacer or x-element axes; -1 if no positional switch |
| evidence | Human-readable axis-by-axis call summary (from cross_feature_detail) |
| spacer_recombination | Spacer axis call: no_change / switch_detected / full_switch / no_data |
| x_element_recombination | X-element axis call: same vocabulary |
| y_prime_recombination_status | Y' axis: No Change / 1st Y' Change / Y' Recombination / Y' Loss / Y' Gain |
| y_prime_observed_array | Y' cluster IDs detected in this read, ordered (e.g. ID1,ID2,ID2) |
| is_complex_event | True if multiple non-self sources are involved |
| cross_feature_consistent | True if all axes point to the same donor |
"""


def derive_switch_pos(row):
    """Earliest non-(-1) switch position across the positional axes."""
    candidates = []
    for col in ('spacer_switch_pos', 'x_element_switch_pos'):
        v = row.get(col, -1)
        try:
            v = int(v)
        except (ValueError, TypeError):
            continue
        if v >= 0:
            candidates.append(v)
    return min(candidates) if candidates else -1


def build_slim(events: pd.DataFrame) -> pd.DataFrame:
    """Build the slim 13-column summary from the full event rows."""
    out = pd.DataFrame()
    out['read_id'] = events.get('read_id', '')
    out['chr_end'] = events.get('chr_end', '')
    out['read_length'] = events.get('read_length', 0)
    out['recombination_source'] = events.get('recombination_source', '')
    out['overall_confidence'] = events.get('overall_confidence', 0.0)
    out['switch_pos'] = events.apply(derive_switch_pos, axis=1)
    out['evidence'] = events.get('cross_feature_detail', '')
    out['spacer_recombination'] = events.get('spacer_recombination', '')
    out['x_element_recombination'] = events.get('x_element_recombination', '')
    out['y_prime_recombination_status'] = events.get('y_prime_recombination_status', '')
    out['y_prime_observed_array'] = events.get('y_prime_observed_array', '')
    out['is_complex_event'] = events.get('is_complex_event', '')
    out['cross_feature_consistent'] = events.get('cross_feature_consistent', '')
    # Reorder explicitly so columns are stable across runs even if pandas reorders.
    return out[SUMMARY_COLUMNS]


def parse_args():
    p = argparse.ArgumentParser(description='Extract recombination events from per-chr_end features.tsv files')
    p.add_argument('recombination_dir', help='Directory containing <base>_<chr_end>_features.tsv files')
    p.add_argument('base_name', help='Sample base name')
    p.add_argument('events_dir', help='Output directory for the events TSVs and README')
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.events_dir, exist_ok=True)

    pattern = os.path.join(args.recombination_dir, f'{args.base_name}_*_features.tsv')
    files = sorted(glob.glob(pattern))
    if not files:
        print(f'No features.tsv files found at: {pattern}')
        # Still write README so users see something explanatory.
        (open(os.path.join(args.events_dir, 'README.md'), 'w')).write(README_TEXT)
        return

    all_summary_rows = []
    n_chr_with_events = 0
    n_total_events = 0

    for fpath in files:
        basename = os.path.basename(fpath)
        suffix = basename.replace(f'{args.base_name}_', '').replace('_features.tsv', '')
        chr_end = suffix

        try:
            df = pd.read_csv(fpath, sep='\t')
        except Exception as e:
            print(f'  WARNING: could not read {fpath}: {e}')
            continue
        if df.empty or 'recombination_detected' not in df.columns:
            continue

        # The TSV stores booleans as the strings "True"/"False" — pd.read_csv may give either bool or str.
        flag = df['recombination_detected']
        if flag.dtype == bool:
            mask = flag
        else:
            mask = flag.astype(str).str.lower() == 'true'

        events = df[mask].copy()
        if events.empty:
            continue

        n_chr_with_events += 1
        n_total_events += len(events)

        full_path = os.path.join(args.events_dir, f'{args.base_name}_{chr_end}_events_full.tsv')
        events.to_csv(full_path, sep='\t', index=False)

        slim = build_slim(events)
        slim_path = os.path.join(args.events_dir, f'{args.base_name}_{chr_end}_events_summary.tsv')
        slim.to_csv(slim_path, sep='\t', index=False)

        all_summary_rows.append(slim)
        print(f'  {chr_end}: {len(events)} events  ->  {os.path.basename(slim_path)}')

    if all_summary_rows:
        all_df = pd.concat(all_summary_rows, ignore_index=True)
        all_df = all_df.sort_values(['chr_end', 'read_id']).reset_index(drop=True)
        all_path = os.path.join(args.events_dir, f'{args.base_name}_all_events_summary.tsv')
        all_df.to_csv(all_path, sep='\t', index=False)
        print(f'\n  Aggregate: {len(all_df)} events across {n_chr_with_events} chr_ends -> {os.path.basename(all_path)}')
    else:
        print('\n  No recombination events found in any chr_end.')

    readme_path = os.path.join(args.events_dir, 'README.md')
    with open(readme_path, 'w') as fh:
        fh.write(README_TEXT)


if __name__ == '__main__':
    main()
