#!/usr/bin/env python3
"""
Per-chr_end track plots for recombination events.

For each `<base>_<chr_end>_features.tsv` in `_pipeline/recombination/`:
  - Filter rows where `recombination_detected == True`.
  - If 0 events, skip (no plot for this chr_end — explicit per-design).
  - Else build a stacked-track PNG: one horizontal row per event-read,
    showing colored segments along the read coordinate (anchor → spacer
    → x-element → Y' → telomere).

Usage:
    python scripts/plot_recombination_tracks.py <recomb_dir> <base_name> <output_dir>

Output:
    <output_dir>/<base>_<chr_end>_tracks.png   (one per chr_end with events)
    <output_dir>/.plots_done                   (marker)
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# Categorical color for each Y′ status. Choose visually distinct.
Y_STATUS_COLORS = {
    'No Change':         '#A8D8A8',  # pale green (rare here — these are events)
    "1st Y' Change":     '#9B59B6',  # purple
    "Y' Recombination":  '#E67E22',  # orange
    "Y' Loss":           '#E74C3C',  # red
    "Y' Gain":           '#3498DB',  # blue
    'no_data':           '#BDBDBD',  # grey
    '':                  '#BDBDBD',
}

# Deterministic chr_end → color mapping. Reused across all plots so the same
# source name always shows the same color.
CHR_END_LIST = [f'chr{n}{s}' for n in range(1, 17) for s in ('L', 'R')]

def make_chr_end_palette():
    palette = sns.color_palette('husl', n_colors=len(CHR_END_LIST))
    return {ce: palette[i] for i, ce in enumerate(CHR_END_LIST)}


SOURCE_COLORS = make_chr_end_palette()
SPECIAL_SOURCE_COLORS = {
    '':          '#BDBDBD',
    'ambiguous': '#7F7F7F',
    'self':      '#CCCCCC',
}

ANCHOR_COLOR    = '#404040'
TELOMERE_COLOR  = '#000000'
BACKGROUND_BAR  = '#F0F0F0'


def source_color(name: str, expected_chr_end: str) -> str:
    """Color for a 'source' label."""
    if name in SPECIAL_SOURCE_COLORS:
        return SPECIAL_SOURCE_COLORS[name]
    if name == expected_chr_end:
        return SPECIAL_SOURCE_COLORS['self']
    return SOURCE_COLORS.get(name, '#888888')


def safe_int(v, default=-1):
    try:
        return int(v)
    except (ValueError, TypeError):
        return default


def is_truthy(v) -> bool:
    if isinstance(v, bool):
        return v
    return str(v).lower() == 'true'


def display_segment_ends(row, rl):
    """
    Compute the four cumulative display-end positions (anchor, spacer, x-element,
    Y') for a read so that the anchor always sits on the LEFT of the bar
    regardless of the read's actual sequencing orientation.

    For telo_side='end' reads (anchor at 5', telomere at 3'), the *_end columns
    are already in left-to-right display order; we use them directly.

    For telo_side='beginning' reads (telomere at 5', anchor at 3'), the read
    was sequenced in the reverse-complement orientation. We mirror it for
    display by computing each segment's right edge as `rl - <feature>_start`.
    Returns -1 (= no segment) for any feature where the start/end column is
    -1 in the source data.
    """
    telo_side = str(row.get('telo_side', 'end')).strip().lower()

    if telo_side == 'beginning':
        def mirrored_end(start_col):
            v = safe_int(row.get(start_col, -1))
            return rl - v if v >= 0 else -1
        a_end = mirrored_end('anchor_start')
        s_end = mirrored_end('spacer_start')
        x_end = mirrored_end('x_element_start')
        y_end = mirrored_end('y_prime_start')
    else:
        a_end = safe_int(row.get('anchor_end', -1))
        s_end = safe_int(row.get('spacer_end', -1))
        x_end = safe_int(row.get('x_element_end', -1))
        y_end = safe_int(row.get('y_prime_end', -1))

    # Clamp -1's to 0 and enforce monotonic ordering so segments fall in
    # left-to-right order without overlap.
    a_end = max(a_end, 0)
    s_end = max(s_end, a_end)
    x_end = max(x_end, s_end)
    y_end = max(y_end, x_end)
    return a_end, s_end, x_end, y_end


def display_switch_pos(row, original_pos, rl):
    """Mirror a switch position for telo_side='beginning' reads."""
    if original_pos < 0:
        return -1
    telo_side = str(row.get('telo_side', 'end')).strip().lower()
    if telo_side == 'beginning':
        return max(rl - original_pos, 0)
    return original_pos


def plot_chr_end_tracks(events: pd.DataFrame, chr_end: str, base_name: str,
                       output_path: str) -> None:
    """Render the track plot for one chr_end."""
    # Sort: group by status, then by recombination_source, then by switch position.
    sort_cols = []
    for c in ('y_prime_recombination_status', 'recombination_source', 'spacer_switch_pos'):
        if c in events.columns:
            sort_cols.append(c)
    if sort_cols:
        events = events.sort_values(sort_cols, kind='stable').reset_index(drop=True)

    n = len(events)
    if n == 0:
        return

    max_len = int(events['read_length'].max()) if 'read_length' in events.columns else 1

    # Layout: 0.18 inches per row, capped at 18 inches; 12 inches wide.
    height = max(2.5, min(18.0, 0.18 * n + 1.5))
    fig, ax = plt.subplots(figsize=(12.0, height))

    bar_height = 0.7
    yticks = []
    yticklabels = []

    for i, row in events.iterrows():
        y = n - 1 - i  # top-down ordering
        rl = safe_int(row.get('read_length', 0), 0)
        a_end, s_end, x_end, y_end = display_segment_ends(row, rl)

        # Background full-read bar so even short reads are visible.
        ax.barh(y, rl, left=0, height=bar_height, color=BACKGROUND_BAR, edgecolor='none')

        # Anchor: 0 → a_end
        if a_end > 0:
            ax.barh(y, a_end, left=0, height=bar_height,
                    color=ANCHOR_COLOR, edgecolor='none')

        # Spacer: a_end → s_end
        if s_end > a_end:
            ax.barh(y, s_end - a_end, left=a_end, height=bar_height,
                    color=source_color(str(row.get('spacer_source', '')), chr_end),
                    edgecolor='none')

        # X-element: s_end → x_end
        if x_end > s_end:
            ax.barh(y, x_end - s_end, left=s_end, height=bar_height,
                    color=source_color(str(row.get('x_element_source', '')), chr_end),
                    edgecolor='none')

        # Y' region: x_end → y_end (color by status, since Y' analysis is categorical)
        status = str(row.get('y_prime_recombination_status', ''))
        if y_end > x_end:
            ax.barh(y, y_end - x_end, left=x_end, height=bar_height,
                    color=Y_STATUS_COLORS.get(status, '#BDBDBD'),
                    edgecolor='none')

        # Telomere region: y_end → rl
        if rl > y_end:
            ax.barh(y, rl - y_end, left=y_end, height=bar_height,
                    color=TELOMERE_COLOR, edgecolor='none')

        # Switch markers (vertical ticks at known switch positions).
        for col in ('spacer_switch_pos', 'x_element_switch_pos'):
            sp = display_switch_pos(row, safe_int(row.get(col, -1)), rl)
            if sp > 0:
                ax.vlines(sp, y - 0.45, y + 0.45, colors='white', linewidth=0.7)

        # y-axis label: short read_id + status code
        rid = str(row.get('read_id', ''))
        rid_short = rid[:8] + ('…' if len(rid) > 8 else '')
        yticks.append(y)
        yticklabels.append(f'{rid_short}  [{status[:14]}]' if status else rid_short)

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=6)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_xlim(0, max(max_len, 1) * 1.02)
    ax.set_xlabel('Position on read (bp)')
    ax.set_title(f'{base_name}: {chr_end} — {n} recombination event(s)', fontsize=11)
    ax.grid(axis='x', linestyle=':', alpha=0.3)

    # Build a compact legend showing the colors that actually appear in this plot.
    legend_handles = []
    seen_sources = set()
    for col in ('spacer_source', 'x_element_source'):
        if col in events.columns:
            for v in events[col].dropna().astype(str).unique():
                if v and v not in seen_sources:
                    seen_sources.add(v)
                    legend_handles.append(
                        mpatches.Patch(color=source_color(v, chr_end), label=f'src: {v}'))
    if 'y_prime_recombination_status' in events.columns:
        for v in events['y_prime_recombination_status'].dropna().astype(str).unique():
            if v and v not in Y_STATUS_COLORS.keys() | {'No Change', ''}:
                continue
            if v in ('', 'No Change'):
                continue
            legend_handles.append(mpatches.Patch(color=Y_STATUS_COLORS[v], label=f"y′: {v}"))
    legend_handles.append(mpatches.Patch(color=ANCHOR_COLOR, label='anchor'))
    legend_handles.append(mpatches.Patch(color=TELOMERE_COLOR, label='telomere'))
    ax.legend(handles=legend_handles, loc='center left',
              bbox_to_anchor=(1.01, 0.5), fontsize=7, frameon=False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(description='Per-chr_end recombination track plots')
    p.add_argument('recombination_dir')
    p.add_argument('base_name')
    p.add_argument('output_dir')
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pattern = os.path.join(args.recombination_dir, f'{args.base_name}_*_features.tsv')
    files = sorted(glob.glob(pattern))
    if not files:
        print(f'No features.tsv files found at: {pattern}')
    else:
        n_plots = 0
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
            flag = df['recombination_detected']
            if flag.dtype == bool:
                mask = flag
            else:
                mask = flag.astype(str).str.lower() == 'true'
            events = df[mask].copy()
            if events.empty:
                continue
            out_png = os.path.join(args.output_dir, f'{args.base_name}_{chr_end}_tracks.png')
            plot_chr_end_tracks(events, chr_end, args.base_name, out_png)
            print(f'  {chr_end}: {len(events)} events  ->  {os.path.basename(out_png)}')
            n_plots += 1
        print(f'\n  Wrote {n_plots} track plot(s) to {args.output_dir}')

    # Marker so Snakemake can track completion.
    marker = os.path.join(args.output_dir, '.plots_done')
    with open(marker, 'w') as fh:
        fh.write('done\n')


if __name__ == '__main__':
    main()
