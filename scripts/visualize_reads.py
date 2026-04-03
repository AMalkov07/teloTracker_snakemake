#!/usr/bin/env python3
"""
Visualize reads from the recombination pipeline, showing features as colored boxes
on a read-length line graph.

Feature positions (spacer, x_element, y_prime) are read directly from the
features TSV output by analyze_features.py. No BAM or BED files needed.

Usage:
    # Single read by ID
    python scripts/visualize_reads.py \
        --recomb-dir recombination/ \
        --chr-end chr14R \
        --read-id 721f2dad-6d96-4bb8-85cd-82b227e96b36 \
        -o output.png

    # 10 random reads from a chr_end
    python scripts/visualize_reads.py \
        --recomb-dir recombination/ \
        --chr-end chr14R \
        --n-reads 10 \
        -o output.png
"""

import argparse
import os
import random
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import pandas as pd


# Feature colors
FEATURE_COLORS = {
    'anchor': '#4A4A4A',
    'spacer': '#D3D3D3',
    'x_element': '#2196F3',
    'y_prime': '#E53935',
    'telomere': '#81C784',
}

FEATURE_DISPLAY_NAMES = {
    'anchor': 'Anchor',
    'spacer': 'Spacer',
    'x_element': 'X element',
    'y_prime': "Y'",
    'telomere': 'Telomere',
}

RECOMB_BORDER_COLOR = '#FFD600'
RECOMB_HATCH = '///'


def parse_args():
    p = argparse.ArgumentParser(description='Visualize recombination pipeline reads')
    p.add_argument('--recomb-dir', required=True, help='Recombination output directory')
    p.add_argument('--chr-end', required=True, help='Chromosome end (e.g., chr14R)')
    p.add_argument('--read-id', default=None, help='Specific read ID to visualize')
    p.add_argument('--n-reads', type=int, default=10, help='Number of random reads (if no --read-id)')
    p.add_argument('-o', '--output', default='read_visualization.png', help='Output image path')
    p.add_argument('--seed', type=int, default=42, help='Random seed for read selection')
    return p.parse_args()


def load_features_df(recomb_dir, chr_end):
    """Load features TSV for a chr_end."""
    features_files = [f for f in os.listdir(recomb_dir)
                      if f.endswith(f'_{chr_end}_features.tsv')]
    if not features_files:
        print(f"ERROR: No features file found for {chr_end} in {recomb_dir}")
        sys.exit(1)

    features_path = os.path.join(recomb_dir, features_files[0])
    return pd.read_csv(features_path, sep='\t')


def build_feature_regions(row):
    """Build list of (start, end, feature_type, label, is_recomb, recomb_source)
    from a features TSV row."""

    regions = []
    read_length = row['read_length']
    telo_side = row['telo_side']
    chr_end = row['chr_end']

    # Anchor
    anchor_start = row.get('anchor_start', -1)
    anchor_end = row.get('anchor_end', -1)
    if anchor_start >= 0 and anchor_end > anchor_start:
        regions.append((anchor_start, anchor_end, 'anchor',
                        f"Anchor ({anchor_end - anchor_start}bp)",
                        False, None))

    # Spacer
    spacer_start = row.get('spacer_start', -1)
    spacer_end = row.get('spacer_end', -1)
    if spacer_start >= 0 and spacer_end > spacer_start:
        recomb = row.get('spacer_recombination', 'no_change')
        source = row.get('spacer_source', chr_end)
        is_recomb = recomb in ('full_switch', 'switch_detected') and source != chr_end
        regions.append((spacer_start, spacer_end, 'spacer',
                        f"Spacer ({spacer_end - spacer_start}bp)",
                        is_recomb, source if is_recomb else None))

    # X element
    x_start = row.get('x_element_start', -1)
    x_end = row.get('x_element_end', -1)
    if x_start >= 0 and x_end > x_start:
        recomb = row.get('x_element_recombination', 'no_change')
        source = row.get('x_element_source', chr_end)
        is_recomb = recomb in ('full_switch', 'switch_detected') and source != chr_end
        regions.append((x_start, x_end, 'x_element',
                        f"X elem ({x_end - x_start}bp)",
                        is_recomb, source if is_recomb else None))

    # Y primes — parse individual positions from y_prime_positions field
    yp_positions = row.get('y_prime_positions', '')
    yp_status = row.get('y_prime_recombination_status', 'No Change')
    yp_compatible = row.get('y_prime_compatible_ends', '')
    yp_is_recomb = yp_status not in ('No Change', '')

    if isinstance(yp_positions, str) and yp_positions:
        for i, entry in enumerate(yp_positions.split(';')):
            # Format: ID7:9209-16049
            m = re.match(r'(\w+):(\d+)-(\d+)', entry)
            if m:
                yp_id = m.group(1)
                yp_start = int(m.group(2))
                yp_end = int(m.group(3))
                label = f"Y'{i+1} {yp_id} ({yp_end - yp_start}bp)"
                regions.append((yp_start, yp_end, 'y_prime', label,
                                yp_is_recomb,
                                yp_compatible if yp_is_recomb else None))
    else:
        # Fall back to overall y_prime_start/end if no per-Y-prime positions
        yp_start = row.get('y_prime_start', -1)
        yp_end = row.get('y_prime_end', -1)
        if yp_start >= 0 and yp_end > yp_start:
            yp_obs = row.get('y_prime_observed_array', '')
            label = f"Y' {yp_obs} ({yp_end - yp_start}bp)"
            regions.append((yp_start, yp_end, 'y_prime', label,
                            yp_is_recomb,
                            yp_compatible if yp_is_recomb else None))

    # Telomere — infer from telo_side and the feature boundaries
    if regions:
        all_starts = [s for s, e, *_ in regions]
        all_ends = [e for s, e, *_ in regions]
        if telo_side == 'end':
            telo_start = max(all_ends)
            if telo_start < read_length - 50:
                regions.append((telo_start, read_length, 'telomere',
                                f"Telomere ({read_length - telo_start}bp)",
                                False, None))
        else:
            telo_end = min(all_starts)
            if telo_end > 50:
                regions.append((0, telo_end, 'telomere',
                                f"Telomere ({telo_end}bp)",
                                False, None))

    return regions


def draw_read(ax, regions, read_length, telo_side, read_id, features_info,
              y_pos=0, height=0.6):
    """Draw a single read with feature boxes."""
    chr_end = features_info.get('chr_end', '')

    # Read backbone line
    ax.plot([0, read_length], [y_pos, y_pos], color='#888888', linewidth=1, zorder=1)

    for rstart, rend, feat_type, label, is_recomb, recomb_source in sorted(regions):
        color = FEATURE_COLORS.get(feat_type, '#9E9E9E')
        width = rend - rstart

        edgecolor = RECOMB_BORDER_COLOR if is_recomb else 'black'
        linewidth = 2.5 if is_recomb else 0.5

        box = FancyBboxPatch(
            (rstart, y_pos - height / 2), width, height,
            boxstyle="round,pad=0",
            facecolor=color, edgecolor=edgecolor, linewidth=linewidth,
            alpha=0.85, zorder=2
        )
        ax.add_patch(box)

        if is_recomb:
            hatch_box = FancyBboxPatch(
                (rstart, y_pos - height / 2), width, height,
                boxstyle="round,pad=0",
                facecolor='none', edgecolor=edgecolor, linewidth=0.5,
                hatch=RECOMB_HATCH, alpha=0.4, zorder=2.5
            )
            ax.add_patch(hatch_box)

        # Label inside box
        if width > read_length * 0.03:
            # Use short label for display
            short_label = label.split('(')[0].strip()
            fontsize = 7 if width > read_length * 0.06 else 5.5
            text_color = 'white' if feat_type not in ('spacer', 'telomere') else 'black'
            ax.text(rstart + width / 2, y_pos, short_label,
                    ha='center', va='center', fontsize=fontsize,
                    color=text_color, fontweight='bold', zorder=3)

        # Recombination source annotation
        if is_recomb and width > read_length * 0.02 and recomb_source:
            ax.text(rstart + width / 2, y_pos - height / 2 - 0.08,
                    f'{recomb_source}',
                    ha='center', va='top', fontsize=5, color='#B8860B',
                    fontweight='bold', fontstyle='italic', zorder=3)

    # Info text above read
    short_id = read_id[:16] + '...'
    spacer_src = features_info.get('spacer_source', '?')
    spacer_recomb = features_info.get('spacer_recombination', '?')
    x_src = features_info.get('x_element_source', '?')
    x_recomb = features_info.get('x_element_recombination', '?')
    yp_obs = features_info.get('y_prime_observed_array', '?')
    yp_status = features_info.get('y_prime_recombination_status', '?')
    recomb = features_info.get('recombination_detected', '?')
    conf = features_info.get('overall_confidence', '?')

    info_text = (f"{short_id}  ({read_length:,} bp, telo={telo_side})\n"
                 f"spacer={spacer_src} [{spacer_recomb}]  "
                 f"x_elem={x_src} [{x_recomb}]  "
                 f"y'={yp_obs} [{yp_status}]  "
                 f"recomb={recomb} conf={conf}")
    ax.text(0, y_pos + height / 2 + 0.15, info_text,
            fontsize=5.5, va='bottom', ha='left', color='#333333', family='monospace')


def main():
    args = parse_args()

    features_df = load_features_df(args.recomb_dir, args.chr_end)

    # Select reads
    if args.read_id:
        read_ids = [args.read_id]
    else:
        random.seed(args.seed)
        all_ids = features_df['read_id'].tolist()
        n = min(args.n_reads, len(all_ids))
        read_ids = random.sample(all_ids, n)

    print(f"Visualizing {len(read_ids)} reads from {args.chr_end}")

    # Build feature regions from TSV data
    read_data = []
    for rid in read_ids:
        feat_row = features_df[features_df['read_id'] == rid]
        if len(feat_row) == 0:
            print(f"  WARNING: Read {rid[:16]}... not found in features TSV")
            continue

        row = feat_row.iloc[0]
        feat_info = row.to_dict()
        read_length = int(row['read_length'])
        telo_side = row['telo_side']

        regions = build_feature_regions(row)
        if regions:
            read_data.append((rid, regions, read_length, telo_side, feat_info))
        else:
            print(f"  WARNING: No feature positions for read {rid[:16]}...")

    if not read_data:
        print("ERROR: No reads with feature positions found.")
        print("Make sure the features TSV contains spacer_start/end, x_element_start/end, y_prime_positions columns.")
        sys.exit(1)

    # Create figure
    n_reads = len(read_data)
    row_height = 1.2
    fig_height = max(3, n_reads * row_height + 1.5)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    max_read_len = max(rl for _, _, rl, _, _ in read_data)

    for i, (rid, regions, read_length, telo_side, feat_info) in enumerate(read_data):
        y_pos = (n_reads - 1 - i) * row_height
        draw_read(ax, regions, read_length, telo_side, rid, feat_info, y_pos=y_pos)

    # Configure axes
    ax.set_xlim(-max_read_len * 0.02, max_read_len * 1.02)
    ax.set_ylim(-0.8, n_reads * row_height + 0.3)
    ax.set_xlabel('Position on read (bp)', fontsize=10)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_title(f'Read Features — {args.chr_end} ({n_reads} reads)', fontsize=12, fontweight='bold')

    # Legend
    legend_handles = []
    for feat_type in ['anchor', 'spacer', 'x_element', 'y_prime', 'telomere']:
        display = FEATURE_DISPLAY_NAMES.get(feat_type, feat_type)
        legend_handles.append(mpatches.Patch(
            color=FEATURE_COLORS[feat_type], label=display, alpha=0.85))
    recomb_patch = mpatches.Patch(
        facecolor='white', edgecolor=RECOMB_BORDER_COLOR, linewidth=2,
        hatch=RECOMB_HATCH, label='Recombination', alpha=0.7)
    legend_handles.append(recomb_patch)

    ax.legend(handles=legend_handles, loc='upper right', fontsize=7,
              ncol=2, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(args.output, dpi=200, bbox_inches='tight')
    print(f"Saved to: {args.output}")


if __name__ == '__main__':
    main()
