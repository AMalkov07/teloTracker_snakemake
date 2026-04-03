#!/usr/bin/env python3
"""
Visualize reads from the recombination pipeline, showing features as colored boxes
on a read-length line graph.

Features (spacer, X element, ITS, Y prime, telomere) are mapped from the day0
reference BED onto each read via BAM alignment coordinates.

Usage:
    # Single read by ID
    python scripts/visualize_reads.py \
        --recomb-dir recombination/ \
        --chr-end chr14R \
        --read-id 721f2dad-6d96-4bb8-85cd-82b227e96b36 \
        --bed results/.../pretelomeric_regions_7302_simp.bed \
        -o output.png

    # 10 random reads from a chr_end
    python scripts/visualize_reads.py \
        --recomb-dir recombination/ \
        --chr-end chr14R \
        --n-reads 10 \
        --bed results/.../pretelomeric_regions_7302_simp.bed \
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
import numpy as np
import pandas as pd
import pysam


# Feature colors
FEATURE_COLORS = {
    'anchor': '#4A4A4A',
    'space_between_anchor': '#D3D3D3',
    'x_core_element': '#2196F3',
    'x_variable_element': '#64B5F6',
    'y_prime': '#E53935',
    'ITS': '#FFA726',
    'Telomere_Repeat': '#4CAF50',
    'telomere_unmapped': '#81C784',
}

FEATURE_DISPLAY_NAMES = {
    'anchor': 'Anchor',
    'space_between_anchor': 'Spacer',
    'x_core_element': 'X core',
    'x_variable_element': 'X variable',
    'y_prime': "Y'",
    'ITS': 'ITS',
    'Telomere_Repeat': 'Telomere',
    'telomere_unmapped': 'Telomere',
}


def parse_args():
    p = argparse.ArgumentParser(description='Visualize recombination pipeline reads')
    p.add_argument('--recomb-dir', required=True, help='Recombination output directory')
    p.add_argument('--chr-end', required=True, help='Chromosome end (e.g., chr14R)')
    p.add_argument('--bed', required=True, help='Pretelomeric regions BED file')
    p.add_argument('--read-id', default=None, help='Specific read ID to visualize')
    p.add_argument('--n-reads', type=int, default=10, help='Number of random reads (if no --read-id)')
    p.add_argument('-o', '--output', default='read_visualization.png', help='Output image path')
    p.add_argument('--seed', type=int, default=42, help='Random seed for read selection')
    return p.parse_args()


def load_bed_features(bed_path, chr_end):
    """Load features for a chr_end from BED file."""
    features = []
    with open(bed_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            feat_name = fields[3]
            if not feat_name.startswith(chr_end + '_'):
                # Also match ITS entries like ITS_chr14R_Y_Prime_0-1
                if not (feat_name.startswith('ITS_') and chr_end in feat_name):
                    continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[4]
            length = int(fields[5])

            # Determine feature type
            if 'anchor' in feat_name and 'space' not in feat_name:
                feat_type = 'anchor'
            elif 'space_between_anchor' in feat_name:
                feat_type = 'space_between_anchor'
            elif 'x_core_element' in feat_name:
                feat_type = 'x_core_element'
            elif 'x_variable_element' in feat_name:
                feat_type = 'x_variable_element'
            elif 'Y_Prime' in feat_name and 'ITS' not in feat_name:
                feat_type = 'y_prime'
            elif 'ITS' in feat_name:
                feat_type = 'ITS'
            elif 'Telomere_Repeat' in feat_name:
                feat_type = 'Telomere_Repeat'
            else:
                feat_type = feat_name

            features.append({
                'chrom': chrom, 'start': start, 'end': end,
                'name': feat_name, 'strand': strand, 'length': length,
                'type': feat_type
            })

    return features


def map_features_to_read(bam_path, read_id, features, chr_end):
    """Map reference feature coordinates to read coordinates using BAM alignment.

    Returns list of (read_start, read_end, feature_dict) tuples.
    """
    # Find the reference contig name
    chrom_num = re.match(r'chr(\d+)', chr_end).group(1)

    mapped = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        target_ref = None
        for ref_name in bam.references:
            if f'chr{chrom_num}' in ref_name:
                target_ref = ref_name
                break

        if target_ref is None:
            return [], 0, 'end'

        for read in bam.fetch():
            if read.query_name != read_id:
                continue
            if read.is_unmapped:
                continue

            read_length = read.query_length or len(read.query_sequence)

            # Determine telo_side from alignment orientation
            # For R-arm: telomere at high coords, anchor at low
            # For L-arm: telomere at low coords, anchor at high
            arm = chr_end[-1]

            # Get aligned pairs for this read on the target reference
            # Try primary and supplementary alignments
            all_pairs = []
            if read.reference_name == target_ref:
                pairs = read.get_aligned_pairs()
                all_pairs.extend(pairs)

            if not all_pairs:
                continue

            # Build ref_pos -> read_pos mapping
            ref_to_read = {}
            for read_pos, ref_pos in all_pairs:
                if read_pos is not None and ref_pos is not None:
                    ref_to_read[ref_pos] = read_pos

            if not ref_to_read:
                continue

            # Determine telo_side
            ref_positions = sorted(ref_to_read.keys())
            read_at_min_ref = ref_to_read[ref_positions[0]]
            read_at_max_ref = ref_to_read[ref_positions[-1]]
            if arm == 'R':
                # R-arm: anchor at low ref coords
                # If read pos increases with ref pos, telo_side=end
                telo_side = 'end' if read_at_max_ref > read_at_min_ref else 'beginning'
            else:
                # L-arm: anchor at high ref coords
                telo_side = 'beginning' if read_at_max_ref > read_at_min_ref else 'end'

            # Map each feature
            for feat in features:
                if feat['chrom'] != target_ref.replace('_extended', ''):
                    # Try matching with _extended suffix
                    if feat['chrom'] + '_extended' != target_ref:
                        continue

                # Find read positions for feature start and end
                feat_start_read = None
                feat_end_read = None

                # Search for closest mapped position to feature boundaries
                for ref_pos in range(feat['start'], feat['end'] + 1):
                    if ref_pos in ref_to_read:
                        if feat_start_read is None:
                            feat_start_read = ref_to_read[ref_pos]
                        feat_end_read = ref_to_read[ref_pos]

                if feat_start_read is not None and feat_end_read is not None:
                    rstart = min(feat_start_read, feat_end_read)
                    rend = max(feat_start_read, feat_end_read)
                    if rend - rstart >= 5:  # minimum size
                        mapped.append((rstart, rend, feat))

            # Add telomere for the unmapped region at the telomere side
            if telo_side == 'end':
                # Telomere at the end of the read
                last_feature_end = max((e for _, e, _ in mapped), default=0)
                if last_feature_end < read_length - 50:
                    mapped.append((last_feature_end, read_length, {
                        'name': 'Telomere', 'type': 'telomere_unmapped',
                        'length': read_length - last_feature_end
                    }))
            else:
                # Telomere at the beginning of the read
                first_feature_start = min((s for s, _, _ in mapped), default=read_length)
                if first_feature_start > 50:
                    mapped.append((0, first_feature_start, {
                        'name': 'Telomere', 'type': 'telomere_unmapped',
                        'length': first_feature_start
                    }))

            return mapped, read_length, telo_side

    return [], 0, 'end'


def get_read_ids(recomb_dir, chr_end, base_pattern=None):
    """Get all read IDs for a chr_end from the features TSV."""
    # Find the features file
    features_files = [f for f in os.listdir(recomb_dir)
                      if f.endswith(f'_{chr_end}_features.tsv')]
    if not features_files:
        print(f"ERROR: No features file found for {chr_end} in {recomb_dir}")
        sys.exit(1)

    features_path = os.path.join(recomb_dir, features_files[0])
    df = pd.read_csv(features_path, sep='\t')
    return df


RECOMB_BORDER_COLOR = '#FFD600'  # yellow border for recombination
RECOMB_HATCH = '///'


def _feature_has_recombination(feat_type, chr_end, features_info):
    """Check if a feature shows recombination based on the pipeline results.

    Returns (is_recomb, source_label):
        is_recomb: True if this feature came from a different end or changed
        source_label: short string describing the source (e.g., "chr9L") or None
    """
    if feat_type == 'space_between_anchor':
        recomb = features_info.get('spacer_recombination', 'no_change')
        source = features_info.get('spacer_source', chr_end)
        if recomb in ('full_switch', 'switch_detected') and source != chr_end:
            return True, source
    elif feat_type in ('x_core_element', 'x_variable_element'):
        recomb = features_info.get('x_element_recombination', 'no_change')
        source = features_info.get('x_element_source', chr_end)
        if recomb in ('full_switch', 'switch_detected') and source != chr_end:
            return True, source
    elif feat_type == 'y_prime':
        yp_status = features_info.get('y_prime_recombination_status', 'No Change')
        if yp_status not in ('No Change', ''):
            compatible = features_info.get('y_prime_compatible_ends', '')
            return True, compatible if compatible else yp_status
    return False, None


def draw_read(ax, mapped_features, read_length, telo_side, read_id, chr_end, features_info,
              y_pos=0, height=0.6):
    """Draw a single read with feature boxes on the given axes.

    Features with detected recombination get a thick yellow border and
    diagonal hatching, plus a label showing the suspected source.
    """

    # Draw the read backbone line
    ax.plot([0, read_length], [y_pos, y_pos], color='#888888', linewidth=1, zorder=1)

    # Draw feature boxes
    for rstart, rend, feat in sorted(mapped_features, key=lambda x: x[0]):
        feat_type = feat['type']
        color = FEATURE_COLORS.get(feat_type, '#9E9E9E')
        display_name = FEATURE_DISPLAY_NAMES.get(feat_type, feat_type)
        width = rend - rstart

        # Check if this feature shows recombination
        is_recomb, recomb_source = _feature_has_recombination(feat_type, chr_end, features_info)

        edgecolor = RECOMB_BORDER_COLOR if is_recomb else 'black'
        linewidth = 2.5 if is_recomb else 0.5

        box = FancyBboxPatch(
            (rstart, y_pos - height / 2), width, height,
            boxstyle="round,pad=0",
            facecolor=color, edgecolor=edgecolor, linewidth=linewidth,
            alpha=0.85, zorder=2
        )
        ax.add_patch(box)

        # Add hatching overlay for recombined features
        if is_recomb:
            hatch_box = FancyBboxPatch(
                (rstart, y_pos - height / 2), width, height,
                boxstyle="round,pad=0",
                facecolor='none', edgecolor=edgecolor, linewidth=0.5,
                hatch=RECOMB_HATCH, alpha=0.4, zorder=2.5
            )
            ax.add_patch(hatch_box)

        # Label inside the box if wide enough
        if width > read_length * 0.03:
            label = display_name
            # Add feature name details for Y primes
            if feat_type == 'y_prime':
                m = re.search(r'Y_Prime_(\d+)', feat.get('name', ''))
                if m:
                    label = f"Y'{m.group(1)}"
            elif feat_type == 'ITS':
                label = f"ITS ({width}bp)"

            fontsize = 7 if width > read_length * 0.06 else 5.5
            ax.text(rstart + width / 2, y_pos, label,
                    ha='center', va='center', fontsize=fontsize,
                    color='white' if feat_type not in ('space_between_anchor', 'ITS', 'telomere_unmapped') else 'black',
                    fontweight='bold', zorder=3)

        # Source annotation below recombined features
        if is_recomb and width > read_length * 0.02 and recomb_source:
            ax.text(rstart + width / 2, y_pos - height / 2 - 0.08, f'→{recomb_source}',
                    ha='center', va='top', fontsize=5, color='#B8860B',
                    fontweight='bold', fontstyle='italic', zorder=3)

    # Read ID and info label
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

    # Load features TSV for this chr_end
    features_df = get_read_ids(args.recomb_dir, args.chr_end)

    # Load reference BED features
    bed_features = load_bed_features(args.bed, args.chr_end)
    if not bed_features:
        print(f"ERROR: No features found for {args.chr_end} in BED file")
        sys.exit(1)

    print(f"Loaded {len(bed_features)} reference features for {args.chr_end}")

    # Find BAM file
    bam_files = [f for f in os.listdir(args.recomb_dir)
                 if f.endswith(f'_{args.chr_end}.bam')]
    if not bam_files:
        print(f"ERROR: No BAM file found for {args.chr_end} in {args.recomb_dir}")
        sys.exit(1)
    bam_path = os.path.join(args.recomb_dir, bam_files[0])

    # Select reads
    if args.read_id:
        read_ids = [args.read_id]
    else:
        random.seed(args.seed)
        all_ids = features_df['read_id'].tolist()
        n = min(args.n_reads, len(all_ids))
        read_ids = random.sample(all_ids, n)

    print(f"Visualizing {len(read_ids)} reads from {args.chr_end}")

    # Map features and collect data for each read
    read_data = []
    for rid in read_ids:
        mapped, read_length, telo_side = map_features_to_read(
            bam_path, rid, bed_features, args.chr_end)

        feat_row = features_df[features_df['read_id'] == rid]
        if len(feat_row) > 0:
            feat_info = feat_row.iloc[0].to_dict()
        else:
            feat_info = {}

        if mapped:
            read_data.append((rid, mapped, read_length, telo_side, feat_info))
        else:
            print(f"  WARNING: Could not map features for read {rid[:16]}...")

    if not read_data:
        print("ERROR: No reads could be mapped. Check BAM file and BED coordinates.")
        sys.exit(1)

    # Create figure
    n_reads = len(read_data)
    row_height = 1.2
    fig_height = max(3, n_reads * row_height + 1.5)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    max_read_len = max(rl for _, _, rl, _, _ in read_data)

    for i, (rid, mapped, read_length, telo_side, feat_info) in enumerate(read_data):
        y_pos = (n_reads - 1 - i) * row_height
        draw_read(ax, mapped, read_length, telo_side, rid, args.chr_end, feat_info, y_pos=y_pos)

    # Configure axes
    ax.set_xlim(-max_read_len * 0.02, max_read_len * 1.02)
    ax.set_ylim(-0.8, n_reads * row_height + 0.3)
    ax.set_xlabel('Position on read (bp)', fontsize=10)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Title
    ax.set_title(f'Read Features — {args.chr_end} ({n_reads} reads)', fontsize=12, fontweight='bold')

    # Legend
    legend_handles = []
    for feat_type in ['anchor', 'space_between_anchor', 'x_core_element', 'x_variable_element',
                       'ITS', 'y_prime', 'Telomere_Repeat', 'telomere_unmapped']:
        if feat_type in FEATURE_COLORS:
            display = FEATURE_DISPLAY_NAMES.get(feat_type, feat_type)
            legend_handles.append(mpatches.Patch(
                color=FEATURE_COLORS[feat_type], label=display, alpha=0.85))

    # Add recombination indicator to legend
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
