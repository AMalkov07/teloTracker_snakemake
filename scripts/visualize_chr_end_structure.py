#!/usr/bin/env python3
"""
Visualize chromosome end structure from labeling pipeline output.

Creates ASCII diagrams and summary statistics showing the spatial relationship
between anchors, X primes, Y primes, and telomere ends for each chromosome end.

Supports two input formats:
1. TSV from labeling pipeline (pretelomeric_regions_*.tsv)
2. BED files with feature annotations (6991_final_features.bed format)
"""

import argparse
import pandas as pd
import sys
import re
from collections import defaultdict


def detect_file_format(file_path):
    """Detect if file is TSV (with header) or BED (no header)."""
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()

    # TSV has a header starting with 'chr_end'
    if first_line.startswith('chr_end\t'):
        return 'tsv'
    else:
        return 'bed'


def load_tsv(tsv_path):
    """Load the labeling TSV file."""
    df = pd.read_csv(tsv_path, sep='\t')
    return df


def load_bed(bed_path):
    """
    Load a BED file with feature annotations.

    Expected format (6991_final_features.bed):
    chr  start  end  name  strand  length

    Feature names contain chr end and feature type, e.g.:
    - chr1L_anchor
    - chr1L_x_core_element
    - chr1L_x_variable_element
    - chr1L_Y_Prime_1
    - chr1L_Telomere_Repeat
    - chr1L_space_between_anchor
    - ITS_chr4R_Y_Prime_0-1
    """
    df = pd.read_csv(bed_path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'name', 'strand', 'length'])
    return df


def parse_feature_name(name):
    """
    Parse feature name to extract chr_end and feature type.

    Returns: (chr_end, feature_type, extra_info)

    Examples:
        chr1L_anchor -> (chr1L, anchor, None)
        chr1L_x_core_element -> (chr1L, x_core, None)
        chr1L_Y_Prime_1 -> (chr1L, y_prime, 1)
        chr4R_Y_Prime_10 -> (chr4R, y_prime, 10)
        ITS_chr4R_Y_Prime_0-1 -> (chr4R, its, 0-1)
        chr1L_Telomere_Repeat -> (chr1L, telomere, None)
        chr1L_space_between_anchor -> (chr1L, space, None)
    """
    # Handle ITS (internal telomeric sequence) entries
    if name.startswith('ITS_'):
        # ITS_chr4R_Y_Prime_0-1
        match = re.match(r'ITS_(chr\d+[LR])_Y_Prime_(\S+)', name)
        if match:
            return match.group(1), 'its', match.group(2)
        return None, 'its', None

    # Standard format: chrXX_feature_name
    match = re.match(r'(chr\d+[LR])_(.+)', name)
    if not match:
        return None, None, None

    chr_end = match.group(1)
    feature = match.group(2)

    # Classify the feature
    if feature == 'anchor':
        return chr_end, 'anchor', None
    elif feature == 'x_core_element':
        return chr_end, 'x_core', None
    elif feature == 'x_variable_element':
        return chr_end, 'x_variable', None
    elif feature == 'x_prime':
        # Simplified BED format uses x_prime directly
        return chr_end, 'x_prime', None
    elif feature.startswith('Y_Prime_'):
        # Y_Prime_1, Y_Prime_10, etc.
        num = feature.replace('Y_Prime_', '')
        return chr_end, 'y_prime', num
    elif feature == 'Telomere_Repeat':
        return chr_end, 'telomere', None
    elif feature == 'space_between_anchor':
        return chr_end, 'space', None
    else:
        return chr_end, 'other', feature

    return None, None, None


def get_chr_end_data_from_bed(df):
    """
    Organize BED data by chromosome end.

    Returns dict: {chr_end: {'anchor': [...], 'x_prime': [...], 'y_prime': [...],
                             'telomere': [...], 'its': [...], 'x_core': [...], 'x_variable': [...]}}
    """
    chr_ends = defaultdict(lambda: {
        'anchor': [], 'x_prime': [], 'y_prime': [],
        'telomere': [], 'its': [], 'x_core': [], 'x_variable': [], 'space': []
    })

    for _, row in df.iterrows():
        chr_end, feature_type, extra = parse_feature_name(row['name'])

        if chr_end is None or feature_type is None:
            continue

        region_info = {
            'chr': row['chr'],
            'start': row['start'],
            'end': row['end'],
            'length': row['length'],
            'strand': row['strand'],
            'source': row['name'],
            'pident': 100.0,  # Reference annotations are 100% by definition
            'extra': extra,
        }

        chr_ends[chr_end][feature_type].append(region_info)

    # Combine x_core and x_variable into x_prime for consistency
    for chr_end, data in chr_ends.items():
        if data['x_core'] or data['x_variable']:
            # Calculate combined X prime region
            all_x = data['x_core'] + data['x_variable']
            if all_x:
                min_start = min(x['start'] for x in all_x)
                max_end = max(x['end'] for x in all_x)
                combined_length = max_end - min_start
                strand = all_x[0]['strand']

                data['x_prime'] = [{
                    'chr': all_x[0]['chr'],
                    'start': min_start,
                    'end': max_end,
                    'length': combined_length,
                    'strand': strand,
                    'source': f"{chr_end}_xprime",
                    'pident': 100.0,
                    'has_core': len(data['x_core']) > 0,
                    'has_variable': len(data['x_variable']) > 0,
                }]

    return chr_ends


def get_chr_end_data(df):
    """
    Organize data by chromosome end from TSV format.

    Returns dict: {chr_end: {'anchor': [...], 'x_prime': [...], 'y_prime': [...]}}
    """
    chr_ends = defaultdict(lambda: {'anchor': [], 'x_prime': [], 'y_prime': [],
                                     'telomere': [], 'its': [], 'x_core': [], 'x_variable': [], 'space': []})

    for _, row in df.iterrows():
        chr_end = row['chr_end']
        region_type = row['region_type']

        region_info = {
            'chr': row['chr'],
            'start': row['start'],
            'end': row['end'],
            'length': row['length'],
            'strand': row['strand'],
            'source': row['source'],
            'pident': row['pident'],
        }

        chr_ends[chr_end][region_type].append(region_info)

    return chr_ends


def determine_arm(chr_end):
    """Determine if this is L (left) or R (right) arm."""
    return chr_end[-1]


def visualize_chr_end(chr_end, data, chr_length=None, scale=100, is_bed_format=False):
    """
    Create ASCII visualization of chromosome end structure.

    For L arm: telomere is at position 0 (left side)
    For R arm: telomere is at chromosome end (right side)

    Args:
        chr_end: Chromosome end name (e.g., 'chr1L')
        data: Dict with 'anchor', 'x_prime', 'y_prime' lists
        chr_length: Optional chromosome length for R arm calculation
        scale: bp per character for visualization
        is_bed_format: If True, include telomere and ITS info from BED data
    """
    arm = determine_arm(chr_end)

    anchors = data['anchor']
    xprimes = data['x_prime']
    yprimes = data['y_prime']
    telomeres = data.get('telomere', [])
    its_regions = data.get('its', [])

    # Collect all regions
    all_regions = []
    for a in anchors:
        all_regions.append(('anchor', a['start'], a['end'], a['length'], a.get('pident', 0)))
    for x in xprimes:
        all_regions.append(('x_prime', x['start'], x['end'], x['length'], x.get('pident', 0)))
    for y in yprimes:
        all_regions.append(('y_prime', y['start'], y['end'], y['length'], y.get('pident', 0)))
    for t in telomeres:
        all_regions.append(('telomere', t['start'], t['end'], t['length'], 100.0))

    if not all_regions:
        return None

    # Sort by position
    all_regions.sort(key=lambda x: x[1])

    # Determine visualization range
    min_pos = min(r[1] for r in all_regions)
    max_pos = max(r[2] for r in all_regions)

    # For L arm, telomere is at 0, so min_pos should include 0
    # For R arm, we need to show distance to chromosome end
    if arm == 'L':
        vis_start = 0
        vis_end = max_pos + 1000  # Add some padding
        telomere_pos = 'left'
    else:
        vis_start = min_pos - 1000  # Add padding before anchor
        vis_end = max_pos + 1000
        telomere_pos = 'right'

    # Build the output
    lines = []
    lines.append(f"\n{'='*80}")
    lines.append(f"  {chr_end}")
    lines.append(f"{'='*80}")

    # Statistics
    lines.append(f"\n  Region Summary:")
    lines.append(f"  {'-'*40}")

    # Telomere info (if from BED format)
    if telomeres:
        t = telomeres[0]
        lines.append(f"  TELOMERE: {t['start']:>10,} - {t['end']:>10,}  ({t['length']:,} bp)")

    # Anchor info
    if anchors:
        a = anchors[0]
        lines.append(f"  ANCHOR:   {a['start']:>10,} - {a['end']:>10,}  ({a['length']:,} bp, {a['pident']:.1f}% identity)")
    else:
        lines.append(f"  ANCHOR:   NOT FOUND")

    # X prime info
    if xprimes:
        x = xprimes[0]
        source_note = f" (from {x['source']})" if x['source'] != f"{chr_end}_xprime" else ""
        # Check if this is from BED and has core/variable info
        if x.get('has_core') is not None:
            if x.get('has_core') and x.get('has_variable'):
                xprime_note = " (core+variable)"
            elif x.get('has_core'):
                xprime_note = " (core only)"
            else:
                xprime_note = " (variable only)"
            source_note = xprime_note
        lines.append(f"  X PRIME:  {x['start']:>10,} - {x['end']:>10,}  ({x['length']:,} bp, {x['pident']:.1f}% identity){source_note}")
    else:
        lines.append(f"  X PRIME:  NOT FOUND")

    # Y prime info
    if yprimes:
        lines.append(f"  Y PRIMES: {len(yprimes)} detected")
        for i, y in enumerate(sorted(yprimes, key=lambda x: x['start'])):
            # Extract Y prime type from source
            source = y['source']
            if '#' in source:
                ytype = source.split('#')[1].split('/')[0]  # Long or Short
            elif 'Y_Prime_' in source:
                # BED format: chr4R_Y_Prime_1
                ytype = y.get('extra', '?')  # Y prime number
            else:
                ytype = "?"
            lines.append(f"      [{i+1}]  {y['start']:>10,} - {y['end']:>10,}  ({y['length']:,} bp, {ytype})")
    else:
        lines.append(f"  Y PRIMES: NONE DETECTED")

    # ITS info (if from BED format)
    if its_regions:
        lines.append(f"  ITS:      {len(its_regions)} internal telomeric sequences")
        total_its_bp = sum(its['length'] for its in its_regions)
        lines.append(f"            Total: {total_its_bp:,} bp")

    # Distance calculations
    lines.append(f"\n  Distance Analysis:")
    lines.append(f"  {'-'*40}")

    if arm == 'L':
        # For L arm: telomere at 0
        # Structure should be: TELOMERE --- Y' --- X' --- ANCHOR --- (chromosome middle)

        # Find the element closest to telomere (smallest start position)
        closest_to_telo = None
        closest_dist = float('inf')

        if yprimes:
            first_yprime = min(yprimes, key=lambda x: x['start'])
            if first_yprime['start'] < closest_dist:
                closest_to_telo = ('y_prime', first_yprime)
                closest_dist = first_yprime['start']

        if xprimes:
            first_xprime = min(xprimes, key=lambda x: x['start'])
            if first_xprime['start'] < closest_dist:
                closest_to_telo = ('x_prime', first_xprime)
                closest_dist = first_xprime['start']

        if closest_to_telo:
            lines.append(f"  Telomere to first element: {closest_dist:,} bp")

        # Anchor to X prime distance
        if anchors and xprimes:
            anchor = anchors[0]
            xprime = xprimes[0]
            dist = anchor['start'] - xprime['end']
            lines.append(f"  X prime to Anchor:         {dist:,} bp")

        # X prime to Y prime distance
        if xprimes and yprimes:
            xprime = xprimes[0]
            # Find Y prime closest to X prime (should be just before it)
            closest_yprime = max([y for y in yprimes if y['end'] <= xprime['start']],
                                 key=lambda x: x['end'], default=None)
            if closest_yprime:
                dist = xprime['start'] - closest_yprime['end']
                lines.append(f"  Y prime to X prime:        {dist:,} bp")

        # Y prime spacing (if multiple)
        if len(yprimes) > 1:
            sorted_yprimes = sorted(yprimes, key=lambda x: x['start'])
            gaps = []
            for i in range(len(sorted_yprimes) - 1):
                gap = sorted_yprimes[i+1]['start'] - sorted_yprimes[i]['end']
                gaps.append(gap)
            lines.append(f"  Y prime gaps:              {', '.join(f'{g:,}' for g in gaps)} bp")

    else:  # R arm
        # For R arm: telomere at chromosome end
        # Structure should be: (chromosome middle) --- ANCHOR --- X' --- Y' --- TELOMERE

        # Find the element closest to telomere (largest end position)
        if yprimes:
            last_yprime = max(yprimes, key=lambda x: x['end'])
            # Estimate distance to telomere (we don't have exact chr length, so estimate)
            lines.append(f"  Last Y prime ends at:      {last_yprime['end']:,} bp")

        # Anchor to X prime distance
        if anchors and xprimes:
            anchor = anchors[0]
            xprime = xprimes[0]
            dist = xprime['start'] - anchor['end']
            lines.append(f"  Anchor to X prime:         {dist:,} bp")

        # X prime to Y prime distance
        if xprimes and yprimes:
            xprime = xprimes[0]
            # Find Y prime closest to X prime (should be just after it)
            closest_yprime = min([y for y in yprimes if y['start'] >= xprime['end']],
                                 key=lambda x: x['start'], default=None)
            if closest_yprime:
                dist = closest_yprime['start'] - xprime['end']
                lines.append(f"  X prime to Y prime:        {dist:,} bp")

        # Y prime spacing (if multiple)
        if len(yprimes) > 1:
            sorted_yprimes = sorted(yprimes, key=lambda x: x['start'])
            gaps = []
            for i in range(len(sorted_yprimes) - 1):
                gap = sorted_yprimes[i+1]['start'] - sorted_yprimes[i]['end']
                gaps.append(gap)
            lines.append(f"  Y prime gaps:              {', '.join(f'{g:,}' for g in gaps)} bp")

    # ASCII diagram
    lines.append(f"\n  Structure Diagram (not to scale):")
    lines.append(f"  {'-'*40}")

    # Build the diagram - always show from (middle)/anchor to telomere
    # For both L and R arms, we want: (middle)---[ANCHOR]---[X']---[Y']---[TELO]
    if arm == 'L':
        # L arm: elements are ordered telomere->anchor on chromosome
        # But we want to display anchor->telomere, so reverse the order
        diagram_parts = ["(middle)---"]

        # Sort all elements by position (descending - anchor first, telomere last)
        elements = []
        for t in telomeres:
            elements.append(('T', t['start'], t['end'], t['length']))
        for y in yprimes:
            elements.append(('Y', y['start'], y['end'], y['length']))
        if xprimes:
            x = xprimes[0]
            elements.append(('X', x['start'], x['end'], x['length']))
        if anchors:
            a = anchors[0]
            elements.append(('A', a['start'], a['end'], a['length']))

        # Sort by position descending (anchor has highest position, telomere lowest)
        elements.sort(key=lambda x: x[1], reverse=True)

        # Calculate gaps between adjacent elements (going from anchor toward telomere)
        prev_start = None
        for i, (elem_type, start, end, length) in enumerate(elements):
            if i > 0 and prev_start is not None:
                gap = prev_start - end
                if gap > 0:
                    gap_str = f"--{gap//1000}kb--" if gap >= 1000 else f"--{gap}bp--"
                    diagram_parts.append(gap_str)

            if elem_type == 'T':
                diagram_parts.append(f"[TELO:{length}bp]")
            elif elem_type == 'Y':
                diagram_parts.append(f"[Y':{length//1000}kb]")
            elif elem_type == 'X':
                diagram_parts.append(f"[X':{length}bp]")
            else:
                diagram_parts.append(f"[ANCHOR:{length//1000}kb]")

            prev_start = start

        # If no telomere annotation, end with TELO label
        if not telomeres:
            # Add gap from last element to position 0 (telomere)
            if prev_start is not None and prev_start > 0:
                gap_str = f"--{prev_start//1000}kb--" if prev_start >= 1000 else f"--{prev_start}bp--"
                diagram_parts.append(gap_str)
            diagram_parts.append("TELO")

    else:  # R arm
        # R arm: (middle) --- [ANCHOR] --- [X'] --- [Y'] --- [TELO]
        diagram_parts = ["(middle)---"]

        # Sort all elements by position
        elements = []
        if anchors:
            a = anchors[0]
            elements.append(('A', a['start'], a['end'], a['length']))
        if xprimes:
            x = xprimes[0]
            elements.append(('X', x['start'], x['end'], x['length']))
        for y in yprimes:
            elements.append(('Y', y['start'], y['end'], y['length']))
        for t in telomeres:
            elements.append(('T', t['start'], t['end'], t['length']))

        elements.sort(key=lambda x: x[1])

        prev_end = elements[0][1] if elements else 0
        for i, (elem_type, start, end, length) in enumerate(elements):
            if i > 0:
                gap = start - prev_end
                if gap > 0:
                    gap_str = f"--{gap//1000}kb--" if gap >= 1000 else f"--{gap}bp--"
                    diagram_parts.append(gap_str)

            if elem_type == 'T':
                diagram_parts.append(f"[TELO:{length}bp]")
            elif elem_type == 'Y':
                diagram_parts.append(f"[Y':{length//1000}kb]")
            elif elem_type == 'X':
                diagram_parts.append(f"[X':{length}bp]")
            else:
                diagram_parts.append(f"[ANCHOR:{length//1000}kb]")

            prev_end = end

        # If no telomere annotation, end with TELO label
        if not telomeres:
            diagram_parts.append("---TELO")

    diagram = "".join(diagram_parts)

    # Word wrap the diagram if too long
    max_width = 76
    if len(diagram) > max_width:
        # Try to break at -- boundaries
        wrapped = []
        current_line = "  "
        for part in diagram_parts:
            if len(current_line) + len(part) > max_width:
                wrapped.append(current_line)
                current_line = "      " + part  # Indent continuation
            else:
                current_line += part
        if current_line.strip():
            wrapped.append(current_line)
        for line in wrapped:
            lines.append(line)
    else:
        lines.append(f"  {diagram}")

    return "\n".join(lines)


def generate_summary_table(chr_ends, is_bed_format=False):
    """Generate a summary table of all chromosome ends."""
    lines = []
    lines.append("\n" + "="*110)
    lines.append("  SUMMARY TABLE")
    lines.append("="*110)

    if is_bed_format:
        header = f"{'Chr End':<10} {'Telo bp':<10} {'Anchor':<8} {'X Prime':<12} {'Y Primes':<10} {'Total Y bp':<12} {'ITS':<6} {'Structure Type':<20}"
    else:
        header = f"{'Chr End':<10} {'Anchor':<8} {'X Prime':<10} {'Y Primes':<10} {'Total Y bp':<12} {'Structure Type':<20}"
    lines.append(header)
    lines.append("-"*110)

    # Sort chr ends naturally
    def chr_sort_key(chr_end):
        # Extract chromosome number and arm
        match = re.match(r'chr(\d+)([LR])', chr_end)
        if match:
            return (int(match.group(1)), match.group(2))
        return (999, chr_end)

    for chr_end in sorted(chr_ends.keys(), key=chr_sort_key):
        data = chr_ends[chr_end]

        has_anchor = "Yes" if data['anchor'] else "No"
        has_xprime = "Yes" if data['x_prime'] else "No"
        n_yprimes = len(data['y_prime'])

        total_yprime_bp = sum(y['length'] for y in data['y_prime'])

        # Telomere and ITS info (for BED format)
        telo_bp = sum(t['length'] for t in data.get('telomere', []))
        n_its = len(data.get('its', []))

        # Determine structure type based on Y prime sizes
        if n_yprimes == 0:
            struct_type = "No Y'"
        elif n_yprimes == 1:
            y = data['y_prime'][0]
            # Check Y prime length to determine Long vs Short
            if 'Long' in y.get('source', ''):
                struct_type = "Solo Long Y'"
            elif 'Short' in y.get('source', ''):
                struct_type = "Solo Short Y'"
            elif y['length'] > 6000:
                struct_type = "Solo Long Y'"
            else:
                struct_type = "Solo Short Y'"
        else:
            # Check if tandem
            struct_type = f"Tandem ({n_yprimes}x Y')"

        # Add X prime notes
        if data['x_prime']:
            x = data['x_prime'][0]
            expected_source = f"{chr_end}_xprime"
            if x['source'] != expected_source and not x.get('has_core'):
                has_xprime = f"Yes*"  # Mark mismatch
            elif x.get('has_core') is not None:
                # BED format with core/variable info
                if x.get('has_core') and x.get('has_variable'):
                    has_xprime = "Yes (c+v)"
                elif x.get('has_core'):
                    has_xprime = "Yes (c)"
                else:
                    has_xprime = "Yes (v)"

        if is_bed_format:
            line = f"{chr_end:<10} {telo_bp:<10} {has_anchor:<8} {has_xprime:<12} {n_yprimes:<10} {total_yprime_bp:>10,}  {n_its:<6} {struct_type:<20}"
        else:
            line = f"{chr_end:<10} {has_anchor:<8} {has_xprime:<10} {n_yprimes:<10} {total_yprime_bp:>10,}  {struct_type:<20}"
        lines.append(line)

    lines.append("-"*110)
    if not is_bed_format:
        lines.append("  * = X prime source mismatch (best match from different chr end)")
    else:
        lines.append("  X Prime: (c+v) = core+variable, (c) = core only, (v) = variable only")

    # Statistics
    total_chr_ends = len(chr_ends)
    with_anchor = sum(1 for d in chr_ends.values() if d['anchor'])
    with_xprime = sum(1 for d in chr_ends.values() if d['x_prime'])
    with_yprime = sum(1 for d in chr_ends.values() if d['y_prime'])
    total_yprimes = sum(len(d['y_prime']) for d in chr_ends.values())

    lines.append(f"\n  Statistics:")
    lines.append(f"  - Total chromosome ends: {total_chr_ends}")
    lines.append(f"  - With anchor: {with_anchor}/{total_chr_ends}")
    lines.append(f"  - With X prime: {with_xprime}/{total_chr_ends}")
    lines.append(f"  - With Y prime(s): {with_yprime}/{total_chr_ends}")
    lines.append(f"  - Total Y prime regions: {total_yprimes}")

    if is_bed_format:
        total_telo_bp = sum(sum(t['length'] for t in d.get('telomere', [])) for d in chr_ends.values())
        total_its = sum(len(d.get('its', [])) for d in chr_ends.values())
        lines.append(f"  - Total telomere bp: {total_telo_bp:,}")
        lines.append(f"  - Total ITS regions: {total_its}")

    return "\n".join(lines)


def extract_xprime_source_chr_end(source_name):
    """
    Extract the chromosome end from an X prime source name.

    Examples:
        chr9L_xprime -> chr9L
        chr10R_xprime -> chr10R
        chr1L_x_prime -> chr1L
    """
    # Handle various formats
    if '_xprime' in source_name:
        return source_name.replace('_xprime', '')
    elif '_x_prime' in source_name:
        return source_name.replace('_x_prime', '')
    return None


def generate_comparison_report(new_data, ref_data, new_file, ref_file):
    """
    Generate a comparison report between two BED datasets.

    Compares:
    1. X prime source mismatches (when X prime comes from different chr end)
    2. Y prime count changes

    Args:
        new_data: Dict of chr_end data from newly created BED file
        ref_data: Dict of chr_end data from reference BED file
        new_file: Path to new BED file (for display)
        ref_file: Path to reference BED file (for display)

    Returns:
        String with comparison report
    """
    lines = []
    lines.append("\n" + "="*100)
    lines.append("  COMPARISON REPORT")
    lines.append("="*100)
    lines.append(f"  New file: {new_file}")
    lines.append(f"  Reference: {ref_file}")
    lines.append("="*100)

    # Sort chr ends naturally
    def chr_sort_key(chr_end):
        match = re.match(r'chr(\d+)([LR])', chr_end)
        if match:
            return (int(match.group(1)), match.group(2))
        return (999, chr_end)

    all_chr_ends = sorted(set(new_data.keys()) | set(ref_data.keys()), key=chr_sort_key)

    # Collect issues
    xprime_mismatches = []
    yprime_changes = []
    missing_in_new = []
    missing_in_ref = []

    for chr_end in all_chr_ends:
        new_chr = new_data.get(chr_end)
        ref_chr = ref_data.get(chr_end)

        # Check if chr end exists in both
        if new_chr is None:
            missing_in_new.append(chr_end)
            continue
        if ref_chr is None:
            missing_in_ref.append(chr_end)
            continue

        # Check X prime source mismatch
        if new_chr['x_prime']:
            x = new_chr['x_prime'][0]
            source = x.get('source', '')
            source_chr_end = extract_xprime_source_chr_end(source)

            if source_chr_end and source_chr_end != chr_end:
                xprime_mismatches.append({
                    'chr_end': chr_end,
                    'expected': chr_end,
                    'actual': source_chr_end,
                    'source': source
                })

        # Check Y prime count changes
        new_yprime_count = len(new_chr.get('y_prime', []))
        ref_yprime_count = len(ref_chr.get('y_prime', []))

        if new_yprime_count != ref_yprime_count:
            yprime_changes.append({
                'chr_end': chr_end,
                'new_count': new_yprime_count,
                'ref_count': ref_yprime_count,
                'diff': new_yprime_count - ref_yprime_count
            })

    # Report X prime source mismatches
    lines.append("\n" + "-"*100)
    lines.append("  X PRIME SOURCE MISMATCHES")
    lines.append("-"*100)

    if xprime_mismatches:
        lines.append(f"\n  Found {len(xprime_mismatches)} chromosome ends where X prime comes from a different chr end:")
        lines.append("")
        lines.append(f"  {'Chr End':<12} {'Expected X Prime':<20} {'Actual X Prime':<20} {'Status':<15}")
        lines.append(f"  {'-'*12} {'-'*20} {'-'*20} {'-'*15}")

        for m in xprime_mismatches:
            status = "MISMATCH"
            lines.append(f"  {m['chr_end']:<12} {m['expected']+'_xprime':<20} {m['actual']+'_xprime':<20} {status:<15}")
    else:
        lines.append("\n  No X prime source mismatches found. All X primes match their expected chromosome ends.")

    # Report Y prime count changes
    lines.append("\n" + "-"*100)
    lines.append("  Y PRIME COUNT CHANGES")
    lines.append("-"*100)

    if yprime_changes:
        lines.append(f"\n  Found {len(yprime_changes)} chromosome ends with different Y prime counts:")
        lines.append("")
        lines.append(f"  {'Chr End':<12} {'New Count':<12} {'Ref Count':<12} {'Change':<15} {'Status':<20}")
        lines.append(f"  {'-'*12} {'-'*12} {'-'*12} {'-'*15} {'-'*20}")

        for c in yprime_changes:
            diff = c['diff']
            if diff > 0:
                change_str = f"+{diff}"
                status = "GAINED Y prime(s)"
            else:
                change_str = str(diff)
                status = "LOST Y prime(s)"

            lines.append(f"  {c['chr_end']:<12} {c['new_count']:<12} {c['ref_count']:<12} {change_str:<15} {status:<20}")
    else:
        lines.append("\n  No Y prime count changes found. All chromosome ends have the same Y prime counts as reference.")

    # Report missing chr ends
    if missing_in_new or missing_in_ref:
        lines.append("\n" + "-"*100)
        lines.append("  MISSING CHROMOSOME ENDS")
        lines.append("-"*100)

        if missing_in_new:
            lines.append(f"\n  Missing in new file (present in reference): {', '.join(missing_in_new)}")
        if missing_in_ref:
            lines.append(f"\n  Missing in reference (present in new file): {', '.join(missing_in_ref)}")

    # Summary statistics
    lines.append("\n" + "-"*100)
    lines.append("  COMPARISON SUMMARY")
    lines.append("-"*100)

    total_new_yprimes = sum(len(d.get('y_prime', [])) for d in new_data.values())
    total_ref_yprimes = sum(len(d.get('y_prime', [])) for d in ref_data.values())
    total_new_xprimes = sum(1 for d in new_data.values() if d.get('x_prime'))
    total_ref_xprimes = sum(1 for d in ref_data.values() if d.get('x_prime'))

    lines.append(f"\n  Total chromosome ends: New={len(new_data)}, Reference={len(ref_data)}")
    lines.append(f"  Total Y primes: New={total_new_yprimes}, Reference={total_ref_yprimes} (diff: {total_new_yprimes - total_ref_yprimes:+d})")
    lines.append(f"  Total X primes detected: New={total_new_xprimes}, Reference={total_ref_xprimes}")
    lines.append(f"  X prime source mismatches: {len(xprime_mismatches)}")
    lines.append(f"  Y prime count changes: {len(yprime_changes)}")

    lines.append("\n" + "="*100)
    lines.append("  END OF COMPARISON REPORT")
    lines.append("="*100)

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Visualize chromosome end structure from labeling pipeline output (TSV or BED format)'
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input file: TSV from labeling pipeline (e.g., pretelomeric_regions_7575.tsv) '
             'or BED with feature annotations (e.g., 6991_final_features.bed)'
    )
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output file for visualization (default: stdout)'
    )
    parser.add_argument(
        '--chr-end', '-c',
        default=None,
        help='Visualize only this chromosome end (e.g., chr4R)'
    )
    parser.add_argument(
        '--summary-only', '-s',
        action='store_true',
        help='Only show summary table, skip individual diagrams'
    )
    parser.add_argument(
        '--format', '-f',
        choices=['auto', 'tsv', 'bed'],
        default='auto',
        help='Input file format (default: auto-detect)'
    )
    parser.add_argument(
        '--compare',
        default=None,
        help='Reference BED file to compare against (e.g., references/6991_final_features.bed). '
             'When provided, generates a comparison report showing X prime mismatches and Y prime count changes.'
    )

    args = parser.parse_args()

    # Detect or use specified format
    if args.format == 'auto':
        file_format = detect_file_format(args.input)
    else:
        file_format = args.format

    # Load data based on format
    print(f"Loading {args.input} (format: {file_format})...", file=sys.stderr)

    if file_format == 'tsv':
        df = load_tsv(args.input)
        chr_ends = get_chr_end_data(df)
        is_bed_format = False
    else:  # bed
        df = load_bed(args.input)
        chr_ends = get_chr_end_data_from_bed(df)
        is_bed_format = True

    print(f"Found {len(chr_ends)} chromosome ends", file=sys.stderr)

    # Load comparison reference if provided
    ref_chr_ends = None
    if args.compare:
        print(f"Loading reference file for comparison: {args.compare}...", file=sys.stderr)
        ref_format = detect_file_format(args.compare)
        if ref_format == 'tsv':
            ref_df = load_tsv(args.compare)
            ref_chr_ends = get_chr_end_data(ref_df)
        else:
            ref_df = load_bed(args.compare)
            ref_chr_ends = get_chr_end_data_from_bed(ref_df)
        print(f"Found {len(ref_chr_ends)} chromosome ends in reference", file=sys.stderr)

    # Generate output
    output_lines = []

    output_lines.append("="*80)
    output_lines.append("  CHROMOSOME END STRUCTURE VISUALIZATION")
    output_lines.append(f"  Input: {args.input}")
    output_lines.append(f"  Format: {file_format.upper()}")
    if args.compare:
        output_lines.append(f"  Comparison: {args.compare}")
    output_lines.append("="*80)

    # Summary table first
    output_lines.append(generate_summary_table(chr_ends, is_bed_format=is_bed_format))

    # Comparison report if reference provided
    if ref_chr_ends:
        output_lines.append(generate_comparison_report(chr_ends, ref_chr_ends, args.input, args.compare))

    # Individual diagrams
    if not args.summary_only:
        if args.chr_end:
            # Single chromosome end
            if args.chr_end in chr_ends:
                vis = visualize_chr_end(args.chr_end, chr_ends[args.chr_end], is_bed_format=is_bed_format)
                if vis:
                    output_lines.append(vis)
            else:
                print(f"ERROR: Chromosome end '{args.chr_end}' not found", file=sys.stderr)
                sys.exit(1)
        else:
            # All chromosome ends
            def chr_sort_key(chr_end):
                match = re.match(r'chr(\d+)([LR])', chr_end)
                if match:
                    return (int(match.group(1)), match.group(2))
                return (999, chr_end)

            output_lines.append("\n" + "="*80)
            output_lines.append("  INDIVIDUAL CHROMOSOME END DETAILS")
            output_lines.append("="*80)

            for chr_end in sorted(chr_ends.keys(), key=chr_sort_key):
                vis = visualize_chr_end(chr_end, chr_ends[chr_end], is_bed_format=is_bed_format)
                if vis:
                    output_lines.append(vis)

    output_lines.append("\n" + "="*80)
    output_lines.append("  END OF REPORT")
    output_lines.append("="*80 + "\n")

    # Output
    output_text = "\n".join(output_lines)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_text)
        print(f"Output written to: {args.output}", file=sys.stderr)
    else:
        print(output_text)


if __name__ == '__main__':
    main()
