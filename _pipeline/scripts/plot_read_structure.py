#!/usr/bin/env python3
"""
plot_read_structure.py — Per-read telomeric structure diagram.

Draws a schematic like the paper figure: left-pointing arrows for Y' elements,
small squares for inter-Y' spacers with bp labels, X element, spacer, and
chromosome end label.

Y' colors are read from the Y' library FASTA (IDn_ColorName encoded in headers),
matching the colors assigned during label_regions.sh clustering. Supply the FASTA
with --yprime-lib, or let the script find it automatically via --config (reads
the Snakemake config.yaml and resolves y_prime_lib / y_prime_lib_override).

Usage (specify exact features TSV + Y' library):
    python scripts/plot_read_structure.py \
        --features results/<base>/recombination/<base>_<chr_end>_features.tsv \
        --read-id <read_id> \
        --config config.yaml \
        --output diagram.png

Usage (auto-search across all chr_ends, explicit Y' library):
    python scripts/plot_read_structure.py \
        --recombination-dir results/<base>/recombination/ \
        --base-name <base> \
        --read-id <read_id> \
        --yprime-lib references/extracted_yprimes_6991.fasta \
        --output diagram.png
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

import pandas as pd


# ─────────────────────────── Color name → hex mapping ─────────────────────────
# Maps the color names embedded in Y' FASTA headers (IDn_ColorName) to hex.
# Covers all names produced by cluster_yprimes_paper_method.py, including
# -Light / -Dark / -Neutral shade variants.

COLOR_NAME_TO_HEX: dict[str, str] = {
    'Gray':            '#909090',
    'Red':             '#E05555',
    'Red-Light':       '#F1948A',
    'Red-Dark':        '#C0392B',
    'Green':           '#5CBF7A',
    'Green-Light':     '#82E0AA',
    'Green-Dark':      '#1E8449',
    'Orange':          '#E8943A',
    'Orange-Light':    '#F0B27A',
    'Orange-Dark':     '#CA6F1E',
    'Purple':          '#A569BD',
    'Purple-Light':    '#C39BD3',
    'Purple-Dark':     '#7D3C98',
    'Purple-Neutral':  '#9B59B6',
    'Blue':            '#5DADE2',
    'Blue-Light':      '#85C1E9',
    'Blue-Dark':       '#1F618D',
    'Yellow':          '#F4D03F',
    'Yellow-Light':    '#F9E79F',
    'Yellow-Dark':     '#D4AC0D',
    'Cyan':            '#48C9B0',
    'Cyan-Light':      '#76D7C4',
    'Cyan-Dark':       '#117A65',
    'Magenta':         '#E874A8',
    'Magenta-Light':   '#F1A9C7',
    'Magenta-Dark':    '#C0396B',
    'Brown':           '#A0522D',
    'Brown-Light':     '#C49A6C',
    'Brown-Dark':      '#6B3619',
    'Pink':            '#F48FB1',
    'Pink-Light':      '#F8BBD0',
    'Pink-Dark':       '#E06080',
    'Teal':            '#26A69A',
    'Teal-Light':      '#80CBC4',
    'Teal-Dark':       '#00695C',
    'Olive':           '#8D9440',
    'Navy':            '#1A5276',
    'Coral':           '#F1755A',
    'Lavender':        '#CE93D8',
    'Maroon':          '#880E4F',
    'Gold':            '#F5C518',
    'Lime':            '#9CCC65',
    'Slate':           '#607D8B',
}

# Fallback hardcoded colors when no Y' FASTA is provided
_FALLBACK_COLORS: dict[str, str] = {
    'ID1':     '#F2AACF',
    'ID2':     '#C9A8E0',
    'ID3':     '#A8C9E0',
    'ID4':     '#A8E0C9',
    'ID5':     '#E0C9A8',
    'ID6':     '#E0A8A8',
    'default': '#F2AACF',
}

INTER_YP_COLOR  = '#1A1A1A'
X_ELEMENT_COLOR = '#444444'
SPACER_COLOR    = '#7399AB'
ANCHOR_COLOR    = '#282828'
RAIL_COLOR      = '#BBBBBB'


# ──────────────────────── Y' library color loader ─────────────────────────────

def load_yprime_colors(fasta_path: str) -> dict[str, str]:
    """
    Parse a Y' FASTA library and return {IDn: hex_color}.

    Headers have format:  >...#Size/Type/IDn_ColorName
    e.g.  >Y_Prime_chr5R1#Long/Solo/ID1_Gray
          >Y_Prime_chr4R1#Long/Tandem/ID2_Red-Light

    The first occurrence of each IDn wins (subsequent sequences in the same
    cluster may have a different shade suffix but share the base color).
    """
    id_to_hex: dict[str, str] = {}
    try:
        with open(fasta_path) as fh:
            for line in fh:
                if not line.startswith('>'):
                    continue
                header = line[1:].strip()
                if '#' not in header:
                    continue
                metadata = header.split('#', 1)[1]   # "Size/Type/IDn_ColorName"
                parts = metadata.split('/')
                if len(parts) < 3:
                    continue
                id_color = parts[2].strip()           # e.g. "ID2_Red-Light"
                if '_' not in id_color:
                    continue
                id_str, color_name = id_color.split('_', 1)   # "ID2", "Red-Light"
                if id_str in id_to_hex:
                    continue  # first occurrence wins
                hex_color = COLOR_NAME_TO_HEX.get(color_name)
                if hex_color is None:
                    # Try base color (strip -Light/-Dark/-Neutral suffix)
                    base = color_name.split('-')[0]
                    hex_color = COLOR_NAME_TO_HEX.get(base, '#F2AACF')
                id_to_hex[id_str] = hex_color
    except OSError as exc:
        print(f'  Warning: could not read Y\' library {fasta_path}: {exc}')
    return id_to_hex


def resolve_yprime_lib(config_path: str) -> str | None:
    """
    Read a Snakemake config.yaml and return the resolved y_prime_lib path.
    Respects y_prime_lib_override when present.
    Path is resolved relative to the config file's directory.
    """
    try:
        import yaml
    except ImportError:
        print('  Warning: PyYAML not available — cannot auto-read config. '
              'Install pyyaml or use --yprime-lib directly.')
        return None

    try:
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
    except Exception as exc:
        print(f'  Warning: could not read config {config_path}: {exc}')
        return None

    refs = cfg.get('references', {})
    override = refs.get('y_prime_lib_override', '')
    lib_template = refs.get('y_prime_lib', '')
    strain = str(cfg.get('strain', ''))

    raw_path = override if override else lib_template.replace('{strain}', strain)
    if not raw_path:
        print('  Warning: no y_prime_lib found in config')
        return None

    # Resolve relative to the config file's directory
    config_dir = os.path.dirname(os.path.abspath(config_path))
    return os.path.join(config_dir, raw_path) if not os.path.isabs(raw_path) else raw_path


# ─────────────────────────────── Parsers ──────────────────────────────────────

def parse_y_prime_positions(raw) -> list[tuple[str, int, int]]:
    """
    Parse 'ID1:start-end;ID2:start-end;...' into [(type, start, end), ...].
    Ordering in the string is innermost-first (chromosome side first).
    """
    if not raw or (isinstance(raw, float) and math.isnan(raw)):
        return []
    result = []
    for part in str(raw).split(';'):
        part = part.strip()
        if ':' not in part or '-' not in part:
            continue
        yp_type, coords = part.split(':', 1)
        try:
            s, e = (int(x) for x in coords.split('-', 1))
            result.append((yp_type.strip(), s, e))
        except ValueError:
            continue
    return result


def safe_int(v, default: int = -1) -> int:
    try:
        return int(v)
    except (ValueError, TypeError):
        return default


# ─────────────────────────────── Read lookup ──────────────────────────────────

def find_read(read_id: str,
              features_tsv: str | None,
              recomb_dir: str | None,
              base_name: str | None) -> tuple[dict | None, str | None]:
    """Return (row_dict, source_file) for the given read_id, or (None, None)."""

    def search_file(fpath: str) -> dict | None:
        try:
            df = pd.read_csv(fpath, sep='\t', dtype=str)
        except Exception as exc:
            print(f'  Warning: could not read {fpath}: {exc}')
            return None
        if 'read_id' not in df.columns:
            return None
        hits = df[df['read_id'] == read_id]
        if not hits.empty:
            return hits.iloc[0].to_dict()
        return None

    if features_tsv:
        row = search_file(features_tsv)
        if row:
            return row, features_tsv
        print(f'  Read {read_id!r} not found in {features_tsv}')
        return None, None

    if recomb_dir and base_name:
        pattern = os.path.join(recomb_dir, f'{base_name}_*_features.tsv')
        files = sorted(glob.glob(pattern))
        if not files:
            print(f'  No features files found: {pattern}')
            return None, None
        for fpath in files:
            row = search_file(fpath)
            if row:
                return row, fpath
        print(f'  Read {read_id!r} not found in any features file under {recomb_dir}')
        return None, None

    print('  Error: provide --features or --recombination-dir + --base-name')
    return None, None


# ─────────────────────────────── Display coords ───────────────────────────────

def to_display(raw_pos: int, telo_side: str, read_length: int) -> int:
    """Convert raw read position → display position (chromosome always on right)."""
    if telo_side == 'beginning':
        # Telomere at raw 0, chromosome at raw read_length → already chromosome-right
        return raw_pos
    # telo_side == 'end': anchor at raw 0, telomere at raw read_length → flip
    return read_length - raw_pos


def get_display_elements(row: dict) -> dict:
    """
    Extract all drawable elements in display coordinates (chromosome on right,
    telomere on left). Returns a dict with keys:
        telo_side, read_length,
        anchor (start, end) | None,
        spacer (start, end) | None, spacer_size, spacer_source,
        x_element (start, end) | None, x_element_size, x_element_source,
        y_primes: [(type, left, right), ...] sorted left→right (telomere→chromosome),
        chr_end, y_prime_recombination_status, y_prime_divergence_idx,
        recombination_source.
    """
    telo_side = str(row.get('telo_side', 'end')).strip().lower()
    rl = safe_int(row.get('read_length', 0), 0)

    def d(col: str) -> int:
        return to_display(safe_int(row.get(col, -1)), telo_side, rl)

    def span(start_col: str, end_col: str) -> tuple[int, int] | None:
        s, e = d(start_col), d(end_col)
        if s < 0 or e < 0:
            return None
        lo, hi = min(s, e), max(s, e)
        return (lo, hi) if hi > lo else None

    anchor    = span('anchor_start',    'anchor_end')
    spacer    = span('spacer_start',    'spacer_end')
    x_element = span('x_element_start', 'x_element_end')

    # Y' elements: parse, convert, sort left→right (outermost/telomere first)
    raw_yps = parse_y_prime_positions(row.get('y_prime_positions', ''))
    yps: list[tuple[str, int, int]] = []
    for yp_type, raw_s, raw_e in raw_yps:
        d_s = to_display(raw_s, telo_side, rl)
        d_e = to_display(raw_e, telo_side, rl)
        lo, hi = min(d_s, d_e), max(d_s, d_e)
        if hi > lo:
            yps.append((yp_type, lo, hi))
    yps.sort(key=lambda x: x[1])   # left = telomere side

    divergence_idx = safe_int(row.get('y_prime_divergence_idx', -1))

    return {
        'telo_side':                   telo_side,
        'read_length':                 rl,
        'anchor':                      anchor,
        'spacer':                      spacer,
        'spacer_size':                 safe_int(row.get('spacer_size', -1)),
        'spacer_source':               str(row.get('spacer_source', '')),
        'x_element':                   x_element,
        'x_element_size':              safe_int(row.get('x_element_size', -1)),
        'x_element_source':            str(row.get('x_element_source', '')),
        'y_primes':                    yps,
        'chr_end':                     str(row.get('chr_end', 'unknown')),
        'y_prime_recombination_status': str(row.get('y_prime_recombination_status', '')),
        'y_prime_divergence_idx':      divergence_idx,
        'recombination_source':        str(row.get('recombination_source', '')),
        'recombination_detected':      str(row.get('recombination_detected', '')).lower() == 'true',
    }


# ─────────────────────────────── Drawing helpers ──────────────────────────────

def draw_arrow(ax, x_left: float, x_right: float, y: float,
               height: float, facecolor: str, edgecolor: str = 'black',
               lw: float = 0.8, direction: str = 'left') -> None:
    """Pentagon arrow. direction='left': tip at x_left; direction='right': tip at x_right."""
    hw = height / 2
    head_w = min(hw * 1.6, (x_right - x_left) * 0.30)
    if direction == 'left':
        verts = [
            (x_left,              y),
            (x_left + head_w,     y + hw),
            (x_right,             y + hw),
            (x_right,             y - hw),
            (x_left + head_w,     y - hw),
        ]
    else:
        verts = [
            (x_right,             y),
            (x_right - head_w,    y + hw),
            (x_left,              y + hw),
            (x_left,              y - hw),
            (x_right - head_w,    y - hw),
        ]
    ax.add_patch(Polygon(verts, closed=True,
                          facecolor=facecolor, edgecolor=edgecolor,
                          linewidth=lw, zorder=3))


def draw_box(ax, x_left: float, x_right: float, y: float,
             height: float, facecolor: str, edgecolor: str = 'black',
             lw: float = 0.6) -> None:
    ax.add_patch(Rectangle(
        (x_left, y - height / 2), x_right - x_left, height,
        facecolor=facecolor, edgecolor=edgecolor, linewidth=lw, zorder=3))


# ─────────────────────────────── Main plot ────────────────────────────────────

def plot_read_structure(row: dict, read_id: str, output_path: str,
                        yprime_colors: dict[str, str] | None = None) -> None:
    """Render the schematic structure diagram for one read and save to output_path.

    yprime_colors: {IDn: hex_color} from load_yprime_colors(). Falls back to
    hardcoded _FALLBACK_COLORS when None.
    """
    colors = yprime_colors if yprime_colors else _FALLBACK_COLORS

    def yp_color(id_str: str) -> str:
        return colors.get(id_str, colors.get('default', '#F2AACF'))

    elems = get_display_elements(row)
    yps   = elems['y_primes']
    n_yp  = len(yps)

    if n_yp == 0:
        print(f'  No Y\' elements in y_prime_positions for {read_id} — nothing to plot')
        return

    chr_end = elems['chr_end']
    div_idx = elems['y_prime_divergence_idx']   # 0-based from innermost; -1 = none

    # ── Schematic layout units ──────────────────────────────────────────────
    YP_W      = 2.5    # Y' arrow width
    YP_H      = 0.55   # Y' arrow height
    IYP_W     = 0.38   # inter-Y' spacer square width
    IYP_H     = 0.38   # inter-Y' spacer square height
    XE_W      = 0.45   # X element box width
    XE_H      = 0.45
    SP_W      = 0.60   # spacer box width
    SP_H      = 0.45
    ANC_W     = 0.38   # anchor box width
    ANC_H     = 0.45
    PAD       = 0.18   # gap between adjacent elements
    LABEL_GAP = 0.08   # extra gap between element edge and text
    YC        = 0.0    # vertical center of all elements

    # ── Build element list left→right (telomere→chromosome) ────────────────
    # Each entry: (kind, x_left, x_right, info_dict)
    elements: list[tuple[str, float, float, dict]] = []
    x = 0.5   # start with left padding (for telomere label)

    # Y' elements and inter-Y' spacers
    for i, (yp_type, yp_left_bp, yp_right_bp) in enumerate(yps):
        # Y' label: Y-N for outermost, Y-1 for innermost
        yp_label = f"Y-{n_yp - i}"
        # Divergence marker: the divergence_idx is 0-based from innermost.
        # Innermost = yps[-1] = index n_yp-1. So the divergent element in display order
        # is at display index = (n_yp - 1 - div_idx).
        is_divergence = (div_idx >= 0 and i == (n_yp - 1 - div_idx))
        elements.append(('y_prime', x, x + YP_W, {
            'label': yp_label,
            'type': yp_type,
            'divergence': is_divergence,
        }))
        x += YP_W

        # Inter-Y' spacer (between this Y' and the next one to the right)
        if i < n_yp - 1:
            _, next_left_bp, _ = yps[i + 1]
            gap_bp = max(0, next_left_bp - yp_right_bp)
            x += PAD
            elements.append(('inter_yp', x, x + IYP_W, {'size_bp': gap_bp}))
            x += IYP_W
            x += PAD

    # X element
    if elems['x_element'] is not None:
        x += PAD
        elements.append(('x_element', x, x + XE_W, {
            'size_bp': elems['x_element_size'],
            'source':  elems['x_element_source'],
        }))
        x += XE_W

    # Spacer (between X element and anchor)
    if elems['spacer'] is not None:
        x += PAD
        elements.append(('spacer', x, x + SP_W, {
            'size_bp': elems['spacer_size'],
            'source':  elems['spacer_source'],
        }))
        x += SP_W

    # Anchor
    if elems['anchor'] is not None:
        x += PAD
        elements.append(('anchor', x, x + ANC_W, {}))
        x += ANC_W

    x += PAD  # trailing padding before chromosome label
    x_content_end = x  # remember end of content before adding label space

    # ── R-arm orientation: mirror elements so chromosome is on the left ─────
    # L arms: telomere ← [Y-N … Y-1] [X] [spacer] [anchor]   chr
    # R arms: chr   [anchor] [spacer] [X] [Y-1 … Y-N] →  telomere
    is_r = str(chr_end).upper().endswith('R')
    if is_r:
        x_start = 0.5  # same leading offset used when building
        elements = [
            (kind, x_start + x_content_end - xr,
                   x_start + x_content_end - xl, info)
            for kind, xl, xr, info in reversed(elements)
        ]
        arrow_dir = 'right'
    else:
        arrow_dir = 'left'

    total_width = x_content_end + 1.2  # extra room for chr label text

    # ── Figure setup ───────────────────────────────────────────────────────
    fig_w = max(9.0, total_width * 1.15)
    fig_h = 3.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(-0.8, total_width)
    ax.set_ylim(-1.6, 1.6)
    ax.axis('off')

    # Background rail
    ax.hlines(YC, 0.0, total_width - 1.1,
              colors=RAIL_COLOR, linewidth=1.5, zorder=1)

    # ── Draw each element ───────────────────────────────────────────────────
    above = YP_H / 2 + LABEL_GAP    # y-offset for labels above elements
    below = -(YP_H / 2 + LABEL_GAP) # y-offset for labels below elements
    sq_below = -(IYP_H / 2 + LABEL_GAP)

    for kind, xl, xr, info in elements:
        mid = (xl + xr) / 2

        if kind == 'y_prime':
            yp_type = info['type']
            label   = info['label']
            is_div  = info['divergence']
            color   = yp_color(yp_type)
            edge    = '#C0392B' if is_div else 'black'
            lw      = 1.8       if is_div else 0.8
            draw_arrow(ax, xl, xr, YC, YP_H, color, edge, lw, direction=arrow_dir)

            # Label above
            ax.text(mid, YC + above, label,
                    ha='center', va='bottom', fontsize=9, fontweight='bold')
            # Y' type below (smaller)
            ax.text(mid, YC + below, yp_type,
                    ha='center', va='top', fontsize=7, color='#555555')
            # Mark divergence point
            if is_div:
                ax.text(mid, YC + above + 0.20, '★',
                        ha='center', va='bottom', fontsize=10, color='#C0392B')

        elif kind == 'inter_yp':
            gap_bp = info['size_bp']
            draw_box(ax, xl, xr, YC, IYP_H, INTER_YP_COLOR)
            if gap_bp > 0:
                ax.text(mid, YC + sq_below, str(gap_bp),
                        ha='center', va='top', fontsize=7, color='#333333')

        elif kind == 'x_element':
            draw_box(ax, xl, xr, YC, XE_H, X_ELEMENT_COLOR)
            src = info['source']
            bp  = info['size_bp']
            label_lines = ['X']
            if src and src not in ('', 'nan'):
                label_lines.append(src)
            ax.text(mid, YC + above, '\n'.join(label_lines),
                    ha='center', va='bottom', fontsize=6.5, linespacing=1.2)
            if bp > 0:
                ax.text(mid, YC + sq_below, f'{bp} bp',
                        ha='center', va='top', fontsize=6.5, color='#333333')

        elif kind == 'spacer':
            draw_box(ax, xl, xr, YC, SP_H, SPACER_COLOR)
            src = info['source']
            bp  = info['size_bp']
            label_lines = ['spacer']
            if src and src not in ('', 'nan'):
                label_lines.append(src)
            ax.text(mid, YC + above, '\n'.join(label_lines),
                    ha='center', va='bottom', fontsize=6.5, linespacing=1.2)
            if bp > 0:
                ax.text(mid, YC + sq_below, f'{bp} bp',
                        ha='center', va='top', fontsize=6.5, color='#333333')

        elif kind == 'anchor':
            draw_box(ax, xl, xr, YC, ANC_H, ANCHOR_COLOR)
            ax.text(mid, YC + above, 'anc',
                    ha='center', va='bottom', fontsize=6.5, color='#333333')

    if is_r:
        # Chromosome label at far left, telomere indicator at far right
        ax.text(0.40, YC, chr_end,
                ha='right', va='center', fontsize=11, fontweight='bold')
        ax.annotate('', xy=(total_width - 0.15, YC), xytext=(total_width - 0.65, YC),
                    arrowprops=dict(arrowstyle='->', color='#888888', lw=1.5))
        ax.text(total_width - 0.10, YC, 'telo', ha='left', va='center',
                fontsize=8, color='#888888')
    else:
        # Chromosome label at far right, telomere indicator at far left
        ax.text(total_width - 1.05, YC, chr_end,
                ha='left', va='center', fontsize=11, fontweight='bold')
        ax.annotate('', xy=(-0.55, YC), xytext=(-0.1, YC),
                    arrowprops=dict(arrowstyle='->', color='#888888', lw=1.5))
        ax.text(-0.60, YC, 'telo', ha='right', va='center',
                fontsize=8, color='#888888')

    # ── Title ──────────────────────────────────────────────────────────────
    status = elems['y_prime_recombination_status']
    src    = elems['recombination_source']
    rid_display = read_id if len(read_id) <= 40 else read_id[:38] + '…'
    title_parts = [f'{rid_display}']
    title_parts.append(f'{chr_end}  ·  {n_yp} Y\' elements  ·  {status}')
    if src and src not in ('', 'nan'):
        title_parts.append(f'source: {src}')
    ax.set_title('\n'.join(title_parts), fontsize=9, pad=8, loc='left')

    # ── Legend ─────────────────────────────────────────────────────────────
    seen_types = {info['type'] for kind, _, _, info in elements if kind == 'y_prime'}
    legend_handles = []
    for t in sorted(seen_types):
        color = yp_color(t)
        legend_handles.append(
            mpatches.Patch(facecolor=color, edgecolor='black', linewidth=0.5, label=f"Y' {t}"))
    legend_handles.append(
        mpatches.Patch(facecolor=INTER_YP_COLOR, label='inter-Y\' spacer'))
    if any(k == 'x_element' for k, *_ in elements):
        legend_handles.append(
            mpatches.Patch(facecolor=X_ELEMENT_COLOR, label='X element'))
    if any(k == 'spacer' for k, *_ in elements):
        legend_handles.append(
            mpatches.Patch(facecolor=SPACER_COLOR, label='spacer'))
    legend_handles.append(
        mpatches.Patch(facecolor=ANCHOR_COLOR, label='anchor'))
    if div_idx >= 0:
        legend_handles.append(
            mpatches.Patch(facecolor='white', edgecolor='#C0392B', linewidth=1.5,
                           label='★ divergence point'))

    ax.legend(handles=legend_handles, loc='lower center',
              ncol=min(len(legend_handles), 6),
              bbox_to_anchor=(0.45, -0.05),
              fontsize=7, frameon=False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {output_path}')


# ─────────────────────────────── CLI ──────────────────────────────────────────

# ─────────────────────────────── Read listing ─────────────────────────────────

def list_reads(features_tsv: str | None,
               recomb_dir: str | None,
               base_name: str | None,
               chr_end_filter: str | None,
               recombination_only: bool,
               min_yprimes: int) -> None:
    """Print a table of available reads that can be plotted."""

    # Collect all features files to scan
    if features_tsv:
        files = [features_tsv]
    elif recomb_dir and base_name:
        pattern = os.path.join(recomb_dir, f'{base_name}_*_features.tsv')
        files = sorted(glob.glob(pattern))
        if not files:
            sys.exit(f'No features files found: {pattern}')
    else:
        sys.exit('Provide --features or --recombination-dir + --base-name')

    rows = []
    for fpath in files:
        try:
            df = pd.read_csv(fpath, sep='\t', dtype=str)
        except Exception as exc:
            print(f'  Warning: could not read {fpath}: {exc}')
            continue
        if df.empty:
            continue

        if chr_end_filter and 'chr_end' in df.columns:
            df = df[df['chr_end'] == chr_end_filter]

        if recombination_only and 'recombination_detected' in df.columns:
            df = df[df['recombination_detected'].str.lower() == 'true']

        if min_yprimes > 0 and 'y_prime_count_on_read' in df.columns:
            df = df[pd.to_numeric(df['y_prime_count_on_read'], errors='coerce').fillna(0) >= min_yprimes]

        rows.append(df)

    if not rows:
        print('No reads match the given filters.')
        return

    all_df = pd.concat(rows, ignore_index=True)

    display_cols = {
        'chr_end':                     'chr_end',
        'read_id':                     'read_id',
        'y_prime_count_on_read':       'y_primes',
        'y_prime_recombination_status':'status',
        'recombination_source':        'source',
        'recombination_detected':      'recomb',
        'telo_side':                   'telo_side',
    }
    present = [c for c in display_cols if c in all_df.columns]
    out = all_df[present].rename(columns=display_cols)
    out = out.fillna('').sort_values(['chr_end', 'read_id'] if 'chr_end' in out.columns else ['read_id'])

    # Column widths
    widths = {col: max(len(col), out[col].astype(str).str.len().max()) for col in out.columns}
    header = '  '.join(col.ljust(widths[col]) for col in out.columns)
    sep    = '  '.join('-' * widths[col] for col in out.columns)

    print(f'\n{len(out)} read(s) available\n')
    print(header)
    print(sep)
    for _, r in out.iterrows():
        print('  '.join(str(r[col]).ljust(widths[col]) for col in out.columns))


# ─────────────────────────────── CLI ──────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Per-read telomeric structure diagram',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Data source (always required)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('--features', metavar='TSV',
                     help='Path to a single <base>_<chr_end>_features.tsv file')
    src.add_argument('--recombination-dir', metavar='DIR',
                     help='Dir containing *_features.tsv (requires --base-name)')
    p.add_argument('--base-name', metavar='NAME',
                   help='Sample base name (used with --recombination-dir)')

    # Mode: list reads or plot one
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--list-reads', action='store_true',
                      help='Print available reads instead of plotting')
    mode.add_argument('--read-id', metavar='ID',
                      help='Read ID to plot')

    p.add_argument('--output', metavar='PNG',
                   help='Output PNG path (required with --read-id)')

    # Filters for --list-reads
    p.add_argument('--chr-end', metavar='CHR_END',
                   help='Only show reads from this chr_end (e.g. chr10R)')
    p.add_argument('--recombination-only', action='store_true',
                   help='Only show reads where recombination_detected is True')
    p.add_argument('--min-yprimes', type=int, default=0, metavar='N',
                   help='Only show reads with >= N Y\' elements on the read')

    # Y' color source
    col = p.add_mutually_exclusive_group()
    col.add_argument('--yprime-lib', metavar='FASTA',
                     help='Y\' library FASTA with IDn_ColorName headers')
    col.add_argument('--config', metavar='YAML',
                     help='Snakemake config.yaml to auto-resolve y_prime_lib')

    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.recombination_dir and not args.base_name:
        sys.exit('Error: --base-name is required when using --recombination-dir')

    # ── List mode ──────────────────────────────────────────────────────────
    if args.list_reads:
        list_reads(
            features_tsv=args.features,
            recomb_dir=args.recombination_dir,
            base_name=args.base_name,
            chr_end_filter=args.chr_end,
            recombination_only=args.recombination_only,
            min_yprimes=args.min_yprimes,
        )
        return

    # ── Plot mode ──────────────────────────────────────────────────────────
    if not args.output:
        sys.exit('Error: --output is required when plotting a read')

    print(f'plot_read_structure.py — read: {args.read_id}')

    row, source_file = find_read(
        args.read_id,
        features_tsv=args.features,
        recomb_dir=args.recombination_dir,
        base_name=args.base_name,
    )
    if row is None:
        sys.exit(1)

    print(f'  Found in: {source_file}')
    print(f'  chr_end={row.get("chr_end")}  '
          f'telo_side={row.get("telo_side")}  '
          f'y_prime_count={row.get("y_prime_count_on_read")}  '
          f'status={row.get("y_prime_recombination_status")}')

    # Resolve Y' colors
    yprime_colors: dict[str, str] | None = None
    if args.yprime_lib:
        yprime_colors = load_yprime_colors(args.yprime_lib)
        print(f'  Y\' colors from: {args.yprime_lib}  '
              f'({len(yprime_colors)} IDs: {", ".join(sorted(yprime_colors))})')
    elif args.config:
        lib_path = resolve_yprime_lib(args.config)
        if lib_path:
            print(f'  Y\' lib resolved from config: {lib_path}')
            yprime_colors = load_yprime_colors(lib_path)
            print(f'  Y\' colors loaded: {len(yprime_colors)} IDs: '
                  f'{", ".join(sorted(yprime_colors))}')
    if yprime_colors is None:
        print('  Y\' colors: using fallback defaults (no --yprime-lib / --config provided)')

    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plot_read_structure(row, args.read_id, args.output, yprime_colors=yprime_colors)


if __name__ == '__main__':
    main()
