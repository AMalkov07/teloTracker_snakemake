#!/usr/bin/env python3
"""
plot_read_structure.py — Per-read telomeric structure diagram.

Draws a schematic like the paper figure: left-pointing arrows for Y' elements,
small squares for inter-Y' spacers with bp labels, X element, spacer, and
chromosome end label.

Usage (specify exact features TSV):
    python scripts/plot_read_structure.py \
        --features results/<base>/recombination/<base>_<chr_end>_features.tsv \
        --read-id <read_id> \
        --output diagram.png

Usage (auto-search across all chr_ends):
    python scripts/plot_read_structure.py \
        --recombination-dir results/<base>/recombination/ \
        --base-name <base> \
        --read-id <read_id> \
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


# ─────────────────────────────── Color palette ────────────────────────────────

# Y' type → fill color (pink family like reference figure)
Y_PRIME_COLORS: dict[str, str] = {
    'ID1':     '#F2AACF',   # pink
    'ID2':     '#C9A8E0',   # light purple
    'ID3':     '#A8C9E0',   # light blue
    'ID4':     '#A8E0C9',   # light teal
    'ID5':     '#E0C9A8',   # light orange
    'ID6':     '#E0A8A8',   # light red
    'default': '#F2AACF',
}

INTER_YP_COLOR = '#1A1A1A'   # near-black squares between Y' elements
X_ELEMENT_COLOR = '#444444'  # dark grey
SPACER_COLOR    = '#7399AB'  # muted teal-blue
ANCHOR_COLOR    = '#282828'  # very dark
RAIL_COLOR      = '#BBBBBB'  # horizontal connecting line


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

def draw_left_arrow(ax, x_left: float, x_right: float, y: float,
                    height: float, facecolor: str, edgecolor: str = 'black',
                    lw: float = 0.8) -> None:
    """Pentagon arrow pointing left (tip at x_left, tail at x_right)."""
    hw = height / 2
    head_w = min(hw * 1.6, (x_right - x_left) * 0.30)
    verts = [
        (x_left,              y),
        (x_left + head_w,     y + hw),
        (x_right,             y + hw),
        (x_right,             y - hw),
        (x_left + head_w,     y - hw),
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

def plot_read_structure(row: dict, read_id: str, output_path: str) -> None:
    """Render the schematic structure diagram for one read and save to output_path."""

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
    total_width = x + 1.2  # extra room for chr label text

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
            color   = Y_PRIME_COLORS.get(yp_type, Y_PRIME_COLORS['default'])
            edge    = '#C0392B' if is_div else 'black'
            lw      = 1.8       if is_div else 0.8
            draw_left_arrow(ax, xl, xr, YC, YP_H, color, edge, lw)

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

    # Chromosome label at far right
    ax.text(total_width - 1.05, YC, chr_end,
            ha='left', va='center', fontsize=11, fontweight='bold')

    # Telomere indicator at far left
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
        color = Y_PRIME_COLORS.get(t, Y_PRIME_COLORS['default'])
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

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Per-read telomeric structure diagram',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('--features', metavar='TSV',
                     help='Path to a single <base>_<chr_end>_features.tsv file')
    src.add_argument('--recombination-dir', metavar='DIR',
                     help='Dir containing *_features.tsv (requires --base-name)')
    p.add_argument('--base-name', metavar='NAME',
                   help='Sample base name (used with --recombination-dir)')
    p.add_argument('--read-id', required=True, metavar='ID',
                   help='Read ID to plot')
    p.add_argument('--output', required=True, metavar='PNG',
                   help='Output PNG file path')
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.recombination_dir and not args.base_name:
        sys.exit('Error: --base-name is required when using --recombination-dir')

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

    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plot_read_structure(row, args.read_id, args.output)


if __name__ == '__main__':
    main()
