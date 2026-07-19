#!/usr/bin/env python3
"""
plot_reference_structure.py — Per-chr_end reference structure diagrams.

Renders a schematic of each chromosome end's pretelomeric structure based on
the labels.bed produced by label_regions.sh.  This is the reference-side
analogue of plot_read_structure.py: the same visual conventions (Y' arrows
with per-ID colors, inter-Y' spacer squares, X element box, spacer box,
anchor box, telomere box) but drawn for the day-0 reference annotations
rather than per-read features.

Inputs:
  --bed         pretelomeric_regions_<strain>_simp.bed (required)
  --yprime-lib  extracted_yprimes_<strain>.fasta (for ID->color mapping)
  --reference   reference FASTA (optional; used to compute gap to telomere
                for each chr_end, which makes assembly truncations visible)
  --output-dir  where to write the PNGs

Outputs:
  <output_dir>/reference_overview.png       — all 32 chr_ends in a single grid
  <output_dir>/per_chr_end/<chr_end>.png    — one PNG per chr_end (detail view)

Auto-detection: if --yprime-lib is omitted, the script looks for it next to
the input bed (the standard label_regions layout).

Usage (intended invocation from label_regions.sh):
    python _pipeline/scripts/plot_reference_structure.py \\
        --bed   ${BED_FILE} \\
        --yprime-lib ${EXTRACTED_YPRIMES} \\
        --reference  ${REFERENCE_FASTA} \\
        --output-dir ${OUTPUT_DIR}/reference_structure_figures
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle


# ─────────────────────────── Color tables (mirror plot_read_structure.py) ─────

COLOR_NAME_TO_HEX: dict[str, str] = {
    'Gray':           '#909090',
    'Red':            '#E05555', 'Red-Light':    '#F1948A', 'Red-Dark':    '#C0392B',
    'Green':          '#5CBF7A', 'Green-Light':  '#82E0AA', 'Green-Dark':  '#1E8449',
    'Orange':         '#E8943A', 'Orange-Light': '#F0B27A', 'Orange-Dark': '#CA6F1E',
    'Purple':         '#A569BD', 'Purple-Light': '#C39BD3', 'Purple-Dark': '#7D3C98',
    'Purple-Neutral': '#9B59B6',
    'Blue':           '#5DADE2', 'Blue-Light':   '#85C1E9', 'Blue-Dark':   '#1F618D',
    'Yellow':         '#F4D03F', 'Yellow-Light': '#F9E79F', 'Yellow-Dark': '#D4AC0D',
    'Cyan':           '#48C9B0', 'Cyan-Light':   '#76D7C4', 'Cyan-Dark':   '#117A65',
    'Magenta':        '#E874A8', 'Magenta-Light':'#F1A9C7', 'Magenta-Dark':'#C0396B',
    'Brown':          '#A0522D', 'Brown-Light':  '#C49A6C', 'Brown-Dark':  '#6B3619',
    'Pink':           '#F48FB1', 'Pink-Light':   '#F8BBD0', 'Pink-Dark':   '#E06080',
    'Teal':           '#26A69A', 'Teal-Light':   '#80CBC4',
}
_FALLBACK_COLORS = {'ID1': '#909090', 'ID2': '#E05555', 'ID3': '#5CBF7A',
                    'ID4': '#E8943A', 'ID5': '#A569BD', 'ID6': '#5DADE2',
                    'ID7': '#F4D03F', 'default': '#F2AACF'}

INTER_YP_COLOR  = '#1A1A1A'
X_ELEMENT_COLOR = '#444444'
SPACER_COLOR    = '#7399AB'
ANCHOR_COLOR    = '#282828'
RAIL_COLOR      = '#BBBBBB'
TELO_COLOR      = '#1A1A1A'


# ─────────────────────────── FASTA / library parsing ──────────────────────────

def parse_yprime_lib(fasta_path: str) -> tuple[dict[str, str], dict[tuple[str, int], str]]:
    """Read the extracted Y' FASTA.  Returns:
      id_to_hex:           {IDn: '#RRGGBB'}                — color per cluster ID
      yprime_id_lookup:    {(chr_end, yp_index): 'IDn'}    — chr_end+index → cluster ID
    Headers look like: >Y_Prime_chr5R1#Short/Solo/ID1_Gray
                                 ^chr_end ^idx          ^id_color
    """
    id_to_hex: dict[str, str] = {}
    lookup: dict[tuple[str, int], str] = {}
    if not fasta_path or not os.path.exists(fasta_path):
        return id_to_hex, lookup
    # Header forms encountered:
    #   Y_Prime_chr14L1#...                  (Solo)
    #   Y_Prime_chr14L2;chr16L2,3#...        (Tandem — primary listed first)
    # Match the primary chr_end+idx prefix only.
    name_pat = re.compile(r'^Y_Prime_(chr\d+[LR])(\d+)')
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            header = line[1:].strip()
            if '#' not in header:
                continue
            seq_id, meta = header.split('#', 1)
            parts = meta.split('/')
            if len(parts) < 3:
                continue
            id_color = parts[2].strip()
            if '_' not in id_color:
                continue
            yp_id, color_name = id_color.split('_', 1)
            if yp_id not in id_to_hex:
                hex_color = COLOR_NAME_TO_HEX.get(color_name)
                if hex_color is None:
                    base = color_name.split('-')[0]
                    hex_color = COLOR_NAME_TO_HEX.get(base, '#F2AACF')
                id_to_hex[yp_id] = hex_color
            m = name_pat.match(seq_id)
            if m:
                chr_end, idx = m.group(1), int(m.group(2))
                lookup[(chr_end, idx)] = yp_id
    return id_to_hex, lookup


def parse_contig_lengths(fasta_path: str) -> dict[str, int]:
    """Return {contig_name: length_bp} from a FASTA."""
    if not fasta_path or not os.path.exists(fasta_path):
        return {}
    lens: dict[str, int] = {}
    cur, cur_len = None, 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur is not None:
                    lens[cur] = cur_len
                cur = line[1:].split()[0]
                cur_len = 0
            else:
                cur_len += len(line.rstrip())
        if cur is not None:
            lens[cur] = cur_len
    # Also add a stripped-of-"_extended" alias since simp.bed often uses chr1
    # while the reference FASTA uses chr1_extended.
    for name in list(lens.keys()):
        stripped = name.replace('_extended', '')
        if stripped not in lens:
            lens[stripped] = lens[name]
    return lens


# ─────────────────────────── Bed parsing ──────────────────────────────────────

def parse_bed(bed_path: str) -> dict[str, dict]:
    """Group features by chr_end.  Returns:
       { chr_end: { 'contig': str,
                    'anchor':   (s, e) | None,
                    'x_var':    (s, e) | None,
                    'x_core':   (s, e) | None,
                    'spacer':   (s, e) | None,
                    'yprimes':  [(idx, s, e)]      sorted by idx,
                    'its':      [(i, j, s, e)]     inter-Y' spacers
                  }
       }
    """
    out: dict[str, dict] = {}
    yprime_pat = re.compile(r'^(chr\d+[LR])_Y_Prime_(\d+)$')
    its_pat    = re.compile(r'^ITS_(chr\d+[LR])_Y_Prime_(\d+)-(\d+)$')
    simple_pat = re.compile(r'^(chr\d+[LR])_(anchor|x_variable_element|x_core_element|space_between_anchor)$')

    with open(bed_path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 4:
                continue
            contig = cols[0]
            s, e = int(cols[1]), int(cols[2])
            name = cols[3]
            m = yprime_pat.match(name)
            if m:
                ce, idx = m.group(1), int(m.group(2))
                rec = out.setdefault(ce, _empty_rec())
                rec['contig'] = contig
                rec['yprimes'].append((idx, s, e))
                continue
            m = its_pat.match(name)
            if m:
                ce = m.group(1)
                i, j = int(m.group(2)), int(m.group(3))
                rec = out.setdefault(ce, _empty_rec())
                rec['contig'] = contig
                rec['its'].append((i, j, s, e))
                continue
            m = simple_pat.match(name)
            if m:
                ce, kind = m.group(1), m.group(2)
                key = {'anchor': 'anchor', 'x_variable_element': 'x_var',
                       'x_core_element': 'x_core',
                       'space_between_anchor': 'spacer'}[kind]
                rec = out.setdefault(ce, _empty_rec())
                rec['contig'] = contig
                rec[key] = (s, e)

    # Sort Y' lists
    for ce in out:
        out[ce]['yprimes'].sort(key=lambda t: t[0])
    return out


def _empty_rec() -> dict:
    return {'contig': None, 'anchor': None, 'x_var': None, 'x_core': None,
            'spacer': None, 'yprimes': [], 'its': []}


# ─────────────────────────── Drawing primitives ───────────────────────────────

def draw_arrow(ax, x_left: float, x_right: float, y: float, height: float,
               facecolor: str, edgecolor='black', lw=0.8, direction='left'):
    hw = height / 2
    head_w = min(hw * 1.6, (x_right - x_left) * 0.30)
    if direction == 'left':
        verts = [
            (x_left,           y),
            (x_left + head_w,  y + hw),
            (x_right,          y + hw),
            (x_right,          y - hw),
            (x_left + head_w,  y - hw),
        ]
    else:
        verts = [
            (x_right,          y),
            (x_right - head_w, y + hw),
            (x_left,           y + hw),
            (x_left,           y - hw),
            (x_right - head_w, y - hw),
        ]
    ax.add_patch(Polygon(verts, closed=True, facecolor=facecolor,
                          edgecolor=edgecolor, linewidth=lw, zorder=3))


def draw_box(ax, x_left: float, x_right: float, y: float, height: float,
             facecolor: str, edgecolor='black', lw=0.6):
    ax.add_patch(Rectangle(
        (x_left, y - height / 2), x_right - x_left, height,
        facecolor=facecolor, edgecolor=edgecolor, linewidth=lw, zorder=3))


# ─────────────────────────── One-chr_end rendering ────────────────────────────

# Schematic units — match plot_read_structure.py so the two outputs look
# visually consistent if opened side by side.
YP_W, YP_H       = 2.5,  0.55
IYP_W, IYP_H     = 0.38, 0.38
XE_W, XE_H       = 0.45, 0.45
SP_W, SP_H       = 0.60, 0.45
ANC_W, ANC_H     = 0.38, 0.45
TELO_W, TELO_H   = 0.55, 0.55
PAD              = 0.18
LABEL_GAP        = 0.08
YC               = 0.0


def yp_color(yp_id: str, id_to_hex: dict) -> str:
    if id_to_hex and yp_id in id_to_hex:
        return id_to_hex[yp_id]
    return _FALLBACK_COLORS.get(yp_id, _FALLBACK_COLORS['default'])


def build_elements(rec: dict, chr_end: str,
                   yprime_id_lookup: dict) -> tuple[list[tuple], int, int]:
    """Convert a chr_end's bed entries into the schematic element list.

    Returns (elements, n_yp, gap_to_telo_bp).
    elements layout: telomere on the LEFT of the schematic for both L and R
    arms; the caller mirrors at draw time if needed.
    """
    side = chr_end[-1]
    yps = rec['yprimes']   # sorted by idx ascending
    n_yp = len(yps)

    # Pick X-element: prefer x_var, fall back to x_core if x_var missing
    x_elem = rec.get('x_var') or rec.get('x_core')

    # Sort Y' on the contig: which Y' is OUTERMOST (telomere side)?
    # In the bed, Y_Prime_1 is innermost (anchor side) and Y_Prime_N is
    # outermost — but verify from the contig coordinates rather than
    # trusting the numbering, because numbering convention can flip.
    if n_yp > 0:
        # Position of each Y' (use start coord)
        yp_positions = [(idx, s, e) for idx, s, e in yps]
        # For L arms: outermost = smallest coord; for R arms: largest coord
        if side == 'L':
            yp_positions.sort(key=lambda t: t[1])  # telomere-side first
        else:
            yp_positions.sort(key=lambda t: -t[1])
    else:
        yp_positions = []

    elements: list[tuple] = []
    x = 0.5

    # Telomere — placed at the outer end of the contig (we don't know
    # the exact telomere repeat length from the bed; the box is a
    # symbolic indicator).
    elements.append(('telo', x, x + TELO_W, {}))
    x += TELO_W
    x += PAD

    # Y' array (outermost → innermost in display order)
    for i, (idx, s_pos, e_pos) in enumerate(yp_positions):
        yp_id = yprime_id_lookup.get((chr_end, idx), 'ID?')
        elements.append(('y_prime', x, x + YP_W, {
            'label':   f'Y-{n_yp - i}',
            'type':    yp_id,
            'idx':     idx,
            'size_bp': e_pos - s_pos,
        }))
        x += YP_W

        # Inter-Y' spacer (between THIS Y' and the next inward)
        if i < n_yp - 1:
            x += PAD
            elements.append(('inter_yp', x, x + IYP_W, {}))
            x += IYP_W
            x += PAD

    # X element
    if x_elem is not None:
        x += PAD
        elements.append(('x_element', x, x + XE_W, {
            'label':   f'X {chr_end}',
            'size_bp': x_elem[1] - x_elem[0],
        }))
        x += XE_W

    # Spacer
    if rec.get('spacer'):
        x += PAD
        elements.append(('spacer', x, x + SP_W, {
            'size_bp': rec['spacer'][1] - rec['spacer'][0],
        }))
        x += SP_W

    # Anchor
    if rec.get('anchor'):
        x += PAD
        elements.append(('anchor', x, x + ANC_W, {}))
        x += ANC_W

    x += PAD

    # Compute gap to telomere (how many bp the contig extends past the
    # outermost annotated feature). This catches assembly truncations.
    gap_to_telo = compute_gap_to_telo(rec, side)
    return elements, n_yp, gap_to_telo


def compute_gap_to_telo(rec: dict, side: str) -> int:
    """How many bp does the contig extend past the OUTERMOST feature
    on the telomere side?  Returns -1 if unknown (no contig length info).
    """
    # Find outermost (telomere-side) feature position
    candidates = []
    for k in ('anchor', 'x_var', 'x_core', 'spacer'):
        if rec.get(k):
            candidates.append(rec[k])
    for idx, s, e in rec.get('yprimes', []):
        candidates.append((s, e))
    if not candidates:
        return -1
    if side == 'L':
        outermost_pos = min(s for s, _ in candidates)
        return outermost_pos
    else:
        outermost_pos = max(e for _, e in candidates)
        return outermost_pos  # gap = clen - outermost; caller fills clen


def render_one_chr_end(ax, chr_end: str, rec: dict, yp_id_lookup: dict,
                       id_to_hex: dict, contig_lens: dict):
    """Render one chr_end's schematic onto a Matplotlib axis."""
    side = chr_end[-1]
    elements, n_yp, outermost_pos = build_elements(rec, chr_end, yp_id_lookup)

    # Compute gap to telomere with contig length
    contig_name = rec.get('contig')
    clen = contig_lens.get(contig_name, 0) if contig_name else 0
    if side == 'L':
        gap_to_telo = outermost_pos if outermost_pos >= 0 else -1
    else:
        gap_to_telo = (clen - outermost_pos) if (clen and outermost_pos >= 0) else -1

    # Mirror for R arms so telomere ends up on the right (matches
    # plot_read_structure.py convention).
    is_r = (side == 'R')
    x_content_end = elements[-1][2] if elements else 0.5
    if is_r:
        x_start = 0.5
        elements = [
            (kind, x_start + x_content_end - xr,
                   x_start + x_content_end - xl, info)
            for kind, xl, xr, info in reversed(elements)
        ]
        arrow_dir = 'right'
    else:
        arrow_dir = 'left'

    total_width = x_content_end + 1.2
    ax.set_xlim(-0.8, total_width)
    ax.set_ylim(-1.5, 1.5)
    ax.axis('off')

    # Rail
    if is_r:
        rail_x0, rail_x1 = 0.40, x_content_end
    else:
        rail_x0, rail_x1 = 0.5, total_width - 1.05
    ax.hlines(YC, rail_x0, rail_x1, colors=RAIL_COLOR, linewidth=1.5, zorder=1)

    above   = YP_H / 2 + LABEL_GAP
    below   = -(YP_H / 2 + LABEL_GAP)
    sq_below = -(IYP_H / 2 + LABEL_GAP)

    for kind, xl, xr, info in elements:
        mid = (xl + xr) / 2
        if kind == 'y_prime':
            draw_arrow(ax, xl, xr, YC, YP_H,
                       yp_color(info['type'], id_to_hex), direction=arrow_dir)
            ax.text(mid, YC + above, info['label'], ha='center', va='bottom',
                    fontsize=9, fontweight='bold')
            ax.text(mid, YC + below, info['type'], ha='center', va='top',
                    fontsize=7, color='#555555')
        elif kind == 'inter_yp':
            draw_box(ax, xl, xr, YC, IYP_H, INTER_YP_COLOR)
        elif kind == 'x_element':
            draw_box(ax, xl, xr, YC, XE_H, X_ELEMENT_COLOR)
            ax.text(mid, YC + above, info['label'], ha='center', va='bottom',
                    fontsize=6.5)
            if info.get('size_bp', 0) > 0:
                ax.text(mid, YC + sq_below, f"{info['size_bp']} bp",
                        ha='center', va='top', fontsize=6.5, color='#333333')
        elif kind == 'spacer':
            draw_box(ax, xl, xr, YC, SP_H, SPACER_COLOR)
            ax.text(mid, YC + above, 'spacer', ha='center', va='bottom',
                    fontsize=6.5)
            if info.get('size_bp', 0) > 0:
                ax.text(mid, YC + sq_below, f"{info['size_bp']} bp",
                        ha='center', va='top', fontsize=6.5, color='#333333')
        elif kind == 'anchor':
            draw_box(ax, xl, xr, YC, ANC_H, ANCHOR_COLOR)
            ax.text(mid, YC + above, 'anc', ha='center', va='bottom',
                    fontsize=6.5, color='#333333')
        elif kind == 'telo':
            draw_box(ax, xl, xr, YC, TELO_H, TELO_COLOR)
            ax.text(mid, YC + above, 'telo', ha='center', va='bottom',
                    fontsize=7.5, color='#333333')
            if gap_to_telo >= 0:
                ax.text(mid, YC + below, f'+{gap_to_telo} bp',
                        ha='center', va='top', fontsize=6.5, color='#666666')

    # chr_end label
    if is_r:
        ax.text(0.40, YC, chr_end, ha='right', va='center',
                fontsize=11, fontweight='bold')
    else:
        ax.text(total_width - 1.05, YC, chr_end, ha='left', va='center',
                fontsize=11, fontweight='bold')

    # QC flag: if no Y' and outermost feature is essentially at the contig
    # end, this might be a truncated chr_end.  Mark in red.
    if n_yp == 0 and gap_to_telo >= 0 and gap_to_telo < 500:
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('#C0392B')
            spine.set_linewidth(2)


# ─────────────────────────── Top-level rendering ──────────────────────────────

ALL_CHR_ENDS = [f'chr{n}{s}' for n in range(1, 17) for s in ('L', 'R')]


def render_per_chr_end(by_chr_end: dict, yp_id_lookup: dict, id_to_hex: dict,
                       contig_lens: dict, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    for ce in ALL_CHR_ENDS:
        if ce not in by_chr_end:
            continue
        fig, ax = plt.subplots(figsize=(12, 2.8))
        render_one_chr_end(ax, ce, by_chr_end[ce], yp_id_lookup,
                           id_to_hex, contig_lens)
        n_yp = len(by_chr_end[ce]['yprimes'])
        ax.set_title(f'{ce}  ·  {n_yp} Y\' element(s)  ·  '
                     f'contig={by_chr_end[ce].get("contig") or "?"}',
                     fontsize=9, pad=8, loc='left')
        out_path = os.path.join(out_dir, f'{ce}.png')
        fig.savefig(out_path, dpi=180, bbox_inches='tight')
        plt.close(fig)


def render_combined(by_chr_end: dict, yp_id_lookup: dict, id_to_hex: dict,
                    contig_lens: dict, out_path: str) -> None:
    """Combined overview: 16 rows × 2 cols, L on left, R on right."""
    fig, axes = plt.subplots(16, 2, figsize=(28, 32),
                              gridspec_kw={'hspace': 0.4, 'wspace': 0.15})
    for chr_num in range(1, 17):
        for col_i, side in enumerate(('L', 'R')):
            ce = f'chr{chr_num}{side}'
            ax = axes[chr_num - 1, col_i]
            if ce in by_chr_end:
                render_one_chr_end(ax, ce, by_chr_end[ce], yp_id_lookup,
                                   id_to_hex, contig_lens)
                n_yp = len(by_chr_end[ce]['yprimes'])
                ax.set_title(f'{ce}  ·  {n_yp} Y\'', fontsize=9, loc='left')
            else:
                ax.axis('off')
                ax.text(0.5, 0.5, f'{ce}: no labels in bed',
                        ha='center', va='center', fontsize=10,
                        color='#999999', transform=ax.transAxes)
    fig.suptitle('Reference structure overview', fontsize=14, fontweight='bold', y=0.995)
    fig.savefig(out_path, dpi=140, bbox_inches='tight')
    plt.close(fig)


# ─────────────────────────── CLI ──────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bed', required=True,
                    help='pretelomeric_regions_<strain>_simp.bed')
    ap.add_argument('--yprime-lib', default=None,
                    help='extracted_yprimes_<strain>.fasta (defaults to a sibling '
                         'of --bed in the same directory)')
    ap.add_argument('--reference', default=None,
                    help='reference FASTA (for gap-to-telomere calculation)')
    ap.add_argument('--output-dir', required=True,
                    help='directory to write reference_overview.png + per_chr_end/')
    args = ap.parse_args()

    if not os.path.exists(args.bed):
        sys.exit(f'ERROR: bed not found: {args.bed}')

    # Auto-detect Y' library next to the bed if not given
    yp_lib = args.yprime_lib
    if yp_lib is None:
        bed_dir = os.path.dirname(os.path.abspath(args.bed))
        candidates = [f for f in os.listdir(bed_dir)
                      if f.startswith('extracted_yprimes_') and f.endswith('.fasta')]
        if candidates:
            yp_lib = os.path.join(bed_dir, candidates[0])

    id_to_hex, yp_id_lookup = parse_yprime_lib(yp_lib) if yp_lib else ({}, {})
    contig_lens = parse_contig_lengths(args.reference) if args.reference else {}
    by_chr_end = parse_bed(args.bed)

    print(f'  bed:           {args.bed}')
    print(f'  yprime-lib:    {yp_lib or "(none — will use fallback colors)"}')
    print(f'  reference:     {args.reference or "(none — gap_to_telomere not shown)"}')
    print(f'  chr_ends in bed: {len(by_chr_end)}/32')
    print(f'  output-dir:    {args.output_dir}')

    os.makedirs(args.output_dir, exist_ok=True)
    per_dir = os.path.join(args.output_dir, 'per_chr_end')
    render_per_chr_end(by_chr_end, yp_id_lookup, id_to_hex, contig_lens, per_dir)
    combined_path = os.path.join(args.output_dir, 'reference_overview.png')
    render_combined(by_chr_end, yp_id_lookup, id_to_hex, contig_lens, combined_path)
    print(f'  Wrote per-chr_end diagrams to: {per_dir}/')
    print(f'  Wrote combined overview to:    {combined_path}')


if __name__ == '__main__':
    main()
