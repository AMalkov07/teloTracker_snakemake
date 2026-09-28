#!/usr/bin/env python3
"""Per-chromosome Y'-element end map for ONE reference, in the lab's figure style:
rows = chromosomes, left arms drawn leftwards and right arms rightwards from the label,
boxes = Y' elements coloured by sequence group, dashed lines = telomeric repeat (the
ITS between elements, labelled with its length, and the terminal telomere).

Groups come from the pipeline's own Y' clustering, encoded in the extracted_yprimes
headers (e.g. Y_Prime_chr12R2,3,4;chr4R1,2#Long/Tandem/ID1_Gray), so the IDs and colours
match the recombination outputs and are consistent across strains.

Styles:
  paper   (default) the published reference figure: blank boxes, dashed telomeres, an end
          without Y' drawn as a solid line then a dashed telomere, a "Without Y' / With Y'"
          key, and a colour key naming each group (ID, size class, number of copies).
  classic the earlier look: ID numbers in the boxes, a heavy telomere bar.

Usage:
  draw_yprime_map.py <simp.bed> <extracted_yprimes.fasta> <strain_label> <out_png>
                     [--style paper|classic] [--labels] [--no-labels] [--variant] [--no-title]

  --labels / --no-labels  print the group number in each box (default: off for paper,
                          on for classic)
  --variant  keep the colour shade of curated libraries (ID2_Red-Light vs ID2_Red-Dark)
             as a separate group: label "2L"/"2D"/"2N", lighter/darker fill
  --no-title leave the title off (for figure panels)

Writes <out_png> plus .svg and .pdf next to it.
"""
import argparse
import re
from collections import defaultdict, Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII",
         "XIV", "XV", "XVI"]
# Every colour name cluster_yprimes_paper_method.py can write (ID_COLOR_NAMES), so no group
# falls back to grey and looks like the grey group.
COLOR_HEX = {
    "Gray": "#8a8d8f", "Grey": "#8a8d8f", "Red": "#e34948", "Green": "#1a9e4a",
    "Orange": "#f07f14", "Purple": "#6a4bbf", "Blue": "#2a78d6", "Yellow": "#e0a500",
    "Cyan": "#00b3c0", "Magenta": "#c13d96", "Brown": "#8a5a44", "Pink": "#f19a9a",
    "Teal": "#2f8f8a", "Olive": "#8f9a2c", "Navy": "#26407a", "Coral": "#ff7f61",
    "Lavender": "#b9a2e0", "Maroon": "#8e1f2a", "Gold": "#c9a227", "Lime": "#9ed84a",
    "Slate": "#5f6f82",
}
SPARE = ["#b15928", "#6a3d9a", "#33a02c", "#1f78b4", "#e7298a", "#a6761d", "#1b9e77",
         "#d95f02", "#7570b3", "#66a61e"]
INK, MUTED = "#0b0b0b", "#52514e"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("simp_bed")
    ap.add_argument("extracted_fa")
    ap.add_argument("strain_label")
    ap.add_argument("out_png")
    ap.add_argument("--style", choices=["paper", "classic"], default="paper")
    ap.add_argument("--labels", dest="labels", action="store_true", default=None)
    ap.add_argument("--no-labels", dest="labels", action="store_false")
    ap.add_argument("--variant", action="store_true")
    ap.add_argument("--no-title", action="store_true")
    a = ap.parse_args()
    if a.labels is None:
        a.labels = a.style == "classic"
    return a


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------

def _shade(hex_color, shade):
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (1, 3, 5))
    if shade == "Light":
        r, g, b = (int(c + (255 - c) * 0.40) for c in (r, g, b))
    elif shade == "Dark":
        r, g, b = (int(c * 0.65) for c in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"


def load_groups(extracted_fa, variant):
    """{element name: group id}, {group id: colour hex}, {group id: size class}."""
    elem, color_name, size = {}, {}, {}
    for line in open(extracted_fa):
        if not line.startswith(">"):
            continue
        members, ann = (line[1:].strip().split("#", 1) + [""])[:2]
        members = re.sub(r"^Y_Prime_", "", members)
        m = re.search(r"ID(\d+)_([A-Za-z]+)(?:-([A-Za-z]+))?", ann)
        if not m:
            continue
        gid, cname, shade = m.group(1), m.group(2), (m.group(3) or "")
        if variant and shade:
            gid += {"Light": "L", "Dark": "D", "Neutral": "N"}.get(shade, shade[0])
            cname += "-" + shade
        color_name[gid] = cname
        size.setdefault(gid, set()).add(ann.split("/")[0] if "/" in ann else "")
        for grp in members.split(";"):
            gm = re.match(r"(chr\d+[LR])(.+)", grp)
            if not gm:
                continue
            for k in gm.group(2).split(","):
                if k.strip():
                    elem[f"{gm.group(1)}_Y_Prime_{k.strip()}"] = gid
    colors, taken = {}, set()
    spare = iter(SPARE)
    for gid in sorted(color_name, key=group_key):
        base, _, shade = color_name[gid].partition("-")
        hexc = COLOR_HEX.get(base)
        hexc = _shade(hexc, shade) if (hexc and shade) else hexc
        if hexc is None or hexc in taken:             # unknown name, or already used
            hexc = next(spare, "#999999")
        colors[gid] = hexc
        taken.add(hexc)
    sizes = {g: "/".join(sorted(s for s in v if s)).lower() for g, v in size.items()}
    return elem, colors, sizes


def load_arms(simp_bed):
    """{(chrom number, 'L'|'R'): [(kind, name, length)]} from the telomere-proximal side
    of the label outwards: ITS before element 1, element 1, ITS 1-2, element 2, ..."""
    arms = defaultdict(list)
    for line in open(simp_bed):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or "_Y_Prime_" not in f[3]:
            continue
        name, length = f[3], int(f[5])
        mm = re.search(r"chr(\d+)([LR])", name)
        if not mm:
            continue
        kind = "ITS" if name.startswith("ITS") else "Y"
        if kind == "ITS" and length > 1000:           # BED artefact
            continue
        arms[(int(mm.group(1)), mm.group(2))].append((int(f[1]), kind, name, length))
    out = {}
    for key, v in arms.items():
        ordered = [(kd, nm, ln) for _, kd, nm, ln in sorted(v)]
        out[key] = ordered[::-1] if key[1] == "L" else ordered
    return out


def group_key(z):
    m = re.match(r"(\d+)([A-Z]?)", z)
    return (z == "?", int(m.group(1)) if m else 0, m.group(2) if m else "")


# ---------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------

def draw_paper(arms, elem, colors, sizes, label, out_png, labels, title):
    BOX_W, BOX_H, GAP_W, FIRST_GAP, TEL_W, SOLID_W = 1.0, 0.62, 0.62, 0.62, 1.0, 1.0
    LABEL_HALF = 1.45
    copies = Counter(elem.values())

    def arm_width(entries):
        if not any(k == "Y" for k, _, _ in entries):
            return SOLID_W + TEL_W
        w, first = 0.0, True
        for kind, _, _ in entries:
            if kind == "ITS":
                w += FIRST_GAP if first else GAP_W
            else:
                w += (0.12 if first and w == 0 else 0) + BOX_W
                first = False
        return w + TEL_W

    left = max(arm_width(arms.get((c, "L"), [])) for c in range(1, 17))
    right = max(arm_width(arms.get((c, "R"), [])) for c in range(1, 17))
    xmin, xmax = -(LABEL_HALF + left + 0.4), LABEL_HALF + right + 0.4
    rows = 16
    key_rows = 2 + (len(colors) + 5) // 6
    fig_w = max(9.0, (xmax - xmin) * 0.52)
    fig_h = (rows + key_rows * 1.3 + 1.2) * 0.42
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=200)
    fig.patch.set_facecolor("white")

    def dashed(x0, x1, y, dash=0.13, gap=0.085):
        """Explicit dashes in data units, spread so both ends are whole dashes (matplotlib's
        own dash pattern leaves a stub wherever the length is not a whole number of periods)."""
        a, b = sorted((x0, x1))
        n = max(1, int(round((b - a + gap) / (dash + gap))))
        g = ((b - a) - n * dash) / (n - 1) if n > 1 else 0
        segs = [(a + i * (dash + g), a + i * (dash + g) + dash) for i in range(n)]
        for s, e in segs:
            ax.plot([s, e], [y, y], color=INK, lw=1.7, solid_capstyle="butt", zorder=1)

    def its_label(x, y, v):
        ax.text(x, y + 0.13, str(v), ha="center", va="bottom", fontsize=7, fontweight="bold", color=INK)

    for row, chrom in enumerate(range(1, 17)):
        y = -row
        ax.text(0, y, f"Chr{ROMAN[chrom - 1]}", ha="center", va="center", fontsize=15,
                fontweight="bold", color=INK)
        for arm, sign in (("L", -1), ("R", +1)):
            entries = arms.get((chrom, arm), [])
            x = sign * LABEL_HALF
            if not any(k == "Y" for k, _, _ in entries):
                ax.plot([x, x + sign * SOLID_W], [y, y], color=INK, lw=1.8, solid_capstyle="butt")
                dashed(x + sign * SOLID_W, x + sign * (SOLID_W + TEL_W), y)
                continue
            pending, first = None, True
            for kind, name, length in entries:
                if kind == "ITS":
                    pending = length
                    continue
                if pending is not None:
                    gw = FIRST_GAP if first else GAP_W
                    dashed(x, x + sign * gw, y)
                    its_label(x + sign * gw / 2, y, pending)
                    x += sign * gw
                    pending = None
                elif first:
                    x += sign * 0.12
                first = False
                gid = elem.get(name, "?")
                x0 = min(x, x + sign * BOX_W)
                ax.add_patch(FancyBboxPatch((x0, y - BOX_H / 2), BOX_W, BOX_H,
                                            boxstyle="round,pad=0,rounding_size=0.06",
                                            facecolor=colors.get(gid, "#cccccc"), edgecolor="none", zorder=2))
                if labels:
                    ax.text(x0 + BOX_W / 2, y, gid, ha="center", va="center", fontsize=7,
                            color="white", fontweight="bold", zorder=3)
                x += sign * BOX_W
            dashed(x, x + sign * TEL_W, y)

    # key: without / with Y'
    ky = -rows - 0.4
    kx = -4.6
    ax.plot([kx - 2.4, kx - 0.6], [ky, ky], color=INK, lw=1.8, solid_capstyle="butt")
    dashed(kx - 0.6, kx + 1.2, ky)
    ax.text(kx + 0.3, ky + 0.15, "telomere", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.text(kx - 0.6, ky - 0.75, "Without Y′", ha="center", va="center", fontsize=14, fontweight="bold")
    kx = 2.0
    dashed(kx - 0.9, kx, ky)
    ax.text(kx - 0.45, ky + 0.15, "ITS", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.add_patch(Rectangle((kx, ky - BOX_H / 2), 1.8, BOX_H, facecolor="white", edgecolor=INK, lw=1.6))
    ax.text(kx + 0.9, ky, "Y′", ha="center", va="center", fontsize=11, fontweight="bold")
    dashed(kx + 1.8, kx + 3.6, ky)
    ax.text(kx + 2.7, ky + 0.15, "telomere", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.text(kx + 1.2, ky - 0.75, "With Y′", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.text(kx + 1.2, ky - 1.25, "(colors represent Y′ sequence groups)", ha="center", va="center", fontsize=9)

    # colour key: one entry per group
    gids = sorted(colors, key=group_key)
    per_row = 6
    col_w = (xmax - xmin) / per_row
    for i, gid in enumerate(gids):
        r, c = divmod(i, per_row)
        gx, gy = xmin + 0.3 + c * col_w, ky - 2.2 - r * 0.75
        ax.add_patch(FancyBboxPatch((gx, gy - 0.22), 0.7, 0.44, boxstyle="round,pad=0,rounding_size=0.05",
                                    facecolor=colors[gid], edgecolor="none"))
        n = copies.get(gid, 0)
        ax.text(gx + 0.85, gy, f"ID{gid}  {sizes.get(gid, '')}  ×{n}", ha="left", va="center", fontsize=8.5)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ky - 2.2 - ((len(gids) - 1) // per_row) * 0.75 - 0.6, 0.9)
    ax.axis("off")
    if title:
        ax.set_title(f"Y′ elements at chromosome ends — {label} reference", fontsize=12, pad=4)
    fig.tight_layout()
    return fig


def draw_classic(arms, elem, colors, label, out_png, labels, title, variant):
    BOX_W, BOX_H, GAP_W, STUB_W, TEL_W = 1.0, 0.62, 0.55, 0.7, 0.9
    fig, ax = plt.subplots(figsize=(11.5, 7.6), dpi=200)
    fig.patch.set_facecolor("#fcfcfb")
    used = {}
    for row, chrom in enumerate(range(1, 17)):
        y = -row
        for arm, sign in (("L", -1), ("R", +1)):
            entries = arms.get((chrom, arm), [])
            x = sign * 1.55
            ax.plot([x, x + sign * STUB_W], [y, y], color=INK, lw=1.4, solid_capstyle="butt", zorder=1)
            x += sign * STUB_W
            pending = None
            for kind, name, length in entries:
                if kind == "ITS":
                    pending = length
                    continue
                gw = GAP_W if pending is not None else 0.35
                ax.plot([x, x + sign * gw], [y, y], color=MUTED, lw=1.1, ls=(0, (2.2, 2.2)), zorder=1)
                if pending is not None:
                    ax.text(x + sign * gw / 2, y + BOX_H / 2 + 0.06, str(pending),
                            ha="center", va="bottom", fontsize=5.6, color=MUTED)
                    pending = None
                x += sign * gw
                gid = elem.get(name, "?")
                ch = colors.get(gid, "#8a8a8a")
                used[gid] = ch
                x0 = min(x, x + sign * BOX_W)
                ax.add_patch(FancyBboxPatch((x0, y - BOX_H / 2), BOX_W, BOX_H,
                                            boxstyle="round,pad=0,rounding_size=0.08",
                                            facecolor=ch, edgecolor=INK, lw=0.5, zorder=2))
                if labels:
                    ax.text(x0 + BOX_W / 2, y, gid, ha="center", va="center",
                            fontsize=6.4, color="white", fontweight="bold", zorder=3)
                x += sign * BOX_W
            ax.plot([x + sign * 0.12, x + sign * (0.12 + TEL_W)], [y, y],
                    color=INK, lw=2.6, solid_capstyle="butt", zorder=1)
        ax.text(0, y, f"Chr{ROMAN[chrom - 1]}", ha="center", va="center",
                fontsize=8.5, fontweight="bold", color=INK)
    ax.set_xlim(-15.5, 15.5)
    ax.set_ylim(-15.9, 1.4)
    ax.axis("off")
    if title:
        ax.set_title(f"Y′ elements at chromosome ends — {label} reference"
                     + (" (curated variants: L/D/N = Light/Dark/Neutral shade)" if variant else ""),
                     fontsize=11, color=INK, pad=6)
    handles = [mpatches.Patch(facecolor=used[i], edgecolor=INK, lw=0.5, label=f"ID{i}")
               for i in sorted(used, key=group_key)]
    handles += [mlines.Line2D([], [], color=MUTED, ls=(0, (2.2, 2.2)), lw=1.1, label="ITS (number = length, bp)"),
                mlines.Line2D([], [], color=INK, lw=2.6, label="telomere")]
    ax.legend(handles=handles, loc="upper center", fontsize=6.8, frameon=False,
              ncol=min(8, len(handles)), bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(bottom=0.10, top=0.93, left=0.02, right=0.98)
    return fig


def main():
    a = parse_args()
    elem, colors, sizes = load_groups(a.extracted_fa, a.variant)
    arms = load_arms(a.simp_bed)
    if a.style == "paper":
        fig = draw_paper(arms, elem, colors, sizes, a.strain_label, a.out_png, a.labels, not a.no_title)
    else:
        fig = draw_classic(arms, elem, colors, a.strain_label, a.out_png, a.labels, not a.no_title, a.variant)
    stem = a.out_png[:-4] if a.out_png.endswith(".png") else a.out_png
    for ext in (".png", ".svg", ".pdf"):
        fig.savefig(stem + ext, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem + ".png (+ .svg, .pdf)")


if __name__ == "__main__":
    main()
