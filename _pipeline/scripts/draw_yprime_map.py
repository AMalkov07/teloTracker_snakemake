#!/usr/bin/env python3
"""Per-chromosome Y'-element end map for ONE reference, in the lab's figure
style (rows = chromosomes, boxes = Y' elements colored+ID'd by sequence-variant
cluster, dashed gaps = ITS with length label, heavy caps = telomere).

Clusters come from the pipeline's own Y' clustering, encoded in the
extracted_yprimes headers (e.g. Y_Prime_chr12R2,3,4;chr4R1,2#Long/Tandem/ID1_Gray),
so the IDs/colors match the read-track plots and are consistent across strains.

Usage: draw_yprime_map.py <simp.bed> <extracted_yprimes.fasta> <strain_label> <out_png>
"""
import re, sys
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.lines as mlines, matplotlib.patches as mpatches

simp_bed, extracted_fa, strain_label, out_png = sys.argv[1:5]
# --variant: keep the colour shade of curated libraries (ID2_Red-Light vs ID2_Red-Dark)
# as a separate group -- label "2L"/"2D"/"2N", lighter/darker fill.
VARIANT = "--variant" in sys.argv[5:]
ROMAN = ["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]
COLOR_HEX = {"Gray":"#8a8a8a","Grey":"#8a8a8a","Red":"#e34948","Green":"#1a9e4a",
             "Orange":"#eb6834","Purple":"#6a4bbf","Blue":"#2a78d6","Yellow":"#e0a500",
             "Cyan":"#00b3c0","Pink":"#e87ba4","Brown":"#8a5a44","Magenta":"#d55181"}
INK, MUTED = "#0b0b0b", "#52514e"
BOX_W, BOX_H, GAP_W, STUB_W, TEL_W = 1.0, 0.62, 0.55, 0.7, 0.9

# ---- 1. element -> (ID number, color name) from extracted_yprimes headers ----
elem = {}
for line in open(extracted_fa):
    if not line.startswith(">"): continue
    h = line[1:].strip()
    members, ann = (h.split("#", 1) + [""])[:2]
    members = re.sub(r"^Y_Prime_", "", members)
    m = re.search(r"ID(\d+)_([A-Za-z]+)(?:-([A-Za-z]+))?", ann)
    if not m: continue
    idnum, colorname, shade = m.group(1), m.group(2), (m.group(3) or "")
    if VARIANT and shade:
        idnum = idnum + {"Light": "L", "Dark": "D", "Neutral": "N"}.get(shade, shade[0])
        colorname = colorname + "-" + shade
    for grp in members.split(";"):
        gm = re.match(r"(chr\d+[LR])(.+)", grp)
        if not gm: continue
        end, nums = gm.group(1), gm.group(2)
        for k in nums.split(","):
            k = k.strip()
            if k: elem[f"{end}_Y_Prime_{k}"] = (idnum, colorname)

# ---- 2. geometry from simp.bed (Y' + ITS only) ----
arms = defaultdict(list)
for line in open(simp_bed):
    f = line.rstrip("\n").split("\t")
    if len(f) < 6: continue
    name, length = f[3], int(f[5])
    if "_Y_Prime_" not in name: continue
    mm = re.search(r"chr(\d+)([LR])", name)
    if not mm: continue
    key = (int(mm.group(1)), mm.group(2))
    kind = "ITS" if name.startswith("ITS") else "Y"
    if kind == "ITS" and length > 1000: continue     # BED artifact
    arms[key].append((int(f[1]), kind, name, length))
arms = {k: [(kd, nm, ln) for _, kd, nm, ln in sorted(v)] for k, v in arms.items()}

def _shade(hex_color, shade):
    r, g, b = (int(hex_color[i:i+2], 16) for i in (1, 3, 5))
    if shade == "Light":   r, g, b = (int(c + (255 - c) * 0.40) for c in (r, g, b))
    elif shade == "Dark":  r, g, b = (int(c * 0.65) for c in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"

def color_of(name):
    idnum, cn = elem.get(name, ("?", "Gray"))
    base, _, shade = cn.partition("-")
    hexc = COLOR_HEX.get(base, "#8a8a8a")
    return (_shade(hexc, shade) if shade else hexc), idnum

# ---- 3. draw ----
fig, ax = plt.subplots(figsize=(11.5, 7.6), dpi=200)
fig.patch.set_facecolor("#fcfcfb")
used = {}   # idnum -> colorhex (for legend)
for row, chrom in enumerate(range(1, 17)):
    y = -row
    for arm, sign in (("L", -1), ("R", +1)):
        entries = arms.get((chrom, arm), [])
        if arm == "L": entries = entries[::-1]
        x = sign * 1.55
        ax.plot([x, x + sign * STUB_W], [y, y], color=INK, lw=1.4, solid_capstyle="butt", zorder=1)
        x += sign * STUB_W
        pending = None
        for kind, name, length in entries:
            if kind == "ITS":
                pending = length; continue
            gw = GAP_W if pending is not None else 0.35
            ax.plot([x, x + sign * gw], [y, y], color=MUTED, lw=1.1, ls=(0, (2.2, 2.2)), zorder=1)
            if pending is not None:
                ax.text(x + sign * gw / 2, y + BOX_H / 2 + 0.06, str(pending),
                        ha="center", va="bottom", fontsize=5.6, color=MUTED)
                pending = None
            x += sign * gw
            ch, idnum = color_of(name)
            used[idnum] = ch
            x0 = min(x, x + sign * BOX_W)
            ax.add_patch(FancyBboxPatch((x0, y - BOX_H / 2), BOX_W, BOX_H,
                boxstyle="round,pad=0,rounding_size=0.08",
                facecolor=ch, edgecolor=INK, lw=0.5, zorder=2))
            ax.text(x0 + BOX_W / 2, y, idnum, ha="center", va="center",
                    fontsize=6.4, color="white", fontweight="bold", zorder=3)
            x += sign * BOX_W
        if pending is not None:
            ax.plot([x, x + sign * GAP_W], [y, y], color=MUTED, lw=1.1, ls=(0, (2.2, 2.2)), zorder=1)
            ax.text(x + sign * GAP_W / 2, y + BOX_H / 2 + 0.06, str(pending),
                    ha="center", va="bottom", fontsize=5.6, color=MUTED)
            x += sign * GAP_W
        ax.plot([x + sign * 0.12, x + sign * (0.12 + TEL_W)], [y, y],
                color=INK, lw=2.6, solid_capstyle="butt", zorder=1)
    ax.text(0, y, f"Chr{ROMAN[chrom-1]}", ha="center", va="center",
            fontsize=8.5, fontweight="bold", color=INK)
ax.set_xlim(-15.5, 15.5); ax.set_ylim(-15.9, 1.4); ax.axis("off")
ax.set_title(f"Y′ elements at chromosome ends — {strain_label} reference" + (" (curated variants: L/D/N = Light/Dark/Neutral shade)" if VARIANT else ""), fontsize=11, color=INK, pad=6)

def _legend_key(z):
    m = re.match(r"(\d+)([A-Z]?)", z)
    return (z == "?", int(m.group(1)) if m else 0, m.group(2) if m else "")
handles = [mpatches.Patch(facecolor=used[i], edgecolor=INK, lw=0.5, label=f"ID{i}")
           for i in sorted(used, key=_legend_key)]
handles += [mlines.Line2D([], [], color=MUTED, ls=(0, (2.2, 2.2)), lw=1.1, label="ITS (number = length, bp)"),
            mlines.Line2D([], [], color=INK, lw=2.6, label="telomere")]
ax.legend(handles=handles, loc="upper center", fontsize=6.8, frameon=False,
          ncol=min(8, len(handles)), bbox_to_anchor=(0.5, -0.01))
fig.subplots_adjust(bottom=0.10, top=0.93, left=0.02, right=0.98)
fig.savefig(out_png, facecolor=fig.get_facecolor(), bbox_inches="tight")
fig.savefig(out_png.replace(".png", ".svg"), facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.close(fig)
print("wrote", out_png)
