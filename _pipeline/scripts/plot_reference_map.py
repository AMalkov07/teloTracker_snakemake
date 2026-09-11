#!/usr/bin/env python3
"""Draw a structural map of a labeled reference: one row per chr end, oriented
anchor-LEFT / telomere-RIGHT, showing anchor, spacer, X-element, each Y' copy
(colored Y'-long vs Y'-short) and the ITS spacers between them. Y' copies are
outlined so they stay countable.

Usage: plot_reference_map.py <simp.bed> <strain_label> <output_png> [--yprime-zoom]
  --yprime-zoom : compress the internal spacer so the Y' array is large/readable
"""
import re, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

bed, label, out_png = sys.argv[1], sys.argv[2], sys.argv[3]
ZOOM = "--yprime-zoom" in sys.argv[4:]
YS_THRESH = 6000  # Y' length >= this = long, else short

C = dict(anchor="#555555", spacer="#C9A66B", x="#3E9E45",
         yl="#1f6fb0", ys="#A8CCE8", its="#F2C744", telo="#111111", bg="#EDEDED")

end_re = re.compile(r'(chr\d+[LR])')
feats = defaultdict(list)   # end -> list of (kind, start, end, length)
for line in open(bed):
    p = line.rstrip("\n").split("\t")
    if len(p) < 4: continue
    chrom, s, e, name = p[0], int(p[1]), int(p[2]), p[3]
    m = end_re.search(name)
    if not m: continue
    end = m.group(1)
    if name.startswith("ITS_"):           kind = "its"
    elif "space_between_anchor" in name:  kind = "spacer"   # must precede _anchor test
    elif name.endswith("_anchor"):        kind = "anchor"
    elif "x_core_element" in name or "x_variable_element" in name: kind = "x"
    elif "_Y_Prime_" in name:             kind = "yl" if (e-s) >= YS_THRESH else "ys"
    else: continue
    feats[end].append((kind, s, e))

def chrnum(end):
    n = int(re.search(r'\d+', end).group()); return (n, 0 if end.endswith("L") else 1)
ends = sorted(feats, key=chrnum)

# per-end: orient anchor-left, local coords (optionally zoom = compress spacer)
rows = {}
for end in ends:
    fs = feats[end]
    anc = [f for f in fs if f[0] == "anchor"]
    if not anc: continue
    a_s, a_e = anc[0][1], anc[0][2]
    lo = min(f[1] for f in fs); hi = max(f[2] for f in fs)
    flip = (a_s + a_e)/2 > (lo + hi)/2
    def loc(s, e):
        return (hi - e, hi - s) if flip else (s - lo, e - lo)
    segs = sorted([(k,)+loc(s, e) for k, s, e in fs], key=lambda x: x[1])
    if ZOOM:   # compress the spacer to a fixed 1500 bp so Y'/X stay large
        out = []; shift = 0; prev_end = 0
        for k, s, e in segs:
            s += shift; e += shift
            if k == "spacer" and (e - s) > 1500:
                new = 1500; shift += new - (e - s); e = s + new
            out.append((k, s, e))
        segs = out
    rows[end] = segs

n = len(rows); maxx = max(e for segs in rows.values() for _, _, e in segs)
fig, ax = plt.subplots(figsize=(16, max(4, 0.34*n + 1.5)))
yt = []; ytl = []
for i, end in enumerate([e for e in ends if e in rows]):
    y = n - 1 - i; segs = rows[end]
    ax.barh(y, maxx, left=0, height=0.7, color=C["bg"], edgecolor="none")
    nyl = sum(1 for k,_,_ in segs if k == "yl"); nys = sum(1 for k,_,_ in segs if k == "ys")
    for k, s, e in segs:
        if k == "its":
            ax.barh(y, e-s, left=s, height=0.7, color=C["its"], edgecolor="none")
        elif k in ("yl", "ys"):
            ax.barh(y, e-s, left=s, height=0.7, color=C[k], edgecolor="black", linewidth=0.6)
        else:
            ax.barh(y, e-s, left=s, height=0.7, color=C[k], edgecolor="none")
    # telomere cap at distal (right) end
    distal = max(e for _, _, e in segs)
    ax.barh(y, maxx*0.006, left=distal, height=0.7, color=C["telo"], edgecolor="none")
    ntot = nyl + nys
    yt.append(y); ytl.append(f"{end}  ({ntot} Y′: {nyl}L,{nys}S)")
ax.set_yticks(yt); ax.set_yticklabels(ytl, fontsize=7)
ax.set_ylim(-0.5, n-0.5); ax.set_xlim(0, maxx*1.02)
ax.set_xlabel("bp     ← anchor (centromere-proximal)          (telomere-distal) →"
              + ("   [spacer compressed]" if ZOOM else "   [to scale]"))
ax.set_title(f"{label} — reference subtelomere map ({n} chr ends)", fontsize=13)
ax.grid(axis="x", linestyle=":", alpha=0.3)
handles = [mpatches.Patch(color=C["anchor"], label="anchor"),
           mpatches.Patch(color=C["spacer"], label="spacer"),
           mpatches.Patch(color=C["x"], label="X element"),
           mpatches.Patch(color=C["yl"], label="Y′ long (≥6 kb)"),
           mpatches.Patch(color=C["ys"], label="Y′ short (<6 kb)"),
           mpatches.Patch(color=C["its"], label="ITS spacer"),
           mpatches.Patch(color=C["telo"], label="telomere")]
ax.legend(handles=handles, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
fig.tight_layout()
fig.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close(fig)
print("wrote", out_png)
