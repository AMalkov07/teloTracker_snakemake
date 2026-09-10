#!/usr/bin/env python3
"""Per-read Y'-COPY plot. Anchor ALWAYS on the left, telomere on the right (reads
flipped for consistency). Each Y' copy = one box colored by identity; the gaps
BETWEEN copies (ITS = inter-Y' telomere spacer) are shaded yellow; anchor = dark
gray, telomere side = black. Right label = copy count + composition and, when the
features TSV carries v2 path columns, the inferred origin path
("gained from: chr13L[1-2]:ID2,ID1,ID2,ID1(circ x2.0 strong) > chr4R[1-2]:ID1,ID1").

Usage: plot_yprime_copies.py <features_dir> <base_name> <output_dir> [--only-gain] [--schematic]
                             [--reads <file with read_ids>] [--id-map <sample>_read_id_map.tsv]
Row labels carry the chr_end when the TSV has a chr_end column (so a combined file of reads
from several ends stays readable).
"""
import os, glob, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from itertools import groupby

feat_dir, base, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
only_gain = "--only-gain" in sys.argv[4:]
# --schematic: not to scale -- every Y' copy the same width, a fixed gap between copies with
# the measured ITS length (bp) printed above it; anchor and telomere side as fixed blocks.
schematic = "--schematic" in sys.argv[4:]
# --id-map <tsv>: translate pipeline read_ids (e.g. SRR33298449.63573) to the original
# ONT read names (UUIDs) in the row labels; TSV with columns read_id, original_name.
id_map = {}
if "--id-map" in sys.argv[4:]:
    with open(sys.argv[sys.argv.index("--id-map") + 1]) as fh:
        next(fh)
        for line in fh:
            a, b = line.rstrip("\n").split("\t")[:2]; id_map[a] = b
def show_id(rid):
    return id_map.get(rid, rid)
keep_ids = None
if "--reads" in sys.argv[4:]:
    keep_ids = {l.strip() for l in open(sys.argv[sys.argv.index("--reads") + 1]) if l.strip()}
os.makedirs(out_dir, exist_ok=True)
ITS_COLOR="#F2C744"; ANCHOR_COLOR="#555555"; TELO_COLOR="#111111"; BG="#ECECEC"
_cmap=plt.get_cmap("tab20"); _idc={}
def id_color(i):
    if i not in _idc: _idc[i]=_cmap(len(_idc)%20)
    return _idc[i]
def parse_positions(s):
    out=[]
    for part in s.split(";"):
        if ":" not in part: continue
        idn,rng=part.rsplit(":",1)
        try: a,b=rng.split("-"); out.append((idn,int(a),int(b)))
        except: pass
    return out
def rle_str(ids):
    return ",".join(f"{k}×{c}" if (c:=len(list(g)))>1 else k for k,g in groupby(ids))

for f in sorted(glob.glob(f"{feat_dir}/{base}_chr*_features.tsv")):
    end=os.path.basename(f).replace(f"{base}_","").replace("_features.tsv","")
    lines=open(f).read().splitlines()
    if len(lines)<2: continue
    h=lines[0].split("\t"); idx={c:i for i,c in enumerate(h)}
    reads=[]
    for line in lines[1:]:
        p=line.split("\t")
        if len(p)<len(h): continue
        status=p[idx["y_prime_recombination_status"]]
        if only_gain and status!="Y' Gain": continue
        if keep_ids is not None and p[idx["read_id"]] not in keep_ids: continue
        copies=sorted(parse_positions(p[idx["y_prime_positions"]]), key=lambda x:x[1])
        if not copies: continue
        rlen=int(p[idx["read_length"]])
        astart=int(float(p[idx["anchor_start"]] or 0)); aend=int(float(p[idx["anchor_end"]] or 0))
        # flip so the anchor is on the LEFT (anchor internal, telomere distal on right)
        amid=(astart+aend)/2; ymid=(copies[0][1]+copies[-1][2])/2
        flip = amid > ymid
        tf = (lambda a,b: (rlen-b, rlen-a)) if flip else (lambda a,b: (a,b))
        astart,aend = tf(astart,aend)
        copies=sorted([(idn,)+tf(a,b) for idn,a,b in copies], key=lambda x:x[1])
        reads.append(dict(rid=p[idx["read_id"]], rlen=rlen, astart=astart, aend=aend,
                          status=status, copies=copies, ce=p[idx["chr_end"]] if "chr_end" in idx else end,
                          path=p[idx["y_prime_path"]] if "y_prime_path" in idx else ""))
    if not reads: continue
    reads.sort(key=lambda r:(r["ce"], r["status"], -len(r["copies"])))
    n=len(reads); maxlen=max(r["rlen"] for r in reads)
    if schematic:
        BW, GAP, AW, TW, PRE = 1.0, 0.55, 0.9, 0.7, 0.8      # box, gap, anchor, telomere, anchor->first-copy widths
        maxcop=max(len(r["copies"]) for r in reads)
        maxx=AW+PRE+maxcop*(BW+GAP)+TW
        fig,ax=plt.subplots(figsize=(max(9, 0.9*maxx+7), max(2.5, min(26, 0.55*n+1.5))))
        yt=[]; ytl=[]
        for i,r in enumerate(reads):
            y=n-1-i; x=0.0
            ax.barh(y, AW, left=x, height=0.62, color=ANCHOR_COLOR, edgecolor="none"); x+=AW
            first=r["copies"][0]; pre_bp=first[1]-r["aend"] if r["aend"]>0 else first[1]
            ax.plot([x, x+PRE], [y, y], color="#777777", lw=1.0, ls=(0,(2,2)))
            ax.text(x+PRE/2, y+0.36, f"{pre_bp:,}", ha="center", va="bottom", fontsize=5.5, color="#555555"); x+=PRE
            for j,(idn,a,b) in enumerate(r["copies"]):
                ax.barh(y, BW, left=x, height=0.62, color=id_color(idn), edgecolor="black", linewidth=0.5)
                ax.text(x+BW/2, y, idn, ha="center", va="center", fontsize=6, color="white", fontweight="bold"); x+=BW
                if j < len(r["copies"])-1:
                    its=r["copies"][j+1][1]-b
                    ax.barh(y, GAP, left=x, height=0.62, color=ITS_COLOR, edgecolor="none")
                    ax.text(x+GAP/2, y+0.36, str(its), ha="center", va="bottom", fontsize=6, color="#333333"); x+=GAP
            tail=r["rlen"]-r["copies"][-1][2]
            ax.plot([x, x+0.25], [y, y], color="#777777", lw=1.0, ls=(0,(2,2))); x+=0.25
            ax.barh(y, TW, left=x, height=0.62, color=TELO_COLOR, edgecolor="none")
            ax.text(x+TW/2, y+0.36, f"{tail:,}", ha="center", va="bottom", fontsize=5.5, color="#555555"); x+=TW
            ncop=len(r["copies"]); ids=[c[0] for c in r["copies"]]
            ax.text(x+0.15, y+0.16, f"{ncop}× ({rle_str(ids)})", va="center", fontsize=6.5, fontweight="bold")
            if r["path"]:
                path=r["path"] if len(r["path"])<=150 else r["path"][:147]+"..."
                ax.text(x+0.15, y-0.2, f"gained from: {path}", va="center", fontsize=6.5, fontweight="bold", color="#222222")
            yt.append(y); ytl.append(f'{r["ce"]} {show_id(r["rid"]) if id_map else r["rid"][-8:]} [{r["status"][:8]}]')
        ax.set_yticks(yt); ax.set_yticklabels(ytl, fontsize=6)
        ax.set_ylim(-0.7, n-0.3); ax.set_xlim(-0.2, maxx+12); ax.set_xticks([])
        for sp in ("top","right","bottom"): ax.spines[sp].set_visible(False)
        ax.set_xlabel("schematic (not to scale): anchor → Y' copies (numbers above gaps = ITS length in bp; first number = bp from anchor to first copy; last = bp to read end) → telomere side")
        ax.set_title(f"{base}: {end} — {n} read(s); each Y' copy = one box (by identity); yellow gap = ITS with its measured length\n"
                     f"right label: copy count (composition) and the inferred origin path: donor[copies]:ids(circ x repeats support) > next donor")
        handles=[mpatches.Patch(color=ANCHOR_COLOR,label="anchor"), mpatches.Patch(color=ITS_COLOR,label="ITS (inter-Y' spacer)"),
                 mpatches.Patch(color=TELO_COLOR,label="telomere side")]
        handles+=[mpatches.Patch(color=id_color(i),label=i) for i in sorted(_idc)]
        ax.legend(handles=handles, fontsize=7, loc="lower right", frameon=False)
        fig.tight_layout()
        op=f"{out_dir}/{base}_{end}_ycopies_schematic.png"
        fig.savefig(op, dpi=130, bbox_inches="tight"); plt.close(fig)
        print(f"  {end}: {n} reads -> {os.path.basename(op)}")
        continue
    fig,ax=plt.subplots(figsize=(13, max(2.5, min(20, 0.32*n+1.5))))
    yt=[]; ytl=[]
    for i,r in enumerate(reads):
        y=n-1-i
        ax.barh(y, r["rlen"], left=0, height=0.72, color=BG, edgecolor="none")           # read
        for j in range(len(r["copies"])-1):                                              # ITS gaps
            g0=r["copies"][j][2]; g1=r["copies"][j+1][1]
            if g1>g0: ax.barh(y, g1-g0, left=g0, height=0.72, color=ITS_COLOR, edgecolor="none")
        if r["aend"]>r["astart"]:                                                        # anchor (left)
            ax.barh(y, r["aend"]-r["astart"], left=r["astart"], height=0.72, color=ANCHOR_COLOR, edgecolor="none")
        last=r["copies"][-1][2]
        if r["rlen"]>last:                                                               # telomere side (right)
            ax.barh(y, r["rlen"]-last, left=last, height=0.72, color=TELO_COLOR, edgecolor="none")
        for idn,a,b in r["copies"]:                                                      # each Y' copy
            ax.barh(y, b-a, left=a, height=0.72, color=id_color(idn), edgecolor="black", linewidth=0.5)
        ncop=len(r["copies"]); ids=[c[0] for c in r["copies"]]
        label=f"{ncop}× ({rle_str(ids)})"
        if r["path"]: label+=f"   gained from: {r['path']}"
        ax.text(r["rlen"]*1.005, y, label, va="center", fontsize=6)
        yt.append(y); ytl.append(f'{r["ce"]} {show_id(r["rid"]) if id_map else r["rid"][-8:]} [{r["status"][:8]}]')
    ax.set_yticks(yt); ax.set_yticklabels(ytl, fontsize=6)
    ax.set_ylim(-0.5, n-0.5); ax.set_xlim(0, maxlen*(1.85 if any(r["path"] for r in reads) else 1.30))
    ax.set_xlabel("Position on read (bp)     ← anchor (centromere-proximal)          (telomere-distal) →")
    ax.set_title(f"{base}: {end} — {n} read(s); each Y' copy = one box (by identity); yellow = ITS spacer\n"
                 f"right label: copy count (composition) and the inferred origin path: donor[copies]:ids(circ x repeats support) > next donor")
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    handles=[mpatches.Patch(color=ANCHOR_COLOR,label="anchor"),
             mpatches.Patch(color=ITS_COLOR,label="ITS (inter-Y' spacer)"),
             mpatches.Patch(color=TELO_COLOR,label="telomere side")]
    handles+=[mpatches.Patch(color=id_color(i),label=i) for i in sorted(_idc)]
    ax.legend(handles=handles, fontsize=7, loc="center left", bbox_to_anchor=(1.0,0.5))
    fig.tight_layout()
    op=f"{out_dir}/{base}_{end}_ycopies.png"
    fig.savefig(op, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"  {end}: {n} reads -> {os.path.basename(op)}")
