#!/usr/bin/env python3
"""Per-read Y'-COPY plot. Anchor ALWAYS on the left, telomere on the right (reads
flipped for consistency). Each Y' copy = one box colored by identity; the gaps
BETWEEN copies (ITS = inter-Y' telomere spacer) are shaded yellow; anchor = dark
gray, telomere side = black. Right label = copy count + composition and, when the
features TSV carries v2 path columns, the inferred origin path
("gained from: chr13L[1-2]:ID2,ID1,ID2,ID1(circ x2.0 strong) > chr4R[1-2]:ID1,ID1").

Usage: plot_yprime_copies.py <features_dir> <base_name> <output_dir> [--only-gain] [--reads <file with read_ids>]
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
        yt.append(y); ytl.append(f'{r["ce"]} {r["rid"][-8:]} [{r["status"][:8]}]')
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
