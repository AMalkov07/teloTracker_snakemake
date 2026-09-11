#!/usr/bin/env python3
"""Y'-grouping comparison: the pipeline's silhouette clustering vs the curated list.

Draws, for one strain's day-0 reference:
  A  every chromosome end's Y' array, each copy split horizontally --
     top half = silhouette group, bottom half = curated variant
  B  contingency table (silhouette cluster x curated variant), counts
  C  group-size distribution for both schemes

Usage: plot_grouping_comparison.py <strain> <ref_name> <out.png>
"""
import os, re, sys
from collections import Counter, defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
from verify_day0_reference import parse_lib, load_bed, end_sort_key

strain, ref, out_png = sys.argv[1], sys.argv[2], sys.argv[3]
SNAP = f'verification/snapshot/{ref}/pretelomeric_labels'
bed  = f'{SNAP}/pretelomeric_regions_{ref}_simp.bed'
sil  = f'{SNAP}/extracted_yprimes_{ref}.fasta'
cur  = f'verification/curated_refs/{strain}_features/repeatmasker_{strain}_all_y_primes.fasta'

sil_e,_,_ = parse_lib(sil, 'variant')           # our clusters:  ID1_Gray, ID2_Red, ...
cur_e,_,_ = parse_lib(cur, 'variant')           # curated:       ID2_Red-Light ...
cur_f,_,_ = parse_lib(cur, 'family')

ROMAN=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]
BASE={"Gray":"#8a8a8a","Grey":"#8a8a8a","Red":"#e34948","Green":"#1a9e4a","Orange":"#eb6834",
      "Purple":"#6a4bbf","Blue":"#2a78d6","Yellow":"#e0a500","Cyan":"#00b3c0","Pink":"#e87ba4","Brown":"#8a5a44"}
def shade(h,s):
    r,g,b=(int(h[i:i+2],16) for i in (1,3,5))
    if s=="Light": r,g,b=(int(c+(255-c)*0.42) for c in (r,g,b))
    elif s=="Dark": r,g,b=(int(c*0.62) for c in (r,g,b))
    return f"#{r:02x}{g:02x}{b:02x}"
def col(group):                                  # 'ID2_Red-Light' / 'ID1_Gray' -> hex
    m=re.match(r'ID\d+_([A-Za-z]+)(?:-([A-Za-z]+))?', group or '')
    if not m: return "#bbbbbb"
    return shade(BASE.get(m.group(1),"#8a8a8a"), m.group(2) or "")
def ink(hexc):                                   # readable text colour on that fill
    r,g,b=(int(hexc[i:i+2],16) for i in (1,3,5))
    return "#111111" if (0.299*r+0.587*g+0.114*b) > 150 else "white"
def short(group, fam=False):
    m=re.match(r'(ID\d+)_([A-Za-z]+)(?:-([A-Za-z]+))?', group or '')
    if not m: return group or '?'
    return m.group(1)[2:] if fam else m.group(1)[2:]+{'Light':'L','Dark':'D','Neutral':'N'}.get(m.group(3) or '','')

# --- element order per arm, from the BED ---
arms=defaultdict(list)
for f in load_bed(bed):
    if f['ftype']!='y_prime': continue
    m=re.match(r'(chr\d+[LR])_Y_Prime_(\d+)$', f['name'])
    if m: arms[m.group(1)].append(int(m.group(2)))
for k in arms: arms[k]=sorted(arms[k])

fig=plt.figure(figsize=(19,11.5)); fig.patch.set_facecolor("#fcfcfb")
gs=fig.add_gridspec(2,2, height_ratios=[2.45,1], width_ratios=[1.35,1], hspace=0.17, wspace=0.16)
axA=fig.add_subplot(gs[0,:]); axB=fig.add_subplot(gs[1,0]); axC=fig.add_subplot(gs[1,1])

BW,GAP,STUB,TEL,H = 1.0,0.30,0.75,0.8,0.70
HILITE = os.environ.get('HILITE','chr13L')
HILITE_NOTE = os.environ.get('HILITE_NOTE',
    "the curated ID7 (yellow) occurs only here, so one copy names chr13L;\nthe silhouette group 2 it sits in also covers chr8L, chr8R, chr12L and chr14L")
for row,chrom in enumerate(range(1,17)):
    y=-row
    for arm,sign in (("L",-1),("R",1)):
        ce=f"chr{chrom}{arm}"; pos=arms.get(ce,[])
        seq=pos[::-1] if arm=="L" else pos
        x=sign*1.5
        axA.plot([x,x+sign*STUB],[y,y],color="#111",lw=1.3,solid_capstyle="butt",zorder=1); x+=sign*STUB
        for p in seq:
            s=sil_e.get((ce,p),{}).get('id',''); c=cur_e.get((ce,p),{}).get('id','')
            x0=min(x,x+sign*BW)
            axA.add_patch(Rectangle((x0,y),BW,H/2,facecolor=col(s),edgecolor="#111",lw=0.4,zorder=2))
            axA.add_patch(Rectangle((x0,y-H/2),BW,H/2,facecolor=col(c),edgecolor="#111",lw=0.4,zorder=2))
            axA.text(x0+BW/2,y+H/4,short(s),ha="center",va="center",fontsize=5.8,color=ink(col(s)),fontweight="bold",zorder=3)
            axA.text(x0+BW/2,y-H/4,short(c),ha="center",va="center",fontsize=5.8,color=ink(col(c)),fontweight="bold",zorder=3)
            x+=sign*BW
            axA.plot([x,x+sign*GAP],[y,y],color="#888",lw=0.9,ls=(0,(2,2)),zorder=1); x+=sign*GAP
        axA.plot([x,x+sign*TEL],[y,y],color="#111",lw=2.4,solid_capstyle="butt",zorder=1)
        if ce == HILITE and pos:
            axA.add_patch(Rectangle((min(sign*1.5, x), y-H/2-0.12), abs(x-sign*1.5), H+0.24,
                                    fill=False, edgecolor="#d62728", lw=1.4, ls=(0,(3,2)), zorder=4))
            axA.text(sign*(abs(x)+0.35), y, HILITE_NOTE, ha="left" if sign>0 else "right",
                     va="center", fontsize=7.0, color="#d62728", style="italic", zorder=4, linespacing=1.4)
    axA.text(0,y,f"Chr{ROMAN[chrom-1]}",ha="center",va="center",fontsize=8.5,fontweight="bold")
axA.set_xlim(-14.5,15.5); axA.set_ylim(-17.0,1.5); axA.axis("off")
axA.set_title(f"Y′ grouping of the {ref} reference — upper half of each copy = silhouette clustering, "
              f"lower half = curated list\n(numbers are the group's ID; L/D/N = the curated Light/Dark/Neutral shade)",
              fontsize=11.5, pad=8)

# --- B: contingency ---
common=sorted(set(sil_e)&set(cur_e), key=lambda e:(end_sort_key(e[0]),e[1]))
S=[sil_e[e]['id'] for e in common]; C=[cur_e[e]['id'] for e in common]
srows=sorted(set(S), key=lambda z:int(re.match(r'ID(\d+)',z).group(1))); ccols=sorted(set(C))
M=[[sum(1 for a,b in zip(S,C) if a==r and b==c) for c in ccols] for r in srows]
axB.imshow([[v if v else float('nan') for v in row] for row in M],cmap="Blues",aspect="auto",vmin=0,vmax=max(max(r) for r in M))
for i,r in enumerate(M):
    for j,v in enumerate(r):
        if v: axB.text(j,i,str(v),ha="center",va="center",fontsize=8,
                       color="white" if v>max(max(x) for x in M)*0.55 else "#222",fontweight="bold")
axB.set_xticks(range(len(ccols))); axB.set_xticklabels(ccols,rotation=45,ha="right",fontsize=7.5)
axB.set_yticks(range(len(srows))); axB.set_yticklabels([f"{s}  (n={S.count(s)})" for s in srows],fontsize=7.5)
axB.set_xlabel("curated variant",fontsize=9); axB.set_ylabel("silhouette cluster",fontsize=9)
from sklearn.metrics import adjusted_rand_score
_ariv=adjusted_rand_score(C,S); _arif=adjusted_rand_score([cur_f[e]['id'] for e in common],S)
axB.set_title(f"How the {len(common)} Y′ copies map between the two schemes\n"
              f"Adjusted Rand Index {_ariv:.2f} vs curated variants, {_arif:.2f} vs curated families "
              f"(1.0 = identical grouping)",fontsize=10)

# --- C: group sizes ---
sc=Counter(S); cc=Counter(C); fc=Counter(cur_f[e]['id'] for e in common)
sets=[("silhouette",sorted(sc.values(),reverse=True),"#2a78d6"),
      ("curated families",sorted(fc.values(),reverse=True),"#1a9e4a"),
      ("curated variants",sorted(cc.values(),reverse=True),"#eb6834")]
w=0.27
for k,(lab,vals,c) in enumerate(sets):
    for i,v in enumerate(vals):
        axC.bar(i+(k-1)*w, v, width=w, color=c, edgecolor="#222", lw=0.4, label=lab if i==0 else None)
        axC.text(i+(k-1)*w, v+0.25, str(v), ha="center", fontsize=6.5)
axC.set_xlabel("group, largest first",fontsize=9); axC.set_ylabel("Y′ copies in the group",fontsize=9)
axC.set_title(f"Group sizes: {len(sc)} silhouette clusters, {len(fc)} curated families, {len(cc)} curated variants",fontsize=10)
axC.legend(fontsize=8,frameon=False); axC.set_xticks(range(max(len(v) for _,v,_ in sets)))
axC.set_xticklabels([str(i+1) for i in range(max(len(v) for _,v,_ in sets))],fontsize=7)
for sp in ("top","right"): axC.spines[sp].set_visible(False)

axA.text(0,-16.6,"each Y′ copy is drawn once, split horizontally:  upper half = silhouette group   |   "
         "lower half = curated variant   |   dashed line = ITS   |   heavy bar = telomere",
         ha="center",va="center",fontsize=8.5,color="#333")
fig.savefig(out_png,dpi=160,facecolor=fig.get_facecolor(),bbox_inches="tight")
fig.savefig(out_png.replace('.png','.svg'),facecolor=fig.get_facecolor(),bbox_inches="tight")
print("wrote",out_png)
