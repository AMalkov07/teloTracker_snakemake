#!/usr/bin/env python3
"""APOBEC motif classification with DAY-0 SUBTRACTION (for semi-repetitive
features like X-elements / Y' / spacers where paralog cross-mapping can create
recurrent non-mutation variants).

Same RAW + RECURRENCE readout as apobec_classify.py, but a recurrence-supported
survivor variant is only kept as an ACQUIRED mutation if it is essentially
absent in the matched day-0 pileup (day0 alt AF < D0_MAX with day0 depth >=
D0_MINDP). Candidates that are also present in day-0 are standing/paralog
variation and are reported separately (their count measures cross-mapping
contamination).

Usage: apobec_classify_subtract.py <baseline.fa> <survivor.pileup> <day0.pileup> <tag>
"""
import sys
from collections import Counter

MINCOUNT=5      # recurrence: min reads supporting the alt in survivor
MINAF=0.2       # recurrence: min alt allele fraction in survivor
D0_MAX=0.05     # day-0 subtraction: keep only if day-0 alt AF < this
D0_MINDP=10     # ...and day-0 depth at the position >= this

base_fa, surv_pileup, day0_pileup, tag = sys.argv[1:5]

def load_fa(path):
    d={}; name=None; seq=[]
    for line in open(path):
        if line.startswith('>'):
            if name: d[name]=''.join(seq)
            name=line[1:].split()[0]; seq=[]
        else: seq.append(line.strip())
    if name: d[name]=''.join(seq)
    return d

def parse_bases(b):
    counts={'A':0,'C':0,'G':0,'T':0}; matches=0; i=0; n=len(b)
    while i<n:
        c=b[i]
        if c=='^': i+=2; continue
        if c=='$': i+=1; continue
        if c in '+-':
            j=i+1; num=''
            while j<n and b[j].isdigit(): num+=b[j]; j+=1
            i=j+(int(num) if num else 0); continue
        if c in '.,': matches+=1; i+=1; continue
        cu=c.upper()
        if cu in counts: counts[cu]+=1
        i+=1
    return matches, counts

def apobec(rb, ab, prev, nxt):
    if rb=='C' and ab=='T' and prev=='T': return ('C>T', True)
    if rb=='G' and ab=='A' and nxt=='A': return ('C>T', True)
    if rb=='C' and ab=='G' and prev=='T': return ('C>G', True)
    if rb=='G' and ab=='C' and nxt=='A': return ('C>G', True)
    return (None, False)

# index day-0: total depth per (chrom,pos) AND alt counts per (chrom,pos,alt).
# Depth MUST be tracked at every covered position, not only positions carrying an
# alt -- a true acquired mutation sits where day-0 matches the reference, so it
# has no alt record but we still need day-0's depth there to validate it.
def index_day0(path):
    dpmap={}; altmap={}
    for line in open(path):
        f=line.rstrip('\n').split('\t')
        if len(f)<5: continue
        chrom=f[0]; pos=int(f[1]); refb=f[2].upper(); bases=f[4]
        m,c=parse_bases(bases); dp=m+sum(c.values())
        if dp==0: continue
        dpmap[(chrom,pos)]=dp
        for ab,cnt in c.items():
            if cnt and ab!=refb:
                altmap[(chrom,pos,ab)]=cnt
    return dpmap, altmap

ref=load_fa(base_fa)
d0dp_map, d0alt_map=index_day0(day0_pileup)

raw_total=raw_apo=raw_apo_ct=raw_apo_cg=raw_deam=0
spec=Counter()
cand=[]   # recurrence candidates (survivor)
for line in open(surv_pileup):
    f=line.rstrip('\n').split('\t')
    if len(f)<5: continue
    chrom=f[0]; pos=int(f[1]); refb=f[2].upper(); bases=f[4]
    if chrom not in ref: continue
    s=ref[chrom]
    prev = s[pos-2] if pos-2>=0 else 'N'
    nxt  = s[pos]   if pos<len(s) else 'N'
    matches, counts=parse_bases(bases)
    depth=matches+sum(counts.values())
    if depth==0: continue
    for ab,cnt in counts.items():
        if cnt==0 or ab==refb: continue
        raw_total+=cnt; spec[f"{refb}>{ab}"]+=cnt
        cat,is_apo=apobec(refb,ab,prev,nxt)
        if (refb=='C' and ab in 'TG') or (refb=='G' and ab in 'AC'): raw_deam+=cnt
        if is_apo:
            raw_apo+=cnt
            if cat=='C>T': raw_apo_ct+=cnt
            else:          raw_apo_cg+=cnt
        af=cnt/depth
        if cnt>=MINCOUNT and af>=MINAF:
            d0dp = d0dp_map.get((chrom,pos),0)
            d0af = (d0alt_map.get((chrom,pos,ab),0)/d0dp) if d0dp else 0.0
            cand.append((chrom,pos,refb,ab,prev,nxt,cnt,depth,af,d0af,d0dp,cat,is_apo))

# bucket candidates: acquired (day-0 covered & clean), standing (day-0 also has
# the alt = paralog/pre-existing), nodata (day-0 depth too low to validate)
acquired=[c for c in cand if (c[10] >= D0_MINDP and c[9] <  D0_MAX)]
standing=[c for c in cand if (c[10] >= D0_MINDP and c[9] >= D0_MAX)]
nodata  =[c for c in cand if (c[10] <  D0_MINDP)]

def pct(a,b): return 100*a/b if b else 0.0
tpc=allC=gpa=allG=0
for s in ref.values():
    for i,ch in enumerate(s):
        if ch=='C':
            allC+=1
            if i>0 and s[i-1]=='T': tpc+=1
        elif ch=='G':
            allG+=1
            if i+1<len(s) and s[i+1]=='A': gpa+=1
opp=100*(tpc+gpa)/(allC+allG) if (allC+allG) else 0

print(f"### {tag}")
print(f"[RAW] mismatch obs: {raw_total};  APOBEC-motif: {raw_apo} ({pct(raw_apo,raw_total):.2f}%)  "
      f"[C>T {raw_apo_ct}, C>G {raw_apo_cg}];  opportunity {opp:.2f}%;  "
      f"enrichment {(pct(raw_apo,raw_total)/opp if opp else 0):.2f}x")
napo=sum(1 for c in acquired if c[12])
print(f"[RECURRENCE] survivor candidates (>= {MINCOUNT} reads, AF>= {MINAF}): {len(cand)}")
print(f"   -> standing/paralog (day-0 also has alt, excluded): {len(standing)}")
print(f"   -> no-day-0-data (day-0 < {D0_MINDP}x, cannot validate): {len(nodata)}")
print(f"   -> ACQUIRED (day-0 covered & clean): {len(acquired)}")
print(f"        APOBEC-motif among acquired: {napo}  ({pct(napo,len(acquired)):.1f}%)")
print(f"[spectrum] {dict(spec.most_common())}")
print("  --- ACQUIRED mutations ---")
for c in acquired:
    chrom,pos,rb,ab,prev,nxt,cnt,dp,af,d0af,d0dp,cat,is_apo=c
    tag2 = f"APOBEC {cat}" if is_apo else "non-motif"
    print(f"    {chrom}:{pos}  {prev}[{rb}>{ab}]{nxt}  reads={cnt}/{dp} AF={af:.2f}  day0AF={d0af:.2f}(dp{d0dp})  {tag2}")
if standing:
    print("  --- excluded standing/paralog (present in day-0) ---")
    for c in standing[:20]:
        chrom,pos,rb,ab,prev,nxt,cnt,dp,af,d0af,d0dp,cat,is_apo=c
        print(f"    {chrom}:{pos}  {prev}[{rb}>{ab}]{nxt}  survAF={af:.2f}  day0AF={d0af:.2f}(dp{d0dp})")
    if len(standing)>20: print(f"    ... (+{len(standing)-20} more)")
print()
