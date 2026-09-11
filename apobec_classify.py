#!/usr/bin/env python3
"""Classify survivor mismatches (vs the day-0 baseline) for the APOBEC motif.

Reports BOTH definitions of a "mutation", side by side:
  RAW           - every read-vs-baseline mismatch observation (ONT-error-laden;
                  meaningful only as A3A-vs-EV differential)
  RECURRENCE    - a baseline position is a mutation only if the alt base recurs
                  across reads (count>=MINCOUNT and AF>=MINAF) -> clonal, ONT
                  error largely removed

APOBEC motif on the reference-forward strand (strand-symmetric, both captured
because mpileup reports reverse-read bases in forward orientation):
  C>T at TpC   (ref C, alt T, 5' base = T)   and mirror  G>A at GpA (ref G, alt A, 3' base = A)
  C>G at TpC   (ref C, alt G, 5' base = T)   and mirror  G>C at GpA (ref G, alt C, 3' base = A)

Usage: apobec_classify.py <baseline.fa> <survivor.pileup> <tag>
"""
import sys
from collections import Counter

MINCOUNT=5      # recurrence: min reads supporting the alt
MINAF=0.2       # recurrence: min alt allele fraction

base_fa, pileup, tag = sys.argv[1:4]

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
    """(category, is_apobec): category in {'C>T','C>G',None}."""
    if rb=='C' and ab=='T' and prev=='T': return ('C>T', True)
    if rb=='G' and ab=='A' and nxt=='A': return ('C>T', True)
    if rb=='C' and ab=='G' and prev=='T': return ('C>G', True)
    if rb=='G' and ab=='C' and nxt=='A': return ('C>G', True)
    return (None, False)

ref=load_fa(base_fa)

raw_total=raw_apo=raw_apo_ct=raw_apo_cg=raw_deam=0
rec_total=rec_apo=rec_apo_ct=rec_apo_cg=0
spec=Counter()

for line in open(pileup):
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
        raw_total+=cnt
        spec[f"{refb}>{ab}"]+=cnt
        cat,is_apo=apobec(refb,ab,prev,nxt)
        if (refb=='C' and ab in 'TG') or (refb=='G' and ab in 'AC'):
            raw_deam+=cnt
        if is_apo:
            raw_apo+=cnt
            if cat=='C>T': raw_apo_ct+=cnt
            else:          raw_apo_cg+=cnt
        af=cnt/depth
        if cnt>=MINCOUNT and af>=MINAF:
            rec_total+=1
            if is_apo:
                rec_apo+=1
                if cat=='C>T': rec_apo_ct+=1
                else:          rec_apo_cg+=1

# opportunity (folded strands): TpC + GpA sites over all C+G
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

def pct(a,b): return 100*a/b if b else 0.0
print(f"### {tag}")
print(f"[RAW] mismatch observations: {raw_total}")
print(f"      APOBEC-motif: {raw_apo}  ({pct(raw_apo,raw_total):.2f}% of all mismatches)")
print(f"        C>T@TpC: {raw_apo_ct}   C>G@TpC: {raw_apo_cg}")
print(f"      APOBEC as % of deamination-type (C>T/C>G at C:G): {pct(raw_apo,raw_deam):.2f}%  (n={raw_deam})")
print(f"      opportunity (TpC+GpA / all C+G in baseline): {opp:.2f}%")
enr = (pct(raw_apo,raw_total)/opp) if opp else 0
print(f"      RAW enrichment (apobec%% of mism / opportunity%%): {enr:.2f}x")
print(f"[RECURRENCE] mutations (count>={MINCOUNT}, AF>={MINAF}): {rec_total}")
print(f"      APOBEC-motif: {rec_apo}  ({pct(rec_apo,rec_total):.2f}%)   [C>T {rec_apo_ct}, C>G {rec_apo_cg}]")
print(f"[spectrum] {dict(spec.most_common())}")
print()
