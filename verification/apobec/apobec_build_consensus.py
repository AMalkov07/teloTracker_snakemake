#!/usr/bin/env python3
"""Build a per-anchor day-0 consensus (the ancestral baseline) from an mpileup.

Majority base per position (requiring MIN_DEPTH); where the majority differs
from the reference anchor we record a strain-specific difference. Output is a
consensus fasta (same coords/length as the reference anchors) plus a TSV of the
differences vs the reference anchors (answers "are the day-0 anchors different
from our reference anchors?").

Usage: apobec_build_consensus.py <ref_anchors.fa> <day0.pileup> <out_consensus.fa> <out_diffs.tsv>
"""
import sys

MIN_DEPTH = 10          # min coverage to call a consensus base
MAJ_FRAC  = 0.5         # majority base must exceed this fraction of depth

ref_fa, pileup, out_fa, out_diffs = sys.argv[1:5]

def load_fa(path):
    d={}; name=None; seq=[]
    for line in open(path):
        if line.startswith('>'):
            if name: d[name]=list(''.join(seq))
            name=line[1:].split()[0]; seq=[]
        else: seq.append(line.strip())
    if name: d[name]=list(''.join(seq))
    return d

def parse_bases(b):
    """Return (matches, {A,C,G,T counts}) from an mpileup read-base string."""
    counts={'A':0,'C':0,'G':0,'T':0}; matches=0; i=0; n=len(b)
    while i<n:
        c=b[i]
        if c=='^': i+=2; continue          # ^ then mapping-qual char
        if c=='$': i+=1; continue
        if c in '+-':                       # indel: [+-]<len><seq>
            j=i+1; num=''
            while j<n and b[j].isdigit(): num+=b[j]; j+=1
            i=j+(int(num) if num else 0); continue
        if c in '.,': matches+=1; i+=1; continue
        cu=c.upper()
        if cu in counts: counts[cu]+=1
        i+=1
    return matches, counts

ref=load_fa(ref_fa)
cons={k:v[:] for k,v in ref.items()}
covered=0; diffs=[]
for line in open(pileup):
    f=line.rstrip('\n').split('\t')
    if len(f)<5: continue
    chrom=f[0]; pos=int(f[1]); refb=f[2].upper(); bases=f[4]
    if chrom not in cons: continue
    matches, counts=parse_bases(bases)
    total=matches+sum(counts.values())
    if total<MIN_DEPTH: continue
    covered+=1
    allc=dict(counts); allc[refb]=allc.get(refb,0)+matches
    maj=max(allc, key=allc.get)
    if maj!=refb and allc[maj] > total*MAJ_FRAC:
        cons[chrom][pos-1]=maj
        diffs.append((chrom,pos,refb,maj,allc[maj],total))

with open(out_fa,'w') as o:
    for k in ref:
        o.write('>'+k+'\n'); s=''.join(cons[k])
        for x in range(0,len(s),80): o.write(s[x:x+80]+'\n')

with open(out_diffs,'w') as o:
    o.write("anchor\tpos\tref\tconsensus\talt_depth\ttotal_depth\n")
    for d in diffs: o.write('\t'.join(map(str,d))+'\n')

sys.stderr.write(f"[consensus] {covered} positions covered >= {MIN_DEPTH}x; "
                 f"{len(diffs)} substitution differences vs reference anchors\n")
print(len(diffs))
