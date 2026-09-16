#!/usr/bin/env python3
"""Run our (ID, ITS) path parser on the PAPER's own per-copy Y' calls and compare the
resulting donor blocks with its template-switching flag. Needs only the Excel plus a
day-0 BED + curated library, so it works for strains whose reads we have not downloaded.

Read orientation is not in the Excel, so it is inferred: of the two orientations, take the
one whose leading copies match more of that end's reference array (ties -> more ITS
agreement). Pass --truth <features-snapshot> to score the heuristic against the real
telo_side where our own output exists.

Usage: compare_paper_paths.py <xlsx> <sheet> <strain> <ref_name> <out_dir> [--truth <snapshot>]
"""
import os, sys, re
from collections import Counter
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
import yprime_path as yp
from verify_day0_reference import parse_lib
from verify_recombination import load_features, md_table

xlsx, sheet, strain, ref, out = sys.argv[1:6]
truth = sys.argv[sys.argv.index('--truth')+1] if '--truth' in sys.argv else None
FAMILY = '--family' in sys.argv     # collapse ID2_Red-Light/ID2_Red-Dark -> ID2 on both sides
os.makedirs(out, exist_ok=True)
strain = int(strain)
TELO_FIRST_MAX = 2000   # bp: first Y' closer than this to the read start => telomere-first listing

bed = f'verification/snapshot/{ref}/pretelomeric_labels/pretelomeric_regions_{ref}_simp.bed'
elem,_,_ = parse_lib(f'verification/curated_refs/{strain}_features/repeatmasker_{strain}_all_y_primes.fasta','variant')
LOC = {k: v['id'] for k,v in elem.items()}     # our curated labels, position -> variant

d = pd.read_excel(xlsx, sheet_name=sheet, header=0)
d.columns = ['strain','PD','read_id','chr_end','yprime_id','yprime_group','sw_score','start','end','switch']
d = d[d.strain.astype(str).str.strip()!='-----'].copy()
for c in ['strain','PD','read_id','chr_end']: d[c]=d[c].ffill()
d = d[d.read_id.astype(str).str.contains('-')]
d['strain']=d.strain.astype(int); d['PD']=d.PD.astype(int)
d = d[d.strain==strain].copy()
d['switch']=d.groupby('read_id')['switch'].ffill()
# the paper's group label -> our variant vocabulary (identical strings in the curated files)
d['ID']=d.yprime_group.str.split('/').str[-1]
if FAMILY: d['ID']=d['ID'].str.split('_').str[0]

# Relabel the reference in the PAPER's vocabulary. The Excel states, for every library entry
# it used (Y'_id, e.g. "Y_Prime_chr13L4"), the group it assigns that entry to. Our copy of the
# curated file sometimes differs -- chr13L's long copies are ID8_Brown to the paper and
# ID2_Red-Light to us -- and a label our reference does not contain makes every copy of it
# unexplainable, which shows up as a spurious template switch. Taking the entry -> group map
# from the Excel removes that whole class of artefact.
def entry_locations(entry):
    body = re.sub(r'^Y_Prime_', '', str(entry))
    body = re.sub(r'(chr\d+[LR])_', r'\1', body)        # the sheet writes chr12R_3,4,5
    out=[]
    for grp in body.split(';'):
        m=re.match(r'(chr\d+[LR])([\d,]+)', grp.strip())
        if m:
            for q in m.group(2).split(','):
                if q: out.append((m.group(1), int(q)))
    return out
relabelled=0
for entry, grp in d[['yprime_id','ID']].drop_duplicates().itertuples(index=False):
    for loc in entry_locations(entry):
        if loc in LOC and LOC[loc]!=grp: relabelled+=1
        LOC[loc]=grp
if FAMILY: LOC={k: v.split('_')[0] for k,v in LOC.items()}
REF = yp.build_reference_tokens(bed, LOC)
print(f"reference relabelled to the paper's vocabulary: {relabelled} positions changed")

def donor_blocks(segs, self_end):
    out=[]
    for sg in segs:
        s={self_end} if sg['tag']=='self' else ({'?'} if sg['tag']=='unk' else set(sg['donor'].split('|')))
        if out and (out[-1] & s): out[-1] &= s
        else: out.append(set(s))
    return out

truthmap={}
if truth:
    for s in sorted(os.listdir(truth)):
        f = load_features(truth, s)
        if not f.empty and 'telo_side' in f.columns:
            truthmap.update(dict(zip(f.read_id, f.telo_side)))

rows=[]
for rid,g in d.groupby('read_id'):
    g=g.sort_values('start'); ids=list(g.ID); ce=g.chr_end.iloc[0]
    gaps=[int(g.start.iloc[i+1])-int(g.end.iloc[i]) for i in range(len(g)-1)]
    rarr=[t[0] for t in REF.get(ce,[])]
    # Orientation. The Excel lists copies in read coordinates. When the telomere is at the
    # start of the read the array begins almost immediately (first copy within ~265 bp in
    # every read we can check); otherwise the anchor, spacer and X element come first and the
    # array starts >= 4.5 kb in. The two are cleanly separated, so a 2 kb threshold recovers
    # the orientation without needing our own output.
    orients = ('rev',) if int(g.start.iloc[0]) < TELO_FIRST_MAX else ('fwd',)
    best=None
    for orient in orients:
        I,G=(ids,gaps) if orient=='fwd' else (ids[::-1],gaps[::-1])
        tok=[(I[i], G[i] if i<len(G) else None) for i in range(len(I))]
        k=0
        while k<len(tok) and k<len(rarr) and tok[k][0]==rarr[k]: k+=1
        p=yp.parse_path(tok[k:], REF, ce)
        score=(k, p['its_verified'], -p['n_segments'])
        if best is None or score>best[0]: best=(score,orient,tok,k,p)
    (k,itsv,_),orient,tok,kept,p = best
    b=donor_blocks(p['segments'], ce)
    rows.append({'read_id':rid,'PD':int(g.PD.iloc[0]),'chr_end':ce,'paper_switch':g.switch.iloc[0],
                 'n_copies':len(ids),'orient_used':orient,'kept_native':kept,
                 'ids':','.join(t[0].split('_')[0] for t in tok),
                 'variants':','.join(t[0] for t in tok),
                 'path':yp.format_path(p),'blocks':' > '.join('|'.join(sorted(x)) for x in b),
                 'n_blocks':len(b),'our_switch':len(b)>=2})
R=pd.DataFrame(rows); tag='_family' if FAMILY else ''
R.to_csv(os.path.join(out,f'{strain}_paper_paths{tag}.tsv'),sep='\t',index=False)

rep=[f'# Our path parser applied to the paper\'s own Y\' calls — strain {strain}\n',
     f'reference: `{ref}`, curated {strain} library; orientation inferred from the data.\n']
if truthmap:
    # score the orientation heuristic where our own telo_side exists
    import json
    hits=json.load(open('/tmp/uuid_hits_all.json')) if os.path.exists('/tmp/uuid_hits_all.json') else {}
    u2our={u:o for fd in hits.values() for u,o in fd.items()}
    ok=tot=0
    for _,r in R.iterrows():
        oid=u2our.get(r.read_id); ts=truthmap.get(oid)
        if ts is None: continue
        tot+=1; ok += int((ts=='beginning') == (r.orient_used=='rev'))
    if tot: rep.append(f'**Orientation heuristic: {ok}/{tot} reads ({100*ok/tot:.1f} %) match the true telo_side.**\n')
    print(f'orientation heuristic: {ok}/{tot} correct' if tot else 'no truth overlap')

summ=[]
for pdv,g in R.groupby('PD'):
    v=g[g.paper_switch.isin(['Y','N'])]
    tp=int(((v.our_switch)&(v.paper_switch=='Y')).sum()); fn=int(((~v.our_switch)&(v.paper_switch=='Y')).sum())
    fp=int(((v.our_switch)&(v.paper_switch=='N')).sum()); tn=int(((~v.our_switch)&(v.paper_switch=='N')).sum())
    summ.append({'PD':pdv,'reads':len(g),'paper_Y':tp+fn,'caught':tp,'missed':fn,'extra':fp,
                 'agreement_pct':round(100*(tp+tn)/max(len(v),1),1)})
    rep += [f'## PD {pdv} — {len(g)} reads',
            f'* paper flags {tp+fn} switching reads; we catch {tp}, miss {fn}, and flag {fp} it calls N '
            f'({round(100*(tp+tn)/max(len(v),1),1)} % agreement)\n',
            "### the paper's switching reads as our parser reads them",
            md_table(g[g.paper_switch=='Y'][['chr_end','ids','blocks','path']],40),
            '### reads we call a switch and the paper does not',
            md_table(g[(g.our_switch)&(g.paper_switch=='N')][['chr_end','ids','blocks']],25), '']
S=pd.DataFrame(summ); S.to_csv(os.path.join(out,f'{strain}_paper_paths{tag}_summary.tsv'),sep='\t',index=False)
open(os.path.join(out,f'{strain}_paper_paths{tag}.md'),'w').write('\n'.join(['# Summary\n',md_table(S),'']+rep))
print(S.to_string(index=False)); print('\n->',out)
