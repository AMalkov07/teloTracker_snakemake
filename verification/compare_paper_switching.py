#!/usr/bin/env python3
"""Compare our curated-library recombination calls with Supplementary Data 6
("mph1 template switching", 41467_2026_72032_MOESM8_ESM.xlsx).

The paper's Y/N flag is per read and means "the donor template changed within the read",
not "several Y' IDs are present": ID7,ID8,ID7,ID8,ID7 at one end is N (chr13L copied
repeatedly) while ID8,ID7,ID3 is Y (chr13L then chr14L). Our equivalent is the number of
donor blocks in y_prime_path -- consecutive path segments whose candidate-donor sets
intersect are merged, and >= 2 blocks means the template changed.

Usage: compare_paper_switching.py <xlsx> <snapshot> <out_dir>
Sample <-> (strain, PD) is resolved by matching ONT read UUIDs against the read-id maps.
"""
import os, re, sys, glob, json
from collections import Counter
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
from verify_recombination import load_features, md_table

xlsx, snap, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)

d = pd.read_excel(xlsx, sheet_name='mph1 template switching', header=0)
d.columns = ['strain','PD','read_id','chr_end','yprime_id','yprime_group','sw_score','start','end','switch']
for c in ['strain','PD','read_id','chr_end']: d[c] = d[c].ffill()
d = d[d.strain.isin([7172,7302])].copy()
d['strain'] = d.strain.astype(int); d['PD'] = d.PD.astype(int)
d['switch'] = d.groupby('read_id')['switch'].ffill()
d['ID'] = d.yprime_group.str.split('/').str[-1].str.split('_').str[0]

# read UUID -> (our sample, our read_id)
uu = set(d.read_id)
owner = {}
for p in sorted(glob.glob('verification/read_id_maps/*_read_id_map.tsv')):
    s = os.path.basename(p).replace('_read_id_map.tsv','')
    with open(p) as fh:
        next(fh)
        for line in fh:
            a,b = line.rstrip('\n').split('\t')[:2]
            if b in uu: owner[b] = (s,a)

def donor_blocks(path, self_end):
    """y_prime_path -> list of candidate-donor sets, consecutive intersecting ones merged."""
    if not isinstance(path,str) or not path: return []
    out=[]
    for seg in path.split(' > '):
        dn = seg.split(':')[0].replace('?','')
        dn = re.sub(r'\[.*','',dn)
        st = {self_end} if dn=='self' else ({'?'} if dn in ('','?') else set(dn.split('|')))
        if out and (out[-1] & st): out[-1] &= st
        else: out.append(set(st))
    return out

rows=[]
for rid, g in d.groupby('read_id'):
    g = g.sort_values('start')
    st, pdv, ce, sw = int(g.strain.iloc[0]), int(g.PD.iloc[0]), g.chr_end.iloc[0], g.switch.iloc[0]
    samp, ourid = owner.get(rid, (None,None))
    rows.append({'strain':st,'PD':pdv,'read_uuid':rid,'sample':samp,'our_read_id':ourid,
                 'paper_chr_end':ce,'paper_n_yprime':len(g),'paper_ids':','.join(g.ID),'paper_switch':sw})
P = pd.DataFrame(rows)

feat={}
for s in P['sample'].dropna().unique():
    f = load_features(snap, s)
    if not f.empty: feat[s] = f.set_index('read_id')

def col(r, name):
    f = feat.get(r['sample'])
    if f is None or r['our_read_id'] not in f.index: return None
    v = f.loc[r['our_read_id'], name]
    return v.iloc[0] if hasattr(v,'iloc') else v

for c,n in [('our_chr_end','chr_end'),('our_n_yprime','y_prime_count_on_read'),('our_ids','y_prime_observed_array'),
            ('our_status','y_prime_recombination_status'),('our_path','y_prime_path'),
            ('our_entries','y_prime_entries'),('our_source','recombination_source'),
            ('our_mechanism','recombination_mechanism')]:
    P[c] = P.apply(lambda r: col(r,n), axis=1)
P['our_blocks'] = [' > '.join('|'.join(sorted(b)) for b in donor_blocks(p,ce)) if p else ''
                   for p,ce in zip(P.our_path, P.paper_chr_end)]
P['our_switch'] = [len(donor_blocks(p,ce))>=2 if isinstance(p,str) else None
                   for p,ce in zip(P.our_path, P.paper_chr_end)]
P.to_csv(os.path.join(out,'paper_vs_ours_per_read.tsv'), sep='\t', index=False)

rep=[f'# Supplementary Data 6 vs our curated-library calls\n', f'snapshot: `{snap}`\n']
summary=[]
for (st,pdv), g in P.groupby(['strain','PD']):
    samples = ', '.join(f"{s} ({n})" for s,n in Counter(g['sample'].dropna()).most_common())
    have = g[g.our_path.notna()]
    ce_ok = int((have.paper_chr_end==have.our_chr_end).sum())
    n_ok  = int((have.paper_n_yprime==have.our_n_yprime).sum())
    v = have[have.our_switch.notna()]
    tp=int(((v.our_switch)&(v.paper_switch=='Y')).sum()); fn=int(((~v.our_switch.astype(bool))&(v.paper_switch=='Y')).sum())
    fp=int(((v.our_switch)&(v.paper_switch=='N')).sum()); tn=int(((~v.our_switch.astype(bool))&(v.paper_switch=='N')).sum())
    summary.append({'strain':st,'PD':pdv,'paper_reads':len(g),'found_in_ours':len(have),
                    'chr_end_agree':ce_ok,'yprime_count_agree':n_ok,
                    'paper_Y':int((g.paper_switch=='Y').sum()),'caught':tp,'missed':fn,
                    'extra_switches':fp,'agreement_pct':round(100*(tp+tn)/max(len(v),1),1),
                    'our_samples':samples})
    rep += [f'## strain {st}, PD {pdv}  ->  {samples}',
            f'* {len(have)} of {len(g)} paper reads found in our output; chr_end agrees {ce_ok}/{len(have)}; '
            f"Y' copy count agrees {n_ok}/{len(have)}",
            f'* template switching: paper flags {int((g.paper_switch=="Y").sum())}; we catch {tp}, miss {fn}, '
            f'and flag {fp} reads the paper calls N ({round(100*(tp+tn)/max(len(v),1),1)} % agreement)\n',
            "### the paper's switching reads, as our pipeline called them",
            md_table(g[g.paper_switch=='Y'][['paper_chr_end','paper_ids','our_ids','our_blocks','our_source','our_mechanism']], 40),
            '### reads we call a switch and the paper does not',
            md_table(v[(v.our_switch)&(v.paper_switch=='N')][['paper_chr_end','paper_ids','our_ids','our_blocks']], 25), '']
S=pd.DataFrame(summary); S.to_csv(os.path.join(out,'paper_vs_ours_summary.tsv'), sep='\t', index=False)
open(os.path.join(out,'paper_comparison.md'),'w').write('\n'.join(['# Summary\n', md_table(S.drop(columns=['our_samples'])), '']+rep))
print(S.drop(columns=['our_samples']).to_string(index=False))
print('\nreports in', out)
