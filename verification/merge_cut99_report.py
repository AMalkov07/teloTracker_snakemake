#!/usr/bin/env python3
"""Merge the per-sample cut99 outputs into one annotated table per strain.

Combines: mismatches.tsv (which copies missed their group), flank_identities.tsv (how each
half of the Y' matches the expected vs donor element), and pair_homology.tsv (how much
homology the pair offers for a crossover). Assigns a single explanation per read.

Usage: merge_cut99_report.py --samples <s1> <s2> ... --strain <name> --out-tsv <f> --out-md <f>
"""
import argparse, csv, os, re
from collections import Counter

p = argparse.ArgumentParser()
p.add_argument('--samples', nargs='+', required=True)
p.add_argument('--strain', required=True)
p.add_argument('--out-tsv', required=True)
p.add_argument('--out-md', required=True)
p.add_argument('--reports-root', default='verification/reports')
a = p.parse_args()

def pid(s):
    m = re.match(r'([\d.]+)%', s or ''); return float(m.group(1)) if m else None

rows, scored_total = [], 0
for s in a.samples:
    D = os.path.join(a.reports_root, f'cut99_{s}')
    mm = list(csv.DictReader(open(f'{D}/mismatches.tsv'), delimiter='\t'))
    fl = {r['read_id']: r for r in csv.DictReader(open(f'{D}/flank_identities.tsv'), delimiter='\t')} \
         if os.path.exists(f'{D}/flank_identities.tsv') else {}
    hom = {}
    if os.path.exists(f'{D}/pair_homology.tsv'):
        for r in csv.DictReader(open(f'{D}/pair_homology.tsv'), delimiter='\t'):
            hom[(r['expected_element'], r['donor_element'])] = r
    for m in mm:
        rid = m['read_id']
        e, d = m['true_element'].replace('E-',''), m['matched_element'].replace('E-','')
        f = fl.get(rid, {})
        ae, ad = pid(f.get('anchorHalf_vs_expected')), pid(f.get('anchorHalf_vs_donor'))
        te, td = pid(f.get('teloHalf_vs_expected')),   pid(f.get('teloHalf_vs_donor'))
        defect = (s == '6991_day0_with_selection' and m['true_element'] == 'E-chr14L-1')
        if defect:
            ev, expl = 'reference defect', \
                "chr14L-1 is mis-assembled in this sample (5,720 bp vs 6,654) -- not recombination"
        elif None in (ae, ad, te, td):
            ev, expl = 'no junction', "no junction found; donor wins across the element or signal too weak"
        else:
            da, dt = ae - ad, td - te
            if da >= 1.5 and dt >= 1.5:  ev, expl = 'strong', "mid-Y' recombination (both halves favour the right reference)"
            elif da > 0 and dt > 0:      ev, expl = 'weak',   "mid-Y' recombination (correct direction, margins < 1.5 %)"
            else:                        ev, expl = 'FAILS',  "not supported: no half is better explained by the donor"
        h = hom.get((e, d), {})
        rows.append({
            'strain': a.strain, 'sample': s, 'read_id': rid, 'chr_end': m['chr_end'],
            'copy': f"{m['copy_index']}/{m['n_copies']}",
            'expected_element': e, 'expected_group': m['true_group'],
            'donor_element': d, 'donor_group': m['matched_group'],
            'anchorHalf_vs_expected': (f.get('anchorHalf_vs_expected') or '').split(' ')[0],
            'anchorHalf_vs_donor':    (f.get('anchorHalf_vs_donor') or '').split(' ')[0],
            'anchor_margin': f'{ae-ad:+.2f}' if None not in (ae, ad) else '',
            'teloHalf_vs_expected':   (f.get('teloHalf_vs_expected') or '').split(' ')[0],
            'teloHalf_vs_donor':      (f.get('teloHalf_vs_donor') or '').split(' ')[0],
            'telo_margin': f'{td-te:+.2f}' if None not in (te, td) else '',
            'pair_longest_homology_bp': h.get('longest_block_bp',''),
            'pair_homology_identity': h.get('longest_block_identity',''),
            'evidence': ev, 'explanation': expl, 'telo_side': m['telo_side'],
        })

order = {'strong':0,'weak':1,'FAILS':2,'no junction':3,'reference defect':4}
rows.sort(key=lambda r: (order.get(r['evidence'],9), r['sample'], r['chr_end']))
cols = list(rows[0].keys())
with open(a.out_tsv,'w',newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
    for r in rows: w.writerow(r)

cnt = Counter(r['evidence'] for r in rows)
md = [f'# {a.strain} day-0: Y\' copies that miss their reference group at the 99 % cutoff', '',
      f'Samples: {", ".join(a.samples)}', '',
      '| evidence | reads | meaning |', '|---|---|---|']
mean = {'strong': "mid-Y' recombination, both halves clearly favour the right reference",
        'weak': "mid-Y' recombination, correct direction but margins < 1.5 %",
        'FAILS': 'no half is better explained by the donor -- unexplained',
        'no junction': 'donor wins across the element, or signal too weak to split',
        'reference defect': 'known mis-assembly, not recombination'}
for k in ['strong','weak','FAILS','no junction','reference defect']:
    if cnt.get(k): md.append(f'| **{k}** | {cnt[k]} | {mean[k]} |')
md += ['', '| sample | read | end | expected | donor | anchor half exp/don | Δ | telo half exp/don | Δ | pair homology | evidence |',
       '|---|---|---|---|---|---|---|---|---|---|---|']
for r in rows:
    ah = f"{r['anchorHalf_vs_expected']}/{r['anchorHalf_vs_donor']}" if r['anchor_margin'] else '—'
    th = f"{r['teloHalf_vs_expected']}/{r['teloHalf_vs_donor']}" if r['telo_margin'] else '—'
    hm = f"{r['pair_longest_homology_bp']}bp @{r['pair_homology_identity']}" if r['pair_longest_homology_bp'] else '—'
    md.append(f"| {r['sample'].replace('6991_day0','').replace('7172_day0_with_selection','7172') or 'day0'} "
              f"| {r['read_id']} | {r['chr_end']} | {r['expected_element']} | {r['donor_element']} | {ah} "
              f"| {r['anchor_margin'] or '—'} | {th} | {r['telo_margin'] or '—'} | {hm} | {r['evidence']} |")
open(a.out_md,'w').write('\n'.join(md) + '\n')
print(f'{a.strain}: {len(rows)} rows  ' + '  '.join(f'{k}={v}' for k,v in cnt.most_common()))
print(f'written {a.out_tsv}\nwritten {a.out_md}')
