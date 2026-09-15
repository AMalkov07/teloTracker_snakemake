#!/usr/bin/env python3
"""Find every Y' copy that misses its reference group at a fixed identity cutoff, for one sample.

Step 1 of the per-sample pipeline (grouping + mismatch detection). Groups come from
build_cut99_groups.py; a copy is a mismatch when the group of the element it matched differs
from the group of the element that positionally belongs there.

Only reads whose Y' copy count equals the anchor end's reference array length are scored,
because only those give each copy a positional truth.

Usage: cut99_mismatch_report.py --sample <name> --groups <groups.json>
                                --results-root <dir> --tag elemYPfix --out <tsv>
"""
import argparse, csv, glob, json, os, re
from collections import defaultdict, Counter

p = argparse.ArgumentParser()
p.add_argument('--sample', required=True)
p.add_argument('--groups', required=True)
p.add_argument('--results-root', default='/home/andrey/argon_scratch/telo_sra_runs/results')
p.add_argument('--tag', default='elemYPfix')
p.add_argument('--out', required=True)
a = p.parse_args()

G = json.load(open(a.groups))['groups']
end_of = {e: re.match(r'^E-(chr\w+?[LR])-\d+$', e).group(1) for e in G}
idx_of = {e: int(e.rsplit('-', 1)[1]) for e in G}
exp = defaultdict(list)
for e in sorted(G, key=lambda x: (end_of[x], idx_of[x])): exp[end_of[e]].append(e)

rows, scored, skipped = [], 0, Counter()
pat = f'{a.results_root}/{a.sample}__{a.tag}/_pipeline/recombination/{a.sample}_chr*_features.tsv'
for f in sorted(glob.glob(pat)):
    for r in csv.DictReader(open(f), delimiter='\t'):
        ce = r['chr_end']
        obs = [x for x in r.get('y_prime_observed_array', '').split(',') if x]
        if ce not in exp: skipped["end has no reference Y'"] += 1; continue
        if not obs: skipped["no Y' copy"] += 1; continue
        if len(obs) != len(exp[ce]): skipped['copy count != reference'] += 1; continue
        if any(o not in G for o in obs): skipped['unknown element'] += 1; continue
        scored += 1
        for i, (o, t) in enumerate(zip(obs, exp[ce])):
            if G[o] == G[t]: continue
            rows.append({'sample': a.sample, 'read_id': r['read_id'], 'chr_end': ce,
                         'copy_index': i + 1, 'n_copies': len(obs),
                         'true_element': t, 'true_group': G[t],
                         'matched_element': o, 'matched_group': G[o],
                         'donor_end': end_of[o], 'same_end': end_of[o] == ce,
                         'telo_side': r.get('telo_side', ''), 'read_length': r.get('read_length', ''),
                         'yp_start': r.get('y_prime_start', ''), 'yp_end': r.get('y_prime_end', '')})

cols = ['sample','read_id','chr_end','copy_index','n_copies','true_element','true_group',
        'matched_element','matched_group','donor_end','same_end','telo_side','read_length',
        'yp_start','yp_end']
with open(a.out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
    for r in rows: w.writerow(r)
print(f'{a.sample}: {scored} reads scored, {len(rows)} mismatched copies in '
      f'{len({r["read_id"] for r in rows})} reads  ({100*len(rows)/max(scored,1):.2f}%)')
