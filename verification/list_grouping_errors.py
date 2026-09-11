#!/usr/bin/env python3
"""List the individual reads a Y'-grouping scheme gets wrong, and what it confuses with what.

Operates on an ELEMENT-LEVEL run (array_yprime_element_lib*.sh), where `y_prime_observed_array`
names, per Y' copy, the reference element RepeatMasker matched best. Truth is positional: on a
day-0 read anchored at end E whose copy count equals E's reference array length, copy i is element
i of E. A copy is an error under scheme S when S's label for the matched element differs from S's
label for the true one.

Writes three files under <out_prefix>:
  _reads.tsv       one row per misread read: the expected and observed element arrays, which copy
                   positions are wrong, and the scheme labels on both sides
  _confusion.tsv   one row per (true element -> matched element) pair, with the scheme groups and
                   how many copies and reads it accounts for
  _summary.md      both, rendered

Usage: list_grouping_errors.py --manifest <tsv: partitions features_dir sample> <scheme> <out_prefix>
                               [--id-maps <dir of <sample>_read_id_map.tsv>]
"""
import csv, glob, json, os, re, sys
from collections import Counter, defaultdict

if '--manifest' not in sys.argv:
    sys.exit(__doc__)
man = sys.argv[sys.argv.index('--manifest') + 1]
rest = [a for i, a in enumerate(sys.argv[1:], 1)
        if a != '--manifest' and sys.argv[i - 1] != '--manifest' and not a.startswith('--')]
SCHEME, OUT = rest[0], rest[1]
ID_DIR = sys.argv[sys.argv.index('--id-maps') + 1] if '--id-maps' in sys.argv else None
os.makedirs(os.path.dirname(OUT) or '.', exist_ok=True)

jobs = [tuple(l.split()) for l in open(man) if len(l.split()) == 3]
reads_out, conf = [], defaultdict(lambda: {'copies': 0, 'reads': set(), 'samples': Counter(), 'ends': Counter()})

for part_path, feat_dir, sample in jobs:
    P = json.load(open(part_path)); ids, parts = P['elements'], P['partitions']
    if SCHEME not in parts:
        sys.exit(f'scheme {SCHEME!r} not in {part_path}; have {sorted(parts)}')
    g = parts[SCHEME]
    end_of = {i: re.match(r'^E-(chr\w+?[LR])-\d+$', i).group(1) for i in ids}
    idx_of = {i: int(i.rsplit('-', 1)[1]) for i in ids}
    exp = defaultdict(list)
    for i in sorted(ids, key=lambda x: (end_of[x], idx_of[x])): exp[end_of[i]].append(i)
    idmap = {}
    if ID_DIR:
        p = os.path.join(ID_DIR, f'{sample}_read_id_map.tsv')
        if os.path.exists(p):
            with open(p) as fh:
                next(fh)
                for line in fh:
                    a, b = line.rstrip('\n').split('\t')[:2]; idmap[a] = b

    for f in sorted(glob.glob(os.path.join(feat_dir, f'{sample}_chr*_features.tsv'))):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                ce = row['chr_end']
                obs = [x for x in row.get('y_prime_observed_array', '').split(',') if x]
                if ce not in exp or not obs or len(obs) != len(exp[ce]): continue
                if any(o not in end_of for o in obs): continue
                truth = exp[ce]
                bad = [k for k, (o, t) in enumerate(zip(obs, truth)) if g[o] != g[t]]
                if not bad: continue
                rid = row['read_id']
                reads_out.append({
                    'sample': sample, 'chr_end': ce, 'read_id': rid,
                    'ont_read_id': idmap.get(rid, ''),
                    'n_copies': len(obs), 'wrong_copies': len(bad),
                    'wrong_at_positions': ','.join(str(k + 1) for k in bad),
                    'expected_elements': ','.join(e.replace('E-', '') for e in truth),
                    'matched_elements': ','.join(e.replace('E-', '') for e in obs),
                    'expected_groups': ','.join(g[e] for e in truth),
                    'assigned_groups': ','.join(g[e] for e in obs),
                    'read_length': row.get('read_length', ''),
                    'telo_side': row.get('telo_side', ''),
                })
                for k in bad:
                    key = (truth[k].replace('E-', ''), obs[k].replace('E-', ''), g[truth[k]], g[obs[k]])
                    c = conf[key]; c['copies'] += 1; c['reads'].add((sample, rid))
                    c['samples'][sample] += 1; c['ends'][ce] += 1

cols = ['sample', 'chr_end', 'read_id', 'ont_read_id', 'n_copies', 'wrong_copies',
        'wrong_at_positions', 'expected_elements', 'matched_elements',
        'expected_groups', 'assigned_groups', 'read_length', 'telo_side']
with open(OUT + '_reads.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
    for r in sorted(reads_out, key=lambda r: (r['chr_end'], r['sample'], r['read_id'])): w.writerow(r)

rows = sorted(conf.items(), key=lambda kv: -kv[1]['copies'])
with open(OUT + '_confusion.tsv', 'w', newline='') as fh:
    fh.write('true_element\tmatched_element\ttrue_group\tassigned_group\tcopies\treads\tchr_ends\tsamples\n')
    for (te, me, tg, mg), c in rows:
        fh.write(f'{te}\t{me}\t{tg}\t{mg}\t{c["copies"]}\t{len(c["reads"])}\t'
                 f'{",".join(f"{k}:{v}" for k, v in c["ends"].most_common())}\t{len(c["samples"])}\n')

md = [f'# Reads the `{SCHEME}` scheme gets wrong', '',
      f'{len(reads_out)} misread reads over {len(jobs)} day-0 populations; '
      f'{sum(c["copies"] for _, c in rows)} wrong Y\' copies.', '',
      '## What is mislabelled as what', '',
      '| true element | matched element | true group | assigned group | copies | reads | chr ends | samples |',
      '|---|---|---|---|---|---|---|---|']
for (te, me, tg, mg), c in rows:
    md.append(f'| {te} | {me} | `{tg}` | `{mg}` | {c["copies"]} | {len(c["reads"])} | '
              f'{", ".join(f"{k} ({v})" for k, v in c["ends"].most_common())} | {len(c["samples"])} |')
md += ['', '## The reads', '',
       '| sample | end | read | copies | wrong at | expected elements | matched elements |',
       '|---|---|---|---|---|---|---|']
for r in sorted(reads_out, key=lambda r: (r['chr_end'], r['sample'])):
    rid = r['ont_read_id'] or r['read_id']
    md.append(f'| {r["sample"]} | {r["chr_end"]} | `{rid}` | {r["n_copies"]} | {r["wrong_at_positions"]} | '
              f'{r["expected_elements"]} | {r["matched_elements"]} |')
open(OUT + '_summary.md', 'w').write('\n'.join(md) + '\n')
print('\n'.join(md[:6 + len(rows) + 1]))
print(f'\nwritten {OUT}_reads.tsv, {OUT}_confusion.tsv, {OUT}_summary.md')
