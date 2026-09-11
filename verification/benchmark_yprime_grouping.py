#!/usr/bin/env python3
"""Score every Y'-grouping scheme on day-0 reads, which should carry no recombination.

Input is a run made with the ELEMENT-LEVEL library (array_yprime_element_lib.sh), where
`y_prime_observed_array` names, per Y' copy, the single reference element RepeatMasker
matched best. That assignment is fixed, so each scheme is scored as a relabelling of it and
the schemes differ only in how they group, never in how they matched.

Truth is positional and scheme-independent: on a day-0 read anchored at end E whose Y' copy
count equals E's reference array length, copy i is reference element i of E. Reads whose
count differs are set aside (real sub-clonal recombination, truncation, a missed copy)
because their copy-to-element correspondence is unknown.

  copy error    a copy's scheme label differs from the label of the element that truly sits
                at that position
  FALSE CALL    any copy in the read is wrong, so the read's array string does not match its
                own anchor's reference array -- this is the read the pipeline would report as
                recombination when nothing happened
  foreign donor a false call where the wrong label belongs to a group carried ONLY by other
                chromosome ends, so the read looks like it gained Y' from somewhere else

Usage: benchmark_yprime_grouping.py --manifest <tsv: partitions_json features_dir sample> <out_dir>
       benchmark_yprime_grouping.py <partitions.json> <features_dir> <sample>
"""
import json, os, re, sys, glob, csv
from collections import Counter, defaultdict

ORDER = ['element', 'condensed', 'cut99', 'cut97', 'curated_variant', 'silhouette', 'curated_family']

def load(part_path):
    P = json.load(open(part_path))
    ids, parts = P['elements'], P['partitions']
    end_of = {i: re.match(r'^E-(chr\w+?[LR])-\d+$', i).group(1) for i in ids}
    idx_of = {i: int(i.rsplit('-', 1)[1]) for i in ids}
    exp = defaultdict(list)
    for i in sorted(ids, key=lambda x: (end_of[x], idx_of[x])): exp[end_of[i]].append(i)
    return ids, parts, end_of, dict(exp)

def rows(feat_dir, sample):
    for f in sorted(glob.glob(os.path.join(feat_dir, f'{sample}_chr*_features.tsv'))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter='\t'): yield r

def score(part_path, feat_dir, sample):
    ids, parts, end_of, exp = load(part_path)
    known = set(ids)
    schemes = [s for s in ORDER if s in parts] + [s for s in parts if s not in ORDER]
    # for each scheme, which chromosome ends carry each group
    ends_of_group = {s: defaultdict(set) for s in schemes}
    for s in schemes:
        for i in ids: ends_of_group[s][parts[s][i]].add(end_of[i])

    st = {s: Counter() for s in schemes}
    confusion = {s: Counter() for s in schemes}
    kept = 0; setaside = Counter()

    for row in rows(feat_dir, sample):
        ce = row['chr_end']
        obs = [x for x in row.get('y_prime_observed_array', '').split(',') if x]
        if ce not in exp: setaside["end has no reference Y' array"] += 1; continue
        if not obs:       setaside["no Y' copy detected"] += 1; continue
        truth = exp[ce]
        if len(obs) != len(truth):
            setaside[f'copy count {len(obs)} vs reference {len(truth)}'] += 1; continue
        if any(o not in known for o in obs): setaside['matched an unknown element'] += 1; continue
        kept += 1
        for s in schemes:
            p = parts[s]; bad = 0; foreign = False
            for o, t in zip(obs, truth):
                st[s]['copies'] += 1
                if p[o] == p[t]: continue
                bad += 1; st[s]['copy_err'] += 1
                if ce not in ends_of_group[s][p[o]]:
                    foreign = True; st[s]['foreign_copies'] += 1
                    confusion[s][(ce, '|'.join(sorted(ends_of_group[s][p[o]])))] += 1
            st[s]['reads'] += 1
            if bad:
                st[s]['false_call'] += 1
                if foreign: st[s]['foreign_reads'] += 1
    for s in schemes: st[s]['ngroups'] = len(set(parts[s].values()))
    return schemes, st, confusion, kept, setaside

HDR = ('| scheme | groups | reads | copies | copy errors | copy err % | '
       'FALSE CALLS | false call % | of those, foreign-donor | foreign % |')
SEP = '|---|---|---|---|---|---|---|---|---|---|'
def table(schemes, st):
    out = [HDR, SEP]
    for s in schemes:
        d = st[s]
        if not d['copies']: continue
        out.append(f"| {s} | {d['ngroups']} | {d['reads']} | {d['copies']} | {d['copy_err']} | "
                   f"{100*d['copy_err']/d['copies']:.2f} | {d['false_call']} | "
                   f"{100*d['false_call']/max(d['reads'],1):.2f} | {d['foreign_reads']} | "
                   f"{100*d['foreign_reads']/max(d['reads'],1):.2f} |")
    return out

if __name__ == '__main__':
    a = sys.argv[1:]
    jobs, out_dir = [], None
    if a[0] == '--manifest':
        out_dir = a[2]
        jobs = [tuple(l.split()) for l in open(a[1]) if len(l.split()) == 3]
    else:
        jobs = [(a[0], a[1], a[2])]

    md = ['# Y\' grouping schemes scored on 6991 day-0 reads', '',
          'Day-0 populations should carry no recombination, so every mismatch below is the',
          'grouping getting it wrong. A FALSE CALL is a read whose Y\' array string does not',
          'match its own anchor\'s reference array -- what the pipeline would report as',
          'recombination on a read where nothing happened.', '']
    pooled = defaultdict(Counter); pooled_conf = defaultdict(Counter); schemes = None; tot_kept = 0
    for part_path, feat_dir, sample in jobs:
        schemes, st, conf, kept, setaside = score(part_path, feat_dir, sample)
        tot_kept += kept
        for s in schemes:
            pooled[s].update({k: v for k, v in st[s].items() if k != 'ngroups'})
            pooled[s]['ngroups'] = st[s]['ngroups']
            pooled_conf[s].update(conf[s])
        blk = [f'## {sample}', f'{kept} non-recombinant reads scored; '
               f'{sum(setaside.values())} set aside ({", ".join(f"{k}: {v}" for k, v in setaside.most_common(3))})',
               ''] + table(schemes, st) + ['']
        md += blk; print('\n'.join(blk))
    if len(jobs) > 1:
        blk = ['## POOLED', f'{tot_kept} non-recombinant reads across {len(jobs)} day-0 populations', ''] \
              + table(schemes, pooled) + ['']
        for s in schemes:
            top = pooled_conf[s].most_common(8)
            if not top: continue
            blk += [f'### where {s} sends the foreign-looking copies', '',
                    '| anchored end | apparent donor group lives at | copies |', '|---|---|---|'] \
                   + [f'| {ce} | {dst} | {n} |' for (ce, dst), n in top] + ['']
        md = md[:5] + blk + md[5:]
        print('\n'.join(blk))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        p = os.path.join(out_dir, 'grouping_benchmark.md')
        open(p, 'w').write('\n'.join(md) + '\n')
        print('written', p)
