#!/usr/bin/env python3
"""Resolution of each Y'-grouping scheme, from the reference alone (no reads needed).

For every scheme in partitions.json this reports
  * the number of groups, and how far that is from the sequence ceiling (`condensed`)
  * how many reference elements the scheme leaves indistinguishable from another element
    at a DIFFERENT chromosome end -- the ones that can send a read to the wrong donor
  * how many chromosome ends carry a unique array signature (the ordered string of group
    labels for that end's Y' array), which is what source attribution actually reads
Usage: yprime_grouping_resolution.py <partitions.json> [<out.md>]
"""
import json, re, sys
from collections import Counter, defaultdict

P = json.load(open(sys.argv[1]))
ids, parts = P['elements'], P['partitions']
end_of = {i: re.match(r'^E-(chr\w+?[LR])-\d+$', i).group(1) for i in ids}
idx_of = {i: int(i.rsplit('-', 1)[1]) for i in ids}
ends = sorted({end_of[i] for i in ids}, key=lambda c: (int(re.search(r'\d+', c).group()), c))
ORDER = ['element', 'condensed', 'cut99', 'cut97', 'curated_variant', 'silhouette', 'curated_family']
schemes = [s for s in ORDER if s in parts] + [s for s in parts if s not in ORDER]

def signature(p, ce):
    els = sorted([i for i in ids if end_of[i] == ce], key=lambda i: idx_of[i])
    return ','.join(p[i] for i in els)

rows = []
for name in schemes:
    p = parts[name]
    ngroups = len(set(p.values()))
    # elements sharing a group with an element at another chromosome end
    cross = 0
    bygroup = defaultdict(set)
    for i in ids: bygroup[p[i]].add(end_of[i])
    for i in ids:
        if len(bygroup[p[i]]) > 1: cross += 1
    sigs = {ce: signature(p, ce) for ce in ends}
    c = Counter(sigs.values())
    uniq = sum(1 for ce in ends if c[sigs[ce]] == 1)
    rows.append({'scheme': name, 'groups': ngroups, 'cross_end_ambiguous_elements': cross,
                 'ends_with_unique_signature': uniq, 'of_ends': len(ends),
                 'distinct_signatures': len(c)})

W = [('scheme', 18), ('groups', 7), ('cross_end_ambiguous_elements', 30),
     ('ends_with_unique_signature', 27), ('distinct_signatures', 20)]
print(''.join(h.ljust(w) for h, w in W))
for r in rows:
    print(''.join(str(r[h]).ljust(w) for h, w in W))

md = ['| scheme | groups | elements ambiguous across ends | ends with a unique array signature | distinct signatures |',
      '|---|---|---|---|---|']
for r in rows:
    md.append(f"| {r['scheme']} | {r['groups']} | {r['cross_end_ambiguous_elements']}/{len(ids)} | "
              f"{r['ends_with_unique_signature']}/{r['of_ends']} | {r['distinct_signatures']} |")

md += ['', '## Array signature per chromosome end', '',
       '| end | ' + ' | '.join(schemes) + ' |', '|---' * (len(schemes) + 1) + '|']
for ce in ends:
    md.append('| ' + ce + ' | ' + ' | '.join(f'`{signature(parts[s], ce)}`' for s in schemes) + ' |')

if len(sys.argv) > 2:
    open(sys.argv[2], 'w').write('\n'.join(md) + '\n')
    print('\nwritten', sys.argv[2])
