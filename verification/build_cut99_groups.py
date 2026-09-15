#!/usr/bin/env python3
"""Rebuild the pipeline's fixed-threshold Y' grouping for one element-level library.

Reproduces exactly what cluster_yprimes_paper_method.py does with
`--stop-mode threshold --identity-threshold <X>`: homopolymer-condense, all-vs-all blastn,
the coverage-penalised similarity (pident x aln_len / (aln_len + max(unaligned))), 99.9 %
dedup, average linkage, then fcluster at t = 100 - X. Labels are pushed from the deduped
representatives back onto every element.

Usage: build_cut99_groups.py <elements.fasta> <out.json> [--threshold 99.0] [--threads 8]
"""
import json, os, re, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
from cluster_yprimes_paper_method import (run_blast, calc_similarity_matrix, deduplicate_yprimes,
                                          condense_homopolymers, find_clusters_threshold)

fa, out = sys.argv[1], sys.argv[2]
thr = float(sys.argv[sys.argv.index('--threshold') + 1]) if '--threshold' in sys.argv else 99.0
threads = int(sys.argv[sys.argv.index('--threads') + 1]) if '--threads' in sys.argv else 8

entries, names, ids = [], [], []
hdr, buf = None, []
def flush():
    if hdr is None: return
    entries.append((hdr, ''.join(buf)))
    names.append(hdr.split('#')[0].lstrip('>').replace('Y_Prime_', ''))
    m = re.search(r'/E-(chr\w+?[LR])-(\d+)$', hdr)
    ids.append(f'E-{m.group(1)}-{m.group(2)}' if m else names[-1])
for line in open(fa):
    if line.startswith('>'): flush(); hdr, buf = line.strip(), []
    else: buf.append(line.strip())
flush()
N = len(entries)

cond_path = out + '.cond.fasta'
with open(cond_path, 'w') as fh:
    for (h, s), nm in zip(entries, names):
        fh.write(f'>{nm}\n{condense_homopolymers(s, max_len=4)}\n')
lens = [len(s) for _, s in entries]
sim = calc_similarity_matrix(run_blast(cond_path, threads=threads), names, lens)
os.remove(cond_path)

cond_entries = [(h, condense_homopolymers(s, max_len=4)) for h, s in entries]
_, dedup_names, orig_to_rep = deduplicate_yprimes(cond_entries, sim, names, threshold=99.9)
rep_of = {names[i]: orig_to_rep[names[i]] for i in range(N)}

keep = [names.index(n) for n in dedup_names]
sub = sim[np.ix_(keep, keep)]
dist = 100.0 - sub
np.fill_diagonal(dist, 0.0); dist[dist < 0] = 0.0
ncl, _, _, labels = find_clusters_threshold(dist, linkage_method='average', threshold=thr)
lab = {dedup_names[i]: f'G{labels[i]}' for i in range(len(dedup_names))}
groups = {ids[i]: lab[rep_of[names[i]]] for i in range(N)}

json.dump({'threshold': thr, 'n_groups': len(set(groups.values())), 'groups': groups},
          open(out, 'w'), indent=1)
print(f'{N} elements -> {len(set(groups.values()))} groups at {thr}% (deduped set: {len(dedup_names)})')
from collections import defaultdict
g2e = defaultdict(list)
for e, g in groups.items(): g2e[g].append(e)
for g in sorted(g2e, key=lambda x: -len(g2e[x])):
    print(f'  {g:<5} n={len(g2e[g]):<3} {", ".join(sorted(e.replace("E-","") for e in g2e[g]))}')
print(f'written {out}')
