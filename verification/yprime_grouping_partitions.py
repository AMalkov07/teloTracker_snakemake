#!/usr/bin/env python3
"""Build every Y'-grouping scheme we can compare, as partitions of the SAME element set.

Schemes (all keyed by element `E-<chr_end>-<n>`):
  element         one group per reference Y' element -- the ceiling, nothing can do better
  condensed       the pipeline's own pre-clustering step: homopolymer-condense, then merge
                  anything >= 99.9 % similar (its "entry" level)
  silhouette      the pipeline default: average linkage on the condensed set, k chosen by
                  the best silhouette score
  cut97 / cut99   the same linkage cut at a fixed similarity instead (fcluster t = 100 - x),
                  i.e. what `--stop-mode threshold --identity-threshold 97|99` produces
  curated_variant the curated RepeatMasker library's colour shades (ID2_Red-Light)
  curated_family  the same collapsed to the ID family (ID2)

Similarity is the pipeline's own coverage-penalised formula
(pident x aln_len / (aln_len + max(unaligned_q, unaligned_s))), so a Short Y' scores ~55 %
against a Long one however identical their shared part is -- worth remembering when reading
a "97 % threshold".

Usage: yprime_grouping_partitions.py <elements.fasta> <curated_lib.fasta> <out_dir> [--threads N]
"""
import os, re, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
from cluster_yprimes_paper_method import (run_blast, calc_similarity_matrix, deduplicate_yprimes,
                                          condense_homopolymers, find_clusters_silhouette,
                                          find_clusters_threshold)

elem_fa, cur_fa, out = sys.argv[1], sys.argv[2], sys.argv[3]
threads = int(sys.argv[sys.argv.index('--threads') + 1]) if '--threads' in sys.argv else 8
os.makedirs(out, exist_ok=True)

# ---- read the element library -------------------------------------------------------
entries, names, ids = [], [], []
hdr, buf = None, []
def _flush():
    if hdr is not None:
        entries.append((hdr, ''.join(buf)))
        names.append(hdr.split('#')[0].lstrip('>').replace('Y_Prime_', ''))
        ids.append(hdr.split('/')[-1].strip())
for line in open(elem_fa):
    if line.startswith('>'): _flush(); hdr, buf = line.strip(), []
    else: buf.append(line.strip())
_flush()
N = len(entries)
print(f'{N} elements')

# ---- similarity, exactly as the pipeline computes it --------------------------------
cond = os.path.join(out, '_condensed.fasta')
with open(cond, 'w') as fh:
    for (h, s), nm in zip(entries, names):
        fh.write(f'>{nm}\n{condense_homopolymers(s, max_len=4)}\n')
lens = [len(s) for _, s in entries]
sim = calc_similarity_matrix(run_blast(cond, threads=threads), names, lens)
np.savetxt(os.path.join(out, 'element_similarity.tsv'),
           sim, delimiter='\t', fmt='%.4f', header='\t'.join(names), comments='')

# ---- partitions ---------------------------------------------------------------------
part = {}
part['element'] = {ids[i]: ids[i] for i in range(N)}

# condensed: the pipeline's 99.9 % dedup over the homopolymer-condensed sequences
cond_entries = [(h, condense_homopolymers(s, max_len=4)) for h, s in entries]
_, dedup_names, orig_to_rep = deduplicate_yprimes(cond_entries, sim, names, threshold=99.9)
rep_of = {names[i]: orig_to_rep[names[i]] for i in range(N)}
part['condensed'] = {ids[i]: f'C:{rep_of[names[i]]}' for i in range(N)}

# cluster the deduped set, then push labels back down to every element
keep = [names.index(n) for n in dedup_names]
sub = sim[np.ix_(keep, keep)]
dist = 100.0 - sub
np.fill_diagonal(dist, 0.0)
dist[dist < 0] = 0.0
rep_index = {n: i for i, n in enumerate(dedup_names)}

def push(labels, tag):
    lab = {dedup_names[i]: f'{tag}{labels[i]}' for i in range(len(dedup_names))}
    return {ids[i]: lab[rep_of[names[i]]] for i in range(N)}

k, Z, scores, labels = find_clusters_silhouette(dist, linkage_method='average')
part['silhouette'] = push(labels, 'S')
print(f'silhouette: k={k}  best score={scores.get(k, float("nan")):.3f}')
for thr, tag in ((97.0, 'T97-'), (99.0, 'T99-')):
    nc, _, _, lb = find_clusters_threshold(dist, linkage_method='average', threshold=thr)
    part[f'cut{int(thr)}'] = push(lb, tag)
    print(f'threshold {thr}: {nc} clusters on the deduped set')

# ---- curated, mapped onto the same elements by position -----------------------------
def cur_map(level):
    m = {}
    for line in open(cur_fa):
        if not line.startswith('>'): continue
        name_part, cls = (line.strip().lstrip('>').split('#', 1) + [''])[:2]
        grp = cls.split('/')[2] if len(cls.split('/')) >= 3 else ''
        if level == 'family': grp = grp.split('_', 1)[0]
        for chunk in name_part.replace('Y_Prime_', '').split(';'):
            mm = re.match(r'^(chr\w+?[LR])((?:\d+)(?:,\d+)*)$', chunk)
            if not mm: continue
            for n in mm.group(2).split(','): m[f'E-{mm.group(1)}-{n}'] = grp
    return m
for lvl in ('variant', 'family'):
    m = cur_map(lvl)
    missing = [i for i in ids if i not in m]
    if missing: print(f'  curated_{lvl}: {len(missing)} element(s) absent from the curated library: {missing}')
    part[f'curated_{lvl}'] = {i: m.get(i, 'UNASSIGNED') for i in ids}

json.dump({'elements': ids, 'partitions': part}, open(os.path.join(out, 'partitions.json'), 'w'), indent=1)
print(f'\n{"scheme":<18}{"groups":>7}   sizes')
for nm, p in part.items():
    from collections import Counter
    c = Counter(p.values())
    print(f'{nm:<18}{len(c):>7}   {sorted(c.values(), reverse=True)}')
print('\nwritten to', out)
