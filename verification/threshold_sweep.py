#!/usr/bin/env python3
"""Mismatch rate versus grouping threshold, one reference at a time.

For each day-0 reference: build the Y' grouping at a series of identity cutoffs (the pipeline's
`--stop-mode threshold` path), then score every read whose Y' copy count equals its anchor
end's reference array length. A copy is a MISMATCH when the group of the element it matched
differs from the group of the element that positionally belongs there.

Reported per threshold: number of groups, per-copy mismatch %, and per-read mismatch % (a read
counts once if any of its copies mismatched -- this is what the pipeline would surface as a
recombination call).

The all-vs-all BLAST is done once per reference and reused across thresholds; only the
dendrogram cut changes.

NOTE ON "ERROR": at day 0 a mismatch is not automatically an error. Some are genuine
recombinants or standing variation in the culture. These rates are an UPPER BOUND on the
error rate, not a measurement of it.

Usage: threshold_sweep.py --libs <dir of elem_<sample>_fixed.fasta> --results-root <dir>
                          [--thresholds 95 96 97 98 99 99.5] [--out <tsv>]
"""
import argparse, csv, glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
from cluster_yprimes_paper_method import (run_blast, calc_similarity_matrix, deduplicate_yprimes,
                                          condense_homopolymers, find_clusters_threshold,
                                          find_clusters_silhouette)

p = argparse.ArgumentParser()
p.add_argument('--libs', default='verification/grouping_benchmark')
p.add_argument('--results-root', default='/home/andrey/argon_scratch/telo_sra_runs/results')
p.add_argument('--tag', default='elemYPfix')
p.add_argument('--thresholds', nargs='+', type=float, default=[95, 96, 97, 98, 99, 99.5])
p.add_argument('--threads', type=int, default=8)
p.add_argument('--curated-root', default='verification/curated_refs',
               help='dir holding <strain>_features/repeatmasker_<strain>_all_y_primes.fasta')
p.add_argument('--out', default='verification/reports/threshold_sweep.tsv')
a = p.parse_args()


def curated_map(strain, level, curated_root):
    """element id -> curated group, mapped by array position from the curated library headers.

    Headers look like  >Y_Prime_chr4R1,2,3,6,7;chr12R6,7#Long/Tandem/ID2_Red-Light
    so one entry can name several positions across several ends. `level` is 'variant'
    (ID2_Red-Light) or 'family' (ID2).
    """
    fa = os.path.join(curated_root, f'{strain}_features',
                      f'repeatmasker_{strain}_all_y_primes.fasta')
    if not os.path.exists(fa): return None
    m = {}
    for line in open(fa):
        if not line.startswith('>'): continue
        name_part, cls = (line.strip().lstrip('>').split('#', 1) + [''])[:2]
        parts = cls.split('/')
        grp = parts[2] if len(parts) >= 3 else ''
        if level == 'family': grp = grp.split('_', 1)[0]
        for chunk in name_part.replace('Y_Prime_', '').split(';'):
            mm = re.match(r'^(chr\w+?[LR])((?:\d+)(?:,\d+)*)$', chunk)
            if not mm: continue
            for n in mm.group(2).split(','):
                m[f'E-{mm.group(1)}-{n}'] = grp
    return m


def load_lib(fa):
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
    return entries, names, ids


rows = []
for fa in sorted(glob.glob(os.path.join(a.libs, 'elem_*_fixed.fasta'))):
    sample = os.path.basename(fa).replace('elem_', '').replace('_fixed.fasta', '')
    entries, names, ids = load_lib(fa)
    N = len(entries)

    # one BLAST + similarity matrix per reference, reused for every threshold
    cond = f'/tmp/_sweep_{sample}.fasta'
    with open(cond, 'w') as fh:
        for (h, s), nm in zip(entries, names):
            fh.write(f'>{nm}\n{condense_homopolymers(s, max_len=4)}\n')
    lens = [len(s) for _, s in entries]
    sim = calc_similarity_matrix(run_blast(cond, threads=a.threads), names, lens)
    os.remove(cond)
    cond_entries = [(h, condense_homopolymers(s, max_len=4)) for h, s in entries]
    _, dedup_names, orig_to_rep = deduplicate_yprimes(cond_entries, sim, names, threshold=99.9)
    rep_of = {names[i]: orig_to_rep[names[i]] for i in range(N)}
    keep = [names.index(n) for n in dedup_names]
    dist = 100.0 - sim[np.ix_(keep, keep)]
    np.fill_diagonal(dist, 0.0); dist[dist < 0] = 0.0

    def push(labels, tag):
        lab = {dedup_names[i]: f'{tag}{labels[i]}' for i in range(len(dedup_names))}
        return {ids[i]: lab[rep_of[names[i]]] for i in range(N)}

    schemes = {}
    for t in a.thresholds:
        _, _, _, lb = find_clusters_threshold(dist, linkage_method='average', threshold=t)
        schemes[f'{t:g}'] = push(lb, 'T')
    k, _, _, lb = find_clusters_silhouette(dist, linkage_method='average')
    schemes['silhouette'] = push(lb, 'S')

    # curated schemes, mapped onto the same elements by array position
    strain = sample.split('_')[0]
    for level, label in (('variant', 'curated_variant'), ('family', 'curated_family')):
        cm = curated_map(strain, level, a.curated_root)
        if cm is None: continue
        missing = [i for i in ids if i not in cm]
        if missing:
            print(f'  {sample} {label}: {len(missing)} element(s) absent from the curated '
                  f'library ({", ".join(missing[:4])}) -> grouped as UNASSIGNED', file=sys.stderr)
        schemes[label] = {i: cm.get(i, 'UNASSIGNED') for i in ids}

    # read the sample's element-level calls once
    end_of = {i: re.match(r'^E-(chr\w+?[LR])-\d+$', i).group(1) for i in ids}
    idx_of = {i: int(i.rsplit('-', 1)[1]) for i in ids}
    exp = {}
    for i in sorted(ids, key=lambda x: (end_of[x], idx_of[x])): exp.setdefault(end_of[i], []).append(i)
    reads = []
    for f in glob.glob(f'{a.results_root}/{sample}__{a.tag}/_pipeline/recombination/{sample}_chr*_features.tsv'):
        for r in csv.DictReader(open(f), delimiter='\t'):
            ce = r['chr_end']; obs = [x for x in r.get('y_prime_observed_array', '').split(',') if x]
            if ce not in exp or not obs or len(obs) != len(exp[ce]): continue
            if any(o not in end_of for o in obs): continue
            reads.append((ce, obs, exp[ce]))
    if not reads:
        print(f'  {sample}: no scored reads', file=sys.stderr); continue

    for name, g in schemes.items():
        copies = bad_copies = bad_reads = 0
        for ce, obs, truth in reads:
            bad = 0
            for o, t in zip(obs, truth):
                copies += 1
                if g[o] != g[t]: bad += 1
            bad_copies += bad
            if bad: bad_reads += 1
        rows.append({'sample': sample, 'scheme': name,
                     'n_groups': len(set(g.values())), 'n_elements': N,
                     'reads_scored': len(reads), 'copies_scored': copies,
                     'mismatch_copies': bad_copies,
                     'copy_mismatch_pct': round(100 * bad_copies / copies, 3),
                     'mismatch_reads': bad_reads,
                     'read_mismatch_pct': round(100 * bad_reads / len(reads), 3)})
    print(f'  {sample}: {len(reads)} reads, {N} elements', file=sys.stderr)

os.makedirs(os.path.dirname(a.out), exist_ok=True)
with open(a.out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t'); w.writeheader()
    for r in rows: w.writerow(r)
print(f'written {a.out}', file=sys.stderr)
