#!/usr/bin/env python3
"""
yprime_variants.py -- find recombinant Y' variants: copies built from pieces of two or more day-0 Y's.

RepeatMasker gives every Y' copy in a read one best library element. A copy that is really a block
recombinant -- the 5' part copied from one Y', the rest from another -- is then silently named after
whichever parent covers more of it. This step makes no change to those labels; it re-reads every Y'
copy and reports the recombinants.

Why not "check every copy below 99.9 % identity": a nanopore copy of an unchanged Y' is usually only
98-99.5 % identical to its element, and a recombinant can be 99.2 % identical to one parent, so an
identity gate either checks everything or misses them. And comparing a copy against just two
elements is misleading: a third element that shares some alleles with each looks like an SNP mix of
the two. So every copy is compared against ALL of the strain's Y's at once, cheaply:

  1. Templates   every chr*_Y_Prime_n of the day-0 reference, split into size classes (short / long);
                 templates within a few edits of each other are merged (tandem copies).
  2. Per copy    the copy is aligned to every template of its class (about ten); per copy position
                 each alignment says match or not. Informative positions are where the templates
                 disagree -- a read error costs every template the same, so it never is one -- and
                 homopolymer context is skipped. A switch-penalised path through the templates over
                 those positions gives the copy's segments. A segment needs >= MIN_SUPPORT positions
                 that separate it from its neighbour, else it is merged back; each segment is labelled
                 with every template that explains it equally well ('chr2L1|chr6L1').
  3. Confirm     single reads are noisy, so a recombinant is only reported when the same template
                 path recurs (>= MIN_COPIES copies at >= MIN_ENDS ends); its consensus sequence is then
                 re-called, which gives the switch points without read errors.

Usage:
    python yprime_variants.py <base_name> <recombination_dir> <telomere_reads_dir> \
        --day0-ref REF.fasta --day0-bed REF_simp.bed --out-dir DIR [--threads N]

Writes, in --out-dir:
    <base>_yprime_variants.tsv          one row per confirmed variant (header-only when none)
    <base>_yprime_variant_copies.tsv.gz one row per Y' copy: its call and variant
    <base>_yprime_variants.fasta        variant consensus sequences
"""

import argparse
import collections
import csv
import glob
import gzip
import os
import re
import sys
from multiprocessing import Pool

import edlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recombination_utils import read_fasta
from extract_yprime_fasta import classify_yprime_size

SWITCH_PENALTY = 4        # DP cost of changing template (in mismatching columns)
MIN_SUPPORT = 3           # columns a segment needs that separate it from its neighbour
MERGE_EDITS = 3           # templates this close are one template
MIN_COVERED = 0.5         # a copy covering less of the frame is 'partial'
MAX_DIVERGENCE = 0.15     # a copy further than this from every template is 'noisy'
EVENT_GAP = 12            # separating positions closer than this are one event
EDGE = 50                 # copy ends: boundary offsets, not template signal
MIN_COPIES, MIN_ENDS = 10, 3
MAX_CONSENSUS_COPIES = 400

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
HOMOPOLYMER = re.compile(r'(A{4,}|C{4,}|G{4,}|T{4,})')
YP_NAME = re.compile(r'^(chr\w+?)_Y_Prime_(\d+)$')
YP_POS = re.compile(r'([^;:]+):(\d+)-(\d+)')


def rc(s):
    return s.translate(COMP)[::-1]


def hp_condense(s):
    return re.sub(r'(.)\1+', r'\1', s)


# ---------------------------------------------------------------------------
# Templates and the per-class site table
# ---------------------------------------------------------------------------

def load_templates(ref_fasta, bed):
    """{'chr2L1': seq} for every chr*_Y_Prime_n in the BED, oriented as the BED strand gives it."""
    genome = {k.split()[0]: v.upper() for k, v in read_fasta(ref_fasta).items()}
    out = {}
    with open(bed) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) < 5:
                continue
            m = YP_NAME.match(c[3])
            if not m or c[0] not in genome:
                continue
            s = genome[c[0]][int(c[1]):int(c[2])]
            out[f'{m.group(1)}{m.group(2)}'] = rc(s) if c[4] == '-' else s
    return out


def merge_templates(templates, max_edits=MERGE_EDITS):
    """Group near-identical templates; name each group by its first member (+N more)."""
    groups = []
    for name in sorted(templates, key=lambda n: (-len(templates[n]), n)):
        s = templates[name]
        for g in groups:
            if abs(len(g['seq']) - len(s)) <= max_edits and \
               edlib.align(s, g['seq'], mode='NW', task='distance', k=max_edits)['editDistance'] != -1:
                g['members'].append(name)
                break
        else:
            groups.append({'seq': s, 'members': [name]})
    for g in groups:
        m = sorted(g['members'], key=natural_key)
        g['members'] = m
        g['name'] = m[0] if len(m) == 1 else f'{m[0]}+{len(m) - 1}'
    return groups


def natural_key(name):
    m = re.match(r'chr(\d+)([LR])(\d+)', name)
    return (int(m.group(1)), m.group(2), int(m.group(3))) if m else (999, name, 0)


def error_profile(copy, template):
    """Align copy inside template (template ends free). Per copy position: 1 when that base is a
    mismatch or an insertion against the template, or a template base is missing just before it.
    Returns (profile, edit distance, template start, template end)."""
    r = edlib.align(copy, template, mode='HW', task='path')
    e = [0] * (len(copy) + 1)
    qi = 0
    for n, op in re.findall(r'(\d+)([=XID])', r['cigar']):
        n = int(n)
        if op == '=':
            qi += n
        elif op in 'XI':
            for k in range(n):
                e[qi + k] = 1
            qi += n
        else:
            e[qi] = 1
    t0, t1 = r['locations'][0]
    return e, r['editDistance'], t0, t1 + 1


class SizeClass:
    """The merged templates of one Y' size class; the medoid picks a copy's class cheaply."""

    def __init__(self, name, groups):
        self.name = name
        self.groups = groups
        self.names = [g['name'] for g in groups]
        self.seqs = [g['seq'] for g in groups]
        fi = 0 if len(self.seqs) == 1 else min(range(len(self.seqs)), key=lambda i: sum(
            edlib.align(self.seqs[i], s, mode='NW', task='distance')['editDistance'] for s in self.seqs))
        self.frame = self.seqs[fi]


def profiles(copy, c):
    """Per template: the copy's error profile. Informative positions are copy positions where the
    templates disagree (homopolymer context skipped); a read error costs every template the same
    and so never becomes informative. Returns (positions, cost[template][position], per-template
    (distance, start, end))."""
    prof = [error_profile(copy, s) for s in c.seqs]
    hp = [False] * (len(copy) + 1)
    for m in HOMOPOLYMER.finditer(copy):
        for q in range(max(0, m.start() - 1), min(len(copy) + 1, m.end() + 1)):
            hp[q] = True
    L = len(copy)
    pos = [q for q in range(EDGE, L + 1 - EDGE) if not hp[q] and len({p[0][q] for p in prof}) > 1]
    cost = [[p[0][q] for q in pos] for p in prof]
    return pos, cost, [p[1:] for p in prof]


# ---------------------------------------------------------------------------
# Template path through a copy
# ---------------------------------------------------------------------------

def best_path(cost, switch=SWITCH_PENALTY):
    """Min-cost template per informative position, SWITCH per change of template."""
    K = len(cost)
    n = len(cost[0]) if K else 0
    if not n:
        return [0] * 0
    D = [cost[k][0] for k in range(K)]
    back = []
    for x in range(1, n):
        b = min(range(K), key=D.__getitem__)
        nd, bk = [], []
        for k in range(K):
            if D[k] <= D[b] + switch:
                nd.append(D[k] + cost[k][x]); bk.append(k)
            else:
                nd.append(D[b] + switch + cost[k][x]); bk.append(b)
        D = nd
        back.append(bk)
    k = min(range(K), key=D.__getitem__)
    path = [k]
    for bk in reversed(back):
        k = bk[k]
        path.append(k)
    return path[::-1]


def segments(path, cost, pos, min_support=MIN_SUPPORT, tie_tol=0):
    """Collapse the path into segments (lists of informative-position indices); merge back any
    segment with fewer than min_support separating events against a neighbour (positions within
    EVENT_GAP bp are one event: an indel or a moved segment is one difference, not one per base);
    label each by every template within tie_tol mismatches of the chosen one."""
    segs = []
    for x, k in enumerate(path):
        if segs and segs[-1]['k'] == k:
            segs[-1]['x'].append(x)
        else:
            segs.append({'k': k, 'x': [x]})

    def sep(s, n):     # positions where s's template matches and n's does not
        return [x for x in s['x'] if cost[s['k']][x] == 0 and cost[n][x] == 1]

    def events(xs):
        return sum(1 for a, b in zip([None] + xs, xs) if a is None or pos[b] - pos[a] > EVENT_GAP)

    def support(i):
        nb = [segs[j]['k'] for j in (i - 1, i + 1) if 0 <= j < len(segs)]
        return min(events(sep(segs[i], n)) for n in nb) if nb else 10 ** 9

    while len(segs) > 1:
        i = min(range(len(segs)), key=support)
        if support(i) >= min_support:
            break
        nbrs = [j for j in (i - 1, i + 1) if 0 <= j < len(segs)]
        t = min(nbrs, key=lambda j: sum(cost[segs[j]['k']][x] for x in segs[i]['x']))
        segs[t]['x'] = sorted(segs[t]['x'] + segs[i]['x'])
        del segs[i]
        merged = [segs[0]]
        for s in segs[1:]:
            if s['k'] == merged[-1]['k']:
                merged[-1]['x'] += s['x']
            else:
                merged.append(s)
        segs = merged
    for i, s in enumerate(segs):
        if i:
            p = segs[i - 1]
            s['switch'] = (max(sep(p, s['k']), default=p['x'][-1]), min(sep(s, p['k']), default=s['x'][0]))
    n = len(cost[0])
    for i, s in enumerate(segs):
        # ties are judged only where the segment is certain: between its switch intervals (inside an
        # interval the switch could sit on either side of a position)
        lo = s['switch'][1] if i else 0
        hi = segs[i + 1]['switch'][0] if i + 1 < len(segs) else n - 1
        core = range(lo, hi + 1)
        mm = [sum(cost[k][x] for x in core) for k in range(len(cost))]
        s['mismatch'] = mm[s['k']]
        s['ties'] = [k for k in range(len(cost)) if mm[k] <= mm[s['k']] + tie_tol]
    return segs


def best_class(copy, classes):
    best = None
    for c in classes:
        d = edlib.align(copy, c.frame, mode='HW', task='distance', k=int(0.3 * len(copy)) + 1)['editDistance']
        if d >= 0 and (best is None or d < best[0]):
            best = (d, c)
    return best[1] if best else None


def call_copy(seq, classes, tie_tol=1):
    """Call one oriented Y' copy: pure / recomb / partial / noisy, with its template signature."""
    c = best_class(seq, classes)
    if c is None or not seq:
        return {'call': 'noisy', 'signature': '', 'class': '', 'covered': 0.0}
    pos, cost, stats = profiles(seq, c)
    kb = min(range(len(stats)), key=lambda k: stats[k][0])
    d, t0, t1 = stats[kb]
    res = {'class': c.name, 'covered': round((t1 - t0) / len(c.seqs[kb]), 3), 'informative': len(pos),
           'divergence': round(d / max(1, len(seq)), 4)}
    if res['covered'] < MIN_COVERED:
        return dict(res, call='partial', signature='')
    if res['divergence'] > MAX_DIVERGENCE:
        return dict(res, call='noisy', signature='')
    if not pos:
        return dict(res, call='pure', signature=c.names[kb])
    segs = segments(best_path(cost), cost, pos, tie_tol=tie_tol)
    sig = '>'.join('|'.join(c.names[k] for k in s['ties']) for s in segs)
    return dict(res, call='pure' if len(segs) == 1 else 'recomb', signature=sig)


# ---------------------------------------------------------------------------
# Reads -> copies
# ---------------------------------------------------------------------------

def iter_copies(features_tsv, reads_fasta):
    """Yield (read_id, chr_end, copy_index, label, start, end, oriented copy sequence)."""
    end = re.search(r'_(chr\d+[LR])_features\.tsv$', features_tsv).group(1)
    with open(features_tsv) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter='\t') if r.get('y_prime_positions')]
    if not rows or not os.path.exists(reads_fasta):
        return
    want = {r['read_id'] for r in rows}
    reads = {}
    for h, s in read_fasta(reads_fasta).items():
        rid = h.split()[0]
        if rid in want:
            reads[rid] = s.upper()
    for r in rows:
        s = reads.get(r['read_id'])
        if not s:
            continue
        L = len(s)
        flip = r.get('telo_side') == 'beginning'
        o = rc(s) if flip else s
        for i, m in enumerate(YP_POS.finditer(r['y_prime_positions'])):
            a, b = int(m.group(2)), int(m.group(3))
            oa, ob = (L - b, L - a) if flip else (a, b)
            yield r['read_id'], end, i, m.group(1), a, b, o[oa:ob]


_CLASSES = None


def _init(classes):
    global _CLASSES
    _CLASSES = classes


def _call_end(job):
    features_tsv, reads_fasta = job
    out = []
    for rid, end, i, lab, a, b, seq in iter_copies(features_tsv, reads_fasta):
        c = call_copy(seq, _CLASSES)
        out.append({'read_id': rid, 'chr_end': end, 'copy_index': i, 'label': lab, 'start': a, 'end': b,
                    'call': c['call'], 'signature': c['signature'], 'class': c['class'],
                    'covered': c['covered'], 'informative': c.get('informative', 0),
                    'divergence': c.get('divergence', ''), 'seq': seq if c['call'] == 'recomb' else ''})
    return out


# ---------------------------------------------------------------------------
# Consensus confirmation
# ---------------------------------------------------------------------------

def consensus(seqs, frame):
    """Majority base per frame position; an insertion is kept when > half the copies carry one."""
    C = [collections.Counter() for _ in range(len(frame))]
    I = [collections.Counter() for _ in range(len(frame) + 1)]
    for s in seqs:
        r = edlib.align(s, frame, mode='NW', task='path')
        qi = ti = 0
        ins = ''
        for n, op in re.findall(r'(\d+)([=XID])', r['cigar']):
            n = int(n)
            if op == 'I':
                ins += s[qi:qi + n]
                qi += n
                continue
            for k in range(n):
                I[ti + k][ins if k == 0 else ''] += 1
                C[ti + k][s[qi + k] if op in '=X' else '-'] += 1
            ins = ''
            if op in '=X':
                qi += n
            ti += n
        I[len(frame)][ins] += 1
    n = len(seqs)
    out = []
    for i in range(len(frame) + 1):
        if n - I[i][''] > n / 2:
            out.append(collections.Counter({x: c for x, c in I[i].items() if x}).most_common(1)[0][0])
        if i < len(frame):
            b = C[i].most_common(1)[0][0]
            if b != '-':
                out.append(b)
    return ''.join(out)


def describe(seq, c):
    """Re-call a consensus sequence: segments with their span along it (1-based) and the switch
    intervals (last base that still follows the left template, first base that follows the right)."""
    pos, cost, stats = profiles(seq, c)
    segs = segments(best_path(cost), cost, pos) if pos else []
    parts = []
    for i, s in enumerate(segs):
        lo = 1 if i == 0 else pos[s['switch'][1]] + 1
        hi = len(seq) if i == len(segs) - 1 else pos[segs[i + 1]['switch'][0]] + 1
        parts.append(('|'.join(c.names[k] for k in s['ties']), lo, hi, s['mismatch']))
    switches = [f"{pos[s['switch'][0]] + 1}-{pos[s['switch'][1]] + 1}" for s in segs[1:]]
    near = sorted((edlib.align(hp_condense(seq), hp_condense(g['seq']), mode='NW', task='distance')['editDistance'],
                   g['name']) for g in c.groups)
    return {'signature': '>'.join(p[0] for p in parts),
            'segments': '; '.join(f'{p[0]} [{p[1]}-{p[2]}]' for p in parts),
            'switches': ', '.join(switches), 'n_segments': len(parts),
            'unexplained': sum(p[3] for p in parts),
            'nearest_template': near[0][1], 'nearest_edits_hp': near[0][0]}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

VARIANT_COLS = ['variant', 'class', 'signature', 'segments', 'switch_intervals', 'n_segments', 'copies',
                'reads', 'ends', 'end_list', 'pipeline_labels', 'pct_of_label_copies', 'consensus_len',
                'nearest_template', 'nearest_edits_hp', 'unexplained_columns', 'per_copy_signatures']
COPY_COLS = ['read_id', 'chr_end', 'copy_index', 'label', 'start', 'end', 'class', 'call', 'signature',
             'covered', 'informative', 'divergence', 'variant']


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('base_name')
    ap.add_argument('recombination_dir')
    ap.add_argument('telomere_reads_dir')
    ap.add_argument('--day0-ref', required=True)
    ap.add_argument('--day0-bed', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--min-copies', type=int, default=MIN_COPIES)
    ap.add_argument('--min-ends', type=int, default=MIN_ENDS)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    base = args.base_name
    out_tsv = os.path.join(args.out_dir, f'{base}_yprime_variants.tsv')
    out_copies = os.path.join(args.out_dir, f'{base}_yprime_variant_copies.tsv.gz')
    out_fa = os.path.join(args.out_dir, f'{base}_yprime_variants.fasta')

    templates = load_templates(args.day0_ref, args.day0_bed)
    by_class = collections.defaultdict(dict)
    for n, s in templates.items():
        by_class[classify_yprime_size(len(s))][n] = s
    classes = [SizeClass(k, merge_templates(v)) for k, v in sorted(by_class.items()) if v]
    for c in classes:
        print(f'  {c.name}: {len(c.groups)} templates (from {sum(len(g["members"]) for g in c.groups)} Y\'), '
              f'lengths {min(map(len, c.seqs))}-{max(map(len, c.seqs))} bp')

    jobs = []
    for f in sorted(glob.glob(os.path.join(args.recombination_dir, f'{base}_chr*_features.tsv'))):
        end = re.search(r'_(chr\d+[LR])_features\.tsv$', f).group(1)
        jobs.append((f, os.path.join(args.telomere_reads_dir, f'{base}_{end}_telomere_reads.fasta')))
    copies = []
    if classes and jobs:
        if args.threads > 1:
            with Pool(args.threads, initializer=_init, initargs=(classes,)) as pool:
                for res in pool.imap_unordered(_call_end, jobs):
                    copies += res
        else:
            _init(classes)
            for j in jobs:
                copies += _call_end(j)
    calls = collections.Counter(c['call'] for c in copies)
    print(f'  {len(copies)} Y\' copies: ' + ', '.join(f'{k} {v}' for k, v in calls.most_common()))

    # recurring per-copy signatures -> consensus -> re-call; merge groups whose consensus calls agree
    cls = {c.name: c for c in classes}
    groups = collections.defaultdict(list)
    for c in copies:
        if c['call'] == 'recomb':
            groups[(c['class'], c['signature'])].append(c)
    confirmed = {}
    for (cname, sig), members in groups.items():
        if len(members) < args.min_copies or len({m['chr_end'] for m in members}) < args.min_ends:
            continue
        full = [m['seq'] for m in members if m['covered'] >= 0.95][:MAX_CONSENSUS_COPIES]
        if len(full) < args.min_copies:
            continue
        c = cls[cname]
        frame = min(c.seqs, key=lambda t: edlib.align(full[0], t, mode='HW', task='distance')['editDistance'])
        cons = consensus(full, frame)
        d = describe(cons, c)
        if d['n_segments'] < 2:
            continue      # the consensus is one template: the per-copy switches were read noise
        key = (cname, d['signature'])
        v = confirmed.setdefault(key, {'class': cname, 'desc': d, 'cons': cons, 'members': [], 'sigs': []})
        v['members'] += members
        v['sigs'].append(sig)
        if len(members) > len(v['members']) - len(members):   # keep the description of the biggest group
            v['desc'], v['cons'] = d, cons

    label_copies = collections.Counter(c['label'] for c in copies)
    rows, fasta = [], []
    for i, v in enumerate(sorted(confirmed.values(), key=lambda x: -len(x['members'])), 1):
        vid = f'V{i}'
        m = v['members']
        for c in m:
            c['variant'] = vid
        labs = collections.Counter(c['label'] for c in m)
        pct = 100 * len(m) / sum(label_copies[l] for l in labs)
        ends = collections.Counter(c['chr_end'] for c in m)
        d = v['desc']
        rows.append({'variant': vid, 'class': v['class'], 'signature': d['signature'], 'segments': d['segments'],
                     'switch_intervals': d['switches'], 'n_segments': d['n_segments'], 'copies': len(m),
                     'reads': len({c['read_id'] for c in m}), 'ends': len(ends),
                     'end_list': ','.join(e for e, _ in ends.most_common()),
                     'pipeline_labels': ','.join(f'{l}:{n}' for l, n in labs.most_common()),
                     'pct_of_label_copies': f'{pct:.1f}', 'consensus_len': len(v['cons']),
                     'nearest_template': d['nearest_template'], 'nearest_edits_hp': d['nearest_edits_hp'],
                     'unexplained_columns': d['unexplained'], 'per_copy_signatures': ' ; '.join(v['sigs'])})
        fasta.append((f'{base}|{vid}|{d["signature"]}|copies={len(m)}|ends={len(ends)}', v['cons']))

    with open(out_tsv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, VARIANT_COLS, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    with gzip.open(out_copies, 'wt', newline='') as fh:
        w = csv.DictWriter(fh, COPY_COLS, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for c in sorted(copies, key=lambda c: (c['chr_end'], c['read_id'], c['copy_index'])):
            c.setdefault('variant', '')
            w.writerow(c)
    with open(out_fa, 'w') as fh:
        for h, s in fasta:
            fh.write(f'>{h}\n{s}\n')
    print(f'  {len(rows)} confirmed recombinant variant(s) -> {out_tsv}')
    for r in rows:
        print(f"    {r['variant']}: {r['segments']}  ({r['copies']} copies, {r['ends']} ends, "
              f"{r['pct_of_label_copies']} % of {r['pipeline_labels'].split(':')[0]}... copies)")


if __name__ == '__main__':
    main()
