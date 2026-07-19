#!/usr/bin/env python3
"""
Prototype: (A) split recombination_confidence from donor_confidence
           (B) make Y' a first-class donor source

Read-only. Recomputes scores from an existing run's per-read features.tsv.
Portable into analyze_features.py (see reconcile_v2 / y_prime_donor_candidates).
"""
import argparse, csv, glob, os, re, sys
from collections import defaultdict

# Feature distinctiveness (same priors as the pipeline)
W = {'spacer': 0.9, 'y_prime': 0.6, 'x_element': 0.4, 'supplementary': 0.25}
ITS_RESOLVABLE_BP = 30      # only trust ITS length differences above nanopore noise
USE_ITS_TIEBREAK = False    # TEMPORARILY OFF: compares the wrong gaps (native, not
                            # gained). Evaluate the Y' matcher unaided first.


# ---------------------------------------------------------------- reference
def build_native_arrays(lib_fasta):
    """chr_end -> ordered list of Y' cluster IDs (position 1..n, anchor->telomere)."""
    pos = defaultdict(dict)
    for line in open(lib_fasta):
        if not line.startswith('>'):
            continue
        h = line[1:].strip()
        name, cls = (h.split('#', 1) + [''])[:2]
        parts = cls.split('/')
        yid = parts[2].split('_', 1)[0] if len(parts) >= 3 else ''
        if not yid:
            continue
        name = re.sub(r'^Y_Prime_', '', name.split()[0])
        for tok in name.split(';'):
            m = re.match(r'(chr\d+[LR])([\d,]+)$', tok)
            if not m:
                continue
            end = m.group(1)
            for p in m.group(2).split(','):
                if p:
                    pos[end][int(p)] = yid
    return {e: [d[k] for k in sorted(d)] for e, d in pos.items()}


def build_its_vectors(bed):
    """chr_end -> ordered list of inter-Y' ITS lengths."""
    its = defaultdict(dict)
    for line in open(bed):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6:
            continue
        m = re.match(r'ITS_(chr\d+[LR])_Y_Prime_(\d+)-(\d+)$', f[3])
        if m:
            its[m.group(1)][(int(m.group(2)), int(m.group(3)))] = int(f[5])
    return {e: [d[k] for k in sorted(d)] for e, d in its.items()}


# ---------------------------------------------------------- (B) Y' donor
def contiguous_match(gained, native):
    """Best contiguous match of `gained` inside `native`.
    Returns (score, is_suffix). Suffix matches score highest: a BIR graft
    replaces the donor's telomere-proximal tail onto the recipient."""
    g, n = len(gained), len(native)
    if g == 0 or n == 0:
        return 0.0, False
    best, best_suffix = 0.0, False
    for start in range(n):
        k = 0
        while start + k < n and k < g and native[start + k] == gained[k]:
            k += 1
        if k == 0:
            continue
        frac = k / g
        is_suffix = (start + k == n)          # match runs to donor's telomere end
        score = frac * (1.0 if is_suffix else 0.85)
        if score > best:
            best, best_suffix = score, is_suffix
    return best, best_suffix


def read_its_lengths(y_prime_positions):
    """Inter-Y' gap lengths on the read, from 'ID:start-end;ID:start-end'."""
    spans = []
    for tok in (y_prime_positions or '').split(';'):
        m = re.match(r'[^:]+:(\d+)-(\d+)$', tok.strip())
        if m:
            spans.append((int(m.group(1)), int(m.group(2))))
    spans.sort()
    return [max(0, spans[i + 1][0] - spans[i][1]) for i in range(len(spans) - 1)]


def y_prime_donor_candidates(chr_end, observed, div_idx, native_arrays, its_vecs,
                             read_its):
    """(B) Rank candidate donor ends for the gained Y' sub-array."""
    if div_idx is None or div_idx < 0 or not observed:
        return [], ''
    gained = observed[div_idx:]
    if not gained:
        return [], 'no_gained_yprimes'          # pure loss: Y' carries no donor info

    scored = []
    for cand, nat in native_arrays.items():
        if cand == chr_end or not nat:
            continue
        s, suf = contiguous_match(gained, nat)
        if s > 0:
            scored.append([cand, round(s, 3), suf])
    if not scored:
        return [], 'no_end_carries_gained_array'

    top = max(s for _, s, _ in scored)
    scored = [c for c in scored if c[1] >= top - 1e-9]      # keep the tied leaders
    note = ''

    # ITS tie-break: only when candidates' tracts differ by a resolvable margin
    if USE_ITS_TIEBREAK and len(scored) > 1 and read_its:
        cand_its = {c[0]: its_vecs.get(c[0], []) for c in scored}
        firsts = [v[0] for v in cand_its.values() if v]
        if firsts and (max(firsts) - min(firsts)) >= ITS_RESOLVABLE_BP:
            def dist(c):
                v = cand_its.get(c, [])
                return abs(v[0] - read_its[0]) if v and read_its else 1e9
            best = min(dist(c[0]) for c in scored)
            keep = [c for c in scored if dist(c[0]) <= best + ITS_RESOLVABLE_BP]
            if len(keep) < len(scored):
                note = 'its_tiebreak'
                scored = keep
    return sorted(scored, key=lambda x: (-x[1], x[0])), note


# ------------------------------------------------- (A) two separate scores
def y_prime_recomb_confidence(status, n_yp, downstream_consistent):
    """Replaces the hardcoded 0.3 floor with real Y'-derived evidence."""
    if status in ('No Change', '', 'no_data'):
        return 0.0
    cleanliness = 1.0 if downstream_consistent else 0.7
    support = min(1.0, n_yp / 2.0) if n_yp else 0.5   # more Y's seen = firmer call
    return W['y_prime'] * cleanliness * support


def noisy_or(vals):
    p = 1.0
    for v in vals:
        p *= (1.0 - max(0.0, min(1.0, v)))
    return 1.0 - p


def reconcile_v2(row, native_arrays, its_vecs):
    chr_end = row['chr_end']
    observed = [x for x in (row.get('y_prime_observed_array') or '').split(',') if x]
    try:
        div_idx = int(row.get('y_prime_divergence_idx', -1))
    except ValueError:
        div_idx = -1
    status = row.get('y_prime_recombination_status', '')
    dsc = str(row.get('y_prime_downstream_consistent', 'True')).lower() == 'true'

    sp_rec = row.get('spacer_recombination', 'no_data')
    x_rec = row.get('x_element_recombination', 'no_data')
    sp_sw = sp_rec in ('switch_detected', 'full_switch')
    x_sw = x_rec in ('switch_detected', 'full_switch')
    y_sw = status not in ('No Change', '', 'no_data')

    def f(k):
        try:
            return float(row.get(k, 0) or 0)
        except ValueError:
            return 0.0
    sp_conf, x_conf = f('spacer_confidence'), f('x_element_confidence')

    read_its = read_its_lengths(row.get('y_prime_positions', ''))
    cands, note = y_prime_donor_candidates(chr_end, observed, div_idx,
                                           native_arrays, its_vecs, read_its)

    # ---- (A1) recombination_confidence: did recombination happen? ----
    y_conf = y_prime_recomb_confidence(status, len(observed), dsc) if y_sw else 0.0
    ev = []
    if sp_sw: ev.append(sp_conf)
    if x_sw:  ev.append(x_conf)
    if y_sw:  ev.append(y_conf)
    recomb_conf = noisy_or(ev)
    detected = bool(ev)

    # ---- (A2) donor_confidence: weighted evidence, argmax, margin ----
    votes = defaultdict(float)
    if sp_sw and row.get('spacer_source') and row['spacer_source'] != chr_end:
        votes[row['spacer_source']] += W['spacer'] * sp_conf
    if x_sw and row.get('x_element_source') and row['x_element_source'] != chr_end:
        votes[row['x_element_source']] += W['x_element'] * x_conf
    if cands:                       # split Y' weight across tied candidates
        share = W['y_prime'] * y_conf / len(cands)
        for c, s, _ in cands:
            votes[c] += share * s
    m = re.search(r'supplementary=(chr\d+)', row.get('cross_feature_detail', '') or '')
    if m:                           # chr-level: split across that chr's arms
        arms = [e for e in native_arrays if e.startswith(m.group(1))
                and not e[len(m.group(1))].isdigit()]
        for a in (arms or []):
            if a != chr_end:
                votes[a] += W['supplementary'] / len(arms)

    if votes:
        total = sum(votes.values())
        donor, best = max(votes.items(), key=lambda kv: kv[1])
        donor_conf = best / total
    else:
        donor, donor_conf = 'ambiguous', 0.0

    return {
        'recombination_detected_v2': detected,
        'recombination_confidence': round(recomb_conf, 4),
        'donor': donor,
        'donor_confidence': round(donor_conf, 4),
        'y_prime_gained_array': ','.join(observed[div_idx:]) if div_idx >= 0 else '',
        'y_prime_donor_candidates': ','.join(f'{c}:{s}' for c, s, _ in cands) or '-',
        'y_prime_donor_note': note,
    }


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features-dir', required=True)
    ap.add_argument('--y-prime-lib', required=True)
    ap.add_argument('--bed', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    native = build_native_arrays(a.y_prime_lib)
    its = build_its_vectors(a.bed)
    print('Native Y\' arrays:', {k: v for k, v in sorted(native.items())}, '\n')

    rows = []
    for fp in sorted(glob.glob(os.path.join(a.features_dir, '*_features.tsv'))):
        with open(fp) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                out = dict(row)
                out.update(reconcile_v2(row, native, its))
                rows.append(out)
    if not rows:
        sys.exit('no features.tsv rows found')

    cols = ['read_id', 'chr_end', 'y_prime_recombination_status',
            'y_prime_observed_array', 'y_prime_gained_array',
            'recombination_source', 'overall_confidence',
            'recombination_confidence', 'donor', 'donor_confidence',
            'y_prime_donor_candidates', 'y_prime_donor_note']
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f'wrote {len(rows)} rows -> {a.out}')


if __name__ == '__main__':
    main()
