#!/usr/bin/env python3
"""ITS tie-breaking for Y' donor prediction.

When the Y' cluster array ties across several candidate donor ends, compare the
internal telomeric tract (ITS) lying BETWEEN the gained Y' elements against each
candidate's reference tract, and break the tie only when the evidence is decisive.

WHY LENGTH AND NOT SEQUENCE
---------------------------
Measured on 116 non-recombinant ("No Change") tracts against their own reference:

  edit distance : p50=3  p90=6   p95=12   (97.4% of edits are indels,
                                           77.5% of those in homopolymer context)
  length diff   : p50=1  p90=5   p95=9

Substitutions are only 2.6% of edits and are non-recurrent across reads (15 of 17
sites seen in a single read), i.e. basecalling noise, not real variants. So there
are no SNVs to catch, and sequence comparison mostly measures read quality: an
insertion and a deletion cancel in a length comparison but both count as edits.
Length is therefore the tighter, more honest statistic here.

Run-length collapsing was evaluated and rejected: the reference tracts contain
zero homopolymer runs >4 (yeast TG(1-3)), so a max-4 cap is a no-op, and full
RLE only improved mean edit distance 6.4 -> 4.5 while adding a parameter.

GATES (a donor is named only if all three hold)
----------------------------------------------
  1. best candidate tract within TOL_BP of the read's tract
  2. the two candidate tracts differ from each other by >= SEP_FRAC
     (are these ends distinguishable in principle? -- read-independent)
  3. the runner-up is worse by >= max(SEP_MIN_BP, SEP_FRAC)
     (did THIS read actually discriminate? -- read-dependent)

max() not min() in gate 3: min() would collapse the requirement to <1bp on short
tracts, admitting 1bp "discriminations" that are pure noise.

Never expands the candidate set -- it can only narrow it. Calls are capped at
CEILING confidence and tagged donor_evidence=its so they stay filterable.

Yield on 7172 day-6: 12 of 78 eligible reads (62 rejected because the candidate
tracts are inherently indistinguishable, 4 because the read matched no tract).
All 12 rest on a ~9bp vs ~169bp contrast, so this is a narrow rule, not a general
donor predictor -- check a strain's tract-length spectrum before relying on it.
"""
import argparse, collections, csv, glob, re

TOL_BP     = 10      # best candidate tract must be within +/- this of the read
SEP_FRAC   = 0.10    # candidates must differ by >=10%; runner-up >=10% worse
SEP_MIN_BP = 10      # ...and never less than this many bp (hard safety floor)
CEILING    = 0.75    # ITS-only evidence never looks corroborated


def read_fasta(p):
    n, buf = None, []
    for line in open(p):
        if line.startswith('>'):
            if n: yield n, ''.join(buf)
            n, buf = line[1:].split()[0], []
        else: buf.append(line.strip())
    if n: yield n, ''.join(buf)


def load_reference_its(bed, ref_fasta):
    """chr_end -> ordered list of inter-Y' tract sequences (anchor->telomere)."""
    contigs = dict(read_fasta(ref_fasta))
    def resolve(name):
        return contigs.get(name) or contigs.get(name + '_extended')
    feats = collections.defaultdict(dict)
    for line in open(bed):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        m = re.match(r'ITS_(chr\d+[LR])_Y_Prime_(\d+)-(\d+)$', f[3])
        if not m: continue
        seq = resolve(f[0])
        if seq is None: continue
        feats[m.group(1)][(int(m.group(2)), int(m.group(3)))] = seq[int(f[1]):int(f[2])]
    return {e: [d[k] for k in sorted(d)] for e, d in feats.items()}


def build_native(lib):
    """chr_end -> ordered list of Y' cluster IDs."""
    pos = collections.defaultdict(dict)
    for line in open(lib):
        if not line.startswith('>'): continue
        h = line[1:].strip()
        nm, c = (h.split('#', 1) + [''])[:2]
        q = c.split('/')
        yid = q[2].split('_', 1)[0] if len(q) >= 3 else ''
        if not yid: continue
        for tok in re.sub(r'^Y_Prime_', '', nm.split()[0]).split(';'):
            m = re.match(r'(chr\d+[LR])([\d,]+)$', tok)
            if m:
                for x in m.group(2).split(','):
                    if x: pos[m.group(1)][int(x)] = yid
    return {e: [d[k] for k in sorted(d)] for e, d in pos.items()}


def full_match_starts(gained, native):
    """Start offsets where `gained` matches `native` in full."""
    out, g, n = [], len(gained), len(native)
    for s in range(n):
        k = 0
        while s + k < n and k < g and native[s + k] == gained[k]: k += 1
        if k == g: out.append(s)
    return out


def its_tiebreak(read_tract_bp, candidates, ref_its):
    """candidates: [(chr_end, match_start)]. -> (donor, confidence, note)."""
    scores = {}
    for end, start in candidates:
        v = ref_its.get(end, [])
        j = start + 1                       # tract between gained[0] and gained[1]
        if j < len(v) and v[j]:
            diff = abs(read_tract_bp - len(v[j]))
            if end not in scores or diff < scores[end][0]:
                scores[end] = (diff, len(v[j]))
    if len(scores) < 2:
        return None, 0.0, 'insufficient_its_data'

    rank = sorted(scores.items(), key=lambda kv: kv[1][0])
    (best_end, (bd, bl)), (_, (sd, sl)) = rank[0], rank[1]

    if bd > TOL_BP:
        return None, 0.0, f'no_tract_match(best={bd}bp)'
    if abs(bl - sl) < SEP_FRAC * max(bl, sl):
        return None, 0.0, 'candidates_indistinguishable'
    need = max(SEP_MIN_BP, SEP_FRAC * max(bl, sl))
    if (sd - bd) < need:
        return None, 0.0, f'read_did_not_discriminate({sd-bd}<{need:.0f})'

    quality = 1.0 - bd / (TOL_BP + 1)             # how well the read fits the winner
    margin  = min(1.0, (sd - bd) / max(1.0, need))  # how decisively it beat the rest
    return best_end, round(min(CEILING, quality * margin * CEILING), 3), 'its_resolved'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features-dir', required=True)
    ap.add_argument('--bed', required=True)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--lib', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    ref_its = load_reference_its(a.bed, a.reference)
    native  = build_native(a.lib)
    rows = []
    for fp in sorted(glob.glob(f'{a.features_dir}/*_features.tsv')):
        rows.extend(csv.DictReader(open(fp), delimiter='\t'))

    out, notes, eligible = [], collections.Counter(), 0
    for r in rows:
        obs = [x for x in (r['y_prime_observed_array'] or '').split(',') if x]
        try: d = int(r['y_prime_divergence_idx'])
        except ValueError: d = -1
        if d < 0 or not obs: continue
        gained = obs[d:]
        if len(gained) < 2:                 # need >=2 gained Y' for an internal tract
            continue
        cands = []
        for c, nat in native.items():
            if c == r['chr_end'] or not nat: continue
            for s in full_match_starts(gained, nat): cands.append((c, s))
        ends = sorted(set(c for c, _ in cands))
        if len(ends) < 2: continue
        eligible += 1

        spans = []
        for tok in (r['y_prime_positions'] or '').split(';'):
            m = re.match(r'[^:]+:(\d+)-(\d+)$', tok.strip())
            if m: spans.append(tuple(sorted((int(m.group(1)), int(m.group(2))))))
        if len(spans) <= d + 1:
            notes['no_spans'] += 1; continue
        x, y = spans[d], spans[d + 1]
        tract_bp = (y[0] - x[1]) if y[0] >= x[1] else (x[0] - y[1])
        if tract_bp < 0:
            notes['overlapping_yprimes'] += 1; continue

        donor, conf, note = its_tiebreak(tract_bp, cands, ref_its)
        notes[re.sub(r'\(.*', '', note)] += 1
        if donor:
            out.append((r['read_id'], r['chr_end'], ','.join(gained), ','.join(ends),
                        donor, tract_bp, conf, 'its'))

    print(f'eligible (>=2 gained Y\', multi-end tie): {eligible}')
    print(f'resolved by ITS: {len(out)}\n')
    for k, v in notes.most_common(): print(f'  {v:5d}  {k}')
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['read_id','chr_end','gained','cluster_candidates',
                    'its_donor','read_tract_bp','donor_confidence','donor_evidence'])
        w.writerows(out)
    print(f'\nwrote {a.out}')
    for o in out:
        print(f'  {o[1]:7} gained={o[2][:16]:<16} {len(o[3].split(",")):2d} cands -> '
              f'{o[4]:<7} tract={o[5]:4d}bp conf={o[6]}')


if __name__ == '__main__':
    main()
