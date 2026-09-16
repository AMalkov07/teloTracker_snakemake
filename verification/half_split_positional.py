#!/usr/bin/env python3
"""Position-enforced half-split identity test for mid-Y' partial-junction reads.

Fixes a real gap in flank_identities.py: that script BLASTed each read-half against the
WHOLE reference and took the best-scoring HSP anywhere in it -- so a spuriously good match
to an unrelated internal region (e.g. the shared tandem-repeat region already characterised
in this investigation) could pass as "evidence" even though it says nothing about a real
junction. Here, the read's own-vs-donor best split point J is searched for directly, and at
every candidate J the read's first part is ONLY ever compared to the reference's first
proportional part, and the second part ONLY to the reference's second part -- via edlib
global (NW) alignment of the two same-region slices, so there is no freedom for either half
to "find" a match somewhere else in the reference.
"""
import sys, csv, edlib

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(COMP)[::-1]

def load(f):
    d = {}; n = None; b = []
    for l in open(f):
        if l.startswith('>'):
            if n: d[n] = ''.join(b)
            n = l[1:].split()[0]; b = []
        else: b.append(l.strip())
    if n: d[n] = ''.join(b)
    return d

def nw_identity(a, b):
    if not a or not b: return 0.0
    res = edlib.align(a, b, mode="NW", task="distance")
    return 100 * (1 - res['editDistance'] / max(len(a), len(b)))

def best_strand_sub(sub, own_ref):
    best = None
    for r in (sub, rc(sub)):
        res = edlib.align(own_ref, r, mode="HW", task="distance")
        ident = 100 * (1 - res['editDistance'] / len(own_ref))
        if best is None or ident > best[0]:
            best = (ident, r)
    return best[1]

def scan_best_split(oriented, own_ref, donor_ref, min_arm=600, step=100):
    """Try candidate split points J in `oriented`; at each, compare
    oriented[:J] to the proportional first part of each ref, and oriented[J:] to the
    proportional second part. Return the J and direction that best explains a real
    crossover (one ref wins first part, the OTHER wins second part)."""
    L = len(oriented)
    best = None
    for J in range(min_arm, L - min_arm, step):
        jo = int(J * len(own_ref) / L)
        jd = int(J * len(donor_ref) / L)
        own_first = nw_identity(oriented[:J], own_ref[:jo])
        don_first = nw_identity(oriented[:J], donor_ref[:jd])
        own_second = nw_identity(oriented[J:], own_ref[jo:])
        don_second = nw_identity(oriented[J:], donor_ref[jd:])
        # hypothesis 1: own wins first half, donor wins second half
        score1 = (own_first - don_first) + (don_second - own_second)
        # hypothesis 2: donor wins first half, own wins second half
        score2 = (don_first - own_first) + (own_second - don_second)
        for score, direction in ((score1, 'own->donor'), (score2, 'donor->own')):
            if best is None or score > best[0]:
                best = (score, direction, J, own_first, don_first, own_second, don_second)
    return best

reads = load('all_mismatch_reads.fasta')
refs = load('all_elements_58.fasta')

master = list(csv.DictReader(open('all58/master_results.tsv'), delimiter='\t'))
mismatches = {r['read_id']: r for r in csv.DictReader(open('all58/all_mismatches.tsv'), delimiter='\t')}

targets = [r for r in master if "partial junction" in r['verdict']]

print(f'{"chr_end":<8}{"read_id":<22}{"direction":<12}{"split":>7}  {"1st:own":>9}{"1st:donor":>11}  {"2nd:own":>9}{"2nd:donor":>11}  contrast')
rows_out = []
for r in targets:
    rid = r['read_id']
    m = mismatches[rid]
    own_elem = m['true_element'].replace('E-', '')
    donor_elem = r['donor_elem']
    seq = reads[rid]
    ys, ye = int(m['yp_start']), int(m['yp_end'])
    sub = seq[max(0, ys - 200):ye + 200]
    oriented = best_strand_sub(sub, refs[own_elem])
    best = scan_best_split(oriented, refs[own_elem], refs[donor_elem])
    if best is None:
        print(f'{r["chr_end"]:<8}{rid:<22} -- too short to split --')
        continue
    score, direction, J, of, df, os_, ds = best
    contrast = f'{score:+.1f}'
    print(f'{r["chr_end"]:<8}{rid:<22}{direction:<12}{J:>6}bp  {of:>8.2f}%{df:>10.2f}%  {os_:>8.2f}%{ds:>10.2f}%  {contrast}')
    rows_out.append(dict(chr_end=r['chr_end'], read_id=rid, own_elem=own_elem, donor_elem=donor_elem,
                          direction=direction, split_bp=J, sub_len=len(oriented),
                          first_own=round(of,2), first_donor=round(df,2),
                          second_own=round(os_,2), second_donor=round(ds,2), contrast=round(score,2)))

with open('all58/half_split_positional.tsv', 'w', newline='') as out:
    w = csv.DictWriter(out, fieldnames=list(rows_out[0].keys()), delimiter='\t')
    w.writeheader()
    for row in rows_out: w.writerow(row)
