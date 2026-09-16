#!/usr/bin/env python3
import sys, csv, json, edlib

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(COMP)[::-1]

def load_fasta(f):
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
    L = len(oriented)
    if L < 2 * min_arm + step: return None
    best = None
    for J in range(min_arm, L - min_arm, step):
        jo = int(J * len(own_ref) / L)
        jd = int(J * len(donor_ref) / L)
        own_first = nw_identity(oriented[:J], own_ref[:jo])
        don_first = nw_identity(oriented[:J], donor_ref[:jd])
        own_second = nw_identity(oriented[J:], own_ref[jo:])
        don_second = nw_identity(oriented[J:], donor_ref[jd:])
        score1 = (own_first - don_first) + (don_second - own_second)
        score2 = (don_first - own_first) + (own_second - don_second)
        for score, direction in ((score1, 'own->donor'), (score2, 'donor->own')):
            if best is None or score > best[0]:
                best = (score, direction, J, own_first, don_first, own_second, don_second)
    return best

tag = sys.argv[1]
reads = load_fasta(f'reads_{tag}.fasta')
elems = json.load(open(f'elems_{tag}.json'))
master = list(csv.DictReader(open(f'master_{tag}.tsv'), delimiter='\t'))
targets = [r for r in master if "partial junction" in r['verdict']]

print(f'{"chr_end":<8}{"read_id":<22}{"1st own":>8}{"1st don":>8}{"marg":>7}   {"2nd own":>8}{"2nd don":>8}{"marg":>7}  clean?')
rows_out = []
for r in targets:
    rid = r['read_id']
    own_elem, donor_elem = r['own_elem'], r['donor_elem']
    seq = reads[rid]
    ys, ye = int(r['yp_start']), int(r['yp_end'])
    sub = seq[max(0, ys - 200):ye + 200]
    oriented = best_strand_sub(sub, elems[own_elem])
    best = scan_best_split(oriented, elems[own_elem], elems[donor_elem])
    if best is None:
        print(f'{r["chr_end"]:<8}{rid:<22} -- too short to split --')
        continue
    score, direction, J, of, df, os_, ds = best
    fm, sm = df - of, ds - os_
    clean = 'YES' if (fm > 5 and sm < -5) or (fm < -5 and sm > 5) else 'no'
    print(f"{r['chr_end']:<8}{rid:<22}{of:>7.1f}%{df:>7.1f}%{fm:>+7.1f}   {os_:>7.1f}%{ds:>7.1f}%{sm:>+7.1f}  {clean}")
    rows_out.append(dict(chr_end=r['chr_end'], read_id=rid, own_elem=own_elem, donor_elem=donor_elem,
                          split_bp=J, first_own=round(of,2), first_donor=round(df,2),
                          second_own=round(os_,2), second_donor=round(ds,2),
                          first_margin=round(fm,1), second_margin=round(sm,1), clean_crossover=clean))

if rows_out:
    with open(f'half_split_{tag}.tsv', 'w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=list(rows_out[0].keys()), delimiter='\t')
        w.writeheader()
        for row in rows_out: w.writerow(row)
