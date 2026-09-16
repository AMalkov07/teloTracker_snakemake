#!/usr/bin/env python3
"""Global-alignment identity test (own vs full donor group) for 7172/7302 mismatched reads."""
import sys, csv, json, edlib
from collections import defaultdict

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

def infix_identity(read, ref):
    best = None
    for r in (read, rc(read)):
        res = edlib.align(ref, r, mode="HW", task="distance")
        ident = 100 * (1 - res['editDistance'] / len(ref))
        if best is None or ident > best:
            best = ident
    return best

tag = sys.argv[1]
reads = load_fasta(f'reads_{tag}.fasta')
elems = json.load(open(f'elems_{tag}.json'))
G = json.load(open(f'groups_{tag}.json'))['groups']
group_members = defaultdict(list)
for e, g in G.items():
    group_members[g].append(e.replace('E-', ''))

rows = []
with open(f'mismatches_{tag}.tsv') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        rid = row['read_id']
        own_elem = row['true_element'].replace('E-', '')
        donor_group = row['matched_group']
        seq = reads[rid]
        own_id = infix_identity(seq, elems[own_elem])
        best_donor = max(((infix_identity(seq, elems[m]), m) for m in group_members[donor_group]), key=lambda x: x[0])
        don_id, don_name = best_donor
        gap = don_id - own_id
        if gap > 5: verdict = 'RECOMBINANT (whole element)'
        elif gap > -2: verdict = "mid-Y' partial junction (unconfirmed breakpoint)"
        else: verdict = 'not supported'
        rows.append((row['chr_end'], rid, own_elem, own_id, donor_group, don_name, don_id, gap, verdict, row['yp_start'], row['yp_end']))

rows.sort(key=lambda x: (x[0], -x[7]))
print(f'{"chr_end":<8}{"read_id":<22}{"own_elem":<11}{"own%":>7}{"donor_grp":<10}{"donor_elem":<12}{"donor%":>8}{"gap":>8}  verdict')
for r in rows:
    chr_end, rid, own_elem, own_id, dg, dn, did, gap, verdict, ys, ye = r
    print(f'{chr_end:<8}{rid:<22}{own_elem:<11}{own_id:>6.2f}%{dg:<10}{dn:<12}{did:>7.2f}%{gap:>+7.2f}%  {verdict}')

with open(f'master_{tag}.tsv', 'w') as out:
    out.write('chr_end\tread_id\town_elem\town_pct\tdonor_group\tdonor_elem\tdonor_pct\tgap\tverdict\typ_start\typ_end\n')
    for r in rows:
        chr_end, rid, own_elem, own_id, dg, dn, did, gap, verdict, ys, ye = r
        out.write(f'{chr_end}\t{rid}\t{own_elem}\t{own_id:.2f}\t{dg}\t{dn}\t{did:.2f}\t{gap:+.2f}\t{verdict}\t{ys}\t{ye}\n')
