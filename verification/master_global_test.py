#!/usr/bin/env python3
"""Run the global-alignment identity test (own ref vs donor GROUP) for every
Y'-mismatch read across the 5 identically-grouped 6991 samples, and print one
master table with a recombination verdict per read."""
import sys, edlib, csv

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

def infix_identity(read, ref):
    best = None
    for r in (read, rc(read)):
        res = edlib.align(ref, r, mode="HW", task="distance")
        ident = 100 * (1 - res['editDistance'] / len(ref))
        if best is None or ident > best:
            best = ident
    return best

GROUP_MEMBERS = {
    'G1': ['chr2L-1', 'chr6L-1'],
    'G2': ['chr13L-1', 'chr14L-3', 'chr14L-4', 'chr14L-5'],
    'G3': ['chr8R-1'],
    'G4': ['chr12L-1'],
    'G5': ['chr8L-1'],
    'G6': ['chr16R-1'],
    'G7': ['chr10L-1', 'chr9L-1'],
    'G8': ['chr12R-2', 'chr12R-3', 'chr12R-4', 'chr12R-5', 'chr12R-6', 'chr12R-7',
           'chr14L-1', 'chr14L-2', 'chr15R-1', 'chr16L-1',
           'chr4R-1', 'chr4R-2', 'chr4R-3', 'chr4R-4', 'chr4R-5', 'chr4R-6', 'chr4R-7', 'chr7R-1'],
    'G9': ['chr5R-1'],
    'G10': ['chr14R-1'],
    'G11': ['chr12R-1'],
    'G12': ['chr5L-1'],
}

reads = load('all_mismatch_reads.fasta')
refs = load('all_elements_58.fasta')

ws_class = {}
with open('/home/andrey/teloTracker_snakemake/verification/reports/cut99_summary_6991.tsv') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['read_id'] not in ws_class:
            ws_class[row['read_id']] = row.get('evidence', '')

rows = []
with open('all58/all_mismatches.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        rid = row['read_id']
        chr_end = row['chr_end']
        own_elem = row['true_element'].replace('E-', '')
        donor_group = row['matched_group']
        seq = reads[rid]
        own_id = infix_identity(seq, refs[own_elem])
        best_donor = max(((infix_identity(seq, refs[m]), m) for m in GROUP_MEMBERS[donor_group]), key=lambda x: x[0])
        don_id, don_name = best_donor
        gap = don_id - own_id
        if gap > 5:
            verdict = 'RECOMBINANT (whole element)'
        elif gap > -2:
            verdict = 'likely recombinant (mid-Y\' partial junction)'
        else:
            verdict = 'not supported'
        wsc = ws_class.get(rid, '')
        rows.append((chr_end, rid, row['sample'], own_elem, own_id, donor_group, don_name, don_id, gap, verdict, wsc))

rows.sort(key=lambda x: (x[0], -x[8]))
print(f'{"chr_end":<8}{"read_id":<22}{"own_elem":<11}{"own%":>7}{"donor_grp":<10}{"donor_elem":<11}{"donor%":>8}{"gap":>8}  {"ws_class":<12}verdict')
for chr_end, rid, sample, own_elem, own_id, dg, dn, did, gap, verdict, wsc in rows:
    print(f'{chr_end:<8}{rid:<22}{own_elem:<11}{own_id:>6.2f}%{dg:<10}{dn:<11}{did:>7.2f}%{gap:>+7.2f}%  {wsc:<12}{verdict}')

with open('all58/master_results.tsv', 'w') as out:
    out.write('chr_end\tread_id\tsample\town_elem\town_pct\tdonor_group\tdonor_elem\tdonor_pct\tgap\twindow_scan_class\tverdict\n')
    for chr_end, rid, sample, own_elem, own_id, dg, dn, did, gap, verdict, wsc in rows:
        out.write(f'{chr_end}\t{rid}\t{sample}\t{own_elem}\t{own_id:.2f}\t{dg}\t{dn}\t{did:.2f}\t{gap:+.2f}\t{wsc}\t{verdict}\n')
