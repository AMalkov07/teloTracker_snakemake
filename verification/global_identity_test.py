#!/usr/bin/env python3
"""Global-alignment identity test for whole-Y'-element recombination.

For each read: compute infix global identity (edlib mode=HW -- the reference aligns fully
end-to-end, read is free at both ends to allow flanking anchor/telomere/other sequence)
against its own expected reference and against every candidate donor-group member, best of
both strands. This is what a true global alignment of read-vs-reference degenerates into once
the read carries extra flanking sequence beyond the Y' itself.

A read whose ENTIRE Y' was replaced by a donor copy should show donor identity clearly
exceeding own identity, over the reference's FULL length -- not just a fragment. A read with
a genuine mid-Y' junction (part native, part donor) should show NEITHER reference reaching
full-length high identity; those need the sliding-window junction scan, not this test.

Usage: global_identity_test.py <reads.fasta> <own_ref.fasta> <donor_group.fasta> [--label X]
Prints one row per read: own identity, best donor identity + which member, gap, verdict.
"""
import sys, edlib

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
    for strand, r in (('+', read), ('-', rc(read))):
        res = edlib.align(ref, r, mode="HW", task="distance")
        ed = res['editDistance']
        ident = 100 * (1 - ed / len(ref))
        if best is None or ident > best[0]:
            best = (ident, strand)
    return best

reads_fa, own_fa, donor_fa = sys.argv[1:4]
label = sys.argv[sys.argv.index('--label') + 1] if '--label' in sys.argv else ''

reads = load(reads_fa)
own_refs = load(own_fa)
donor_refs = load(donor_fa)

print(f'{"read":<22}{"own identity":>13}{"best donor identity":>21}{"donor member":<14}{"gap":>8}  verdict')
n_recomb = 0
for rid, seq in reads.items():
    own_id, own_strand = max((infix_identity(seq, r) for r in own_refs.values()), key=lambda x: x[0])
    best_donor = max(((infix_identity(seq, r), name) for name, r in donor_refs.items()), key=lambda x: x[0][0])
    (don_id, don_strand), don_name = best_donor
    gap = don_id - own_id
    verdict = 'RECOMBINANT (whole element)' if gap > 5 else ('not supported' if gap < -2 else 'ambiguous -- needs junction scan')
    if gap > 5: n_recomb += 1
    print(f'{rid:<22}{own_id:>12.2f}%{don_id:>20.2f}%  {don_name:<12}{gap:>+7.2f}%  {verdict}')
print(f'\n{label}: {n_recomb}/{len(reads)} whole-element recombinant (gap > 5 points)')
