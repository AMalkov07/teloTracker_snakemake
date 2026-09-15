#!/usr/bin/env python3
"""Test whether a read that matches NO single Y' element well is a recombinant hybrid of two.

Motivation: a read carrying a recombinant Y' -- anchor-proximal half from one element,
telomere-distal half from another -- will by construction fail to match any single library
entry end to end. Scoring it against the library one element at a time therefore produces a
mediocre best hit and looks like a matching failure. The diagnostic is to SPLICE the two
candidate references at the homology block and ask whether the read then matches the
splice product in one clean full-length alignment.

    read vs reference A   -> fragmented, two blocks, mediocre identity
    read vs reference B   -> fragmented, two blocks, mediocre identity
    read vs A:B hybrid    -> ONE full-length block at high identity   <- recombinant

A non-recombinant control read does the opposite: it matches its own reference in one
full-length block and matches the hybrid only in fragments.

Usage:
  test_recombinant_hybrid.py --reads <fasta> --ref-a <fasta> --ref-b <fasta>
                             --splice <int>  [--controls <fasta>] [--out <tsv>]

The splice takes TWO coordinates, not one: ref-a[:splice_a] + ref-b[splice_b:]. Y' elements
differ in length, so the point in B equivalent to position splice_a in A is at a different
coordinate -- align A against B first and read off the corresponding position, or the hybrid
will be nonsense and the test will silently report "single-element" for everything.
"""
import argparse, os, subprocess, sys, tempfile


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name: seqs[name] = ''.join(buf)
            name = line[1:].split()[0].strip(); buf = []
        else:
            buf.append(line.strip())
    if name: seqs[name] = ''.join(buf)
    return seqs


def best_hit(query_fa, subject_fa, min_len=1000):
    """Best HSP per query against a single subject: (pident, aln_len, bitscore)."""
    with tempfile.TemporaryDirectory() as td:
        db = os.path.join(td, 'db')
        subprocess.run(['makeblastdb', '-in', subject_fa, '-dbtype', 'nucl', '-out', db],
                       check=True, capture_output=True)
        out = subprocess.run(
            ['blastn', '-query', query_fa, '-db', db, '-evalue', '1e-10',
             '-outfmt', '6 qseqid pident length bitscore'],
            check=True, capture_output=True, text=True).stdout
    best = {}
    for line in out.splitlines():
        rid, pid, ln, bs = line.split('\t')
        ln, bs, pid = int(ln), float(bs), float(pid)
        if ln < min_len: continue
        if rid not in best or bs > best[rid][2]:
            best[rid] = (pid, ln, bs)
    return best


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--reads', required=True)
    p.add_argument('--ref-a', required=True)
    p.add_argument('--ref-b', required=True)
    p.add_argument('--splice-a', type=int, required=True,
                   help='keep ref-a[:splice_a] (1-based inclusive end of the A-derived part)')
    p.add_argument('--splice-b', type=int, required=True,
                   help='then append ref-b[splice_b:] -- the EQUIVALENT point in B, which is '
                        'not the same coordinate as splice_a when A and B differ in length')
    p.add_argument('--controls')
    p.add_argument('--out')
    a = p.parse_args()

    sa = list(read_fasta(a.ref_a).values())[0]
    sb = list(read_fasta(a.ref_b).values())[0]
    name_a = list(read_fasta(a.ref_a).keys())[0]
    name_b = list(read_fasta(a.ref_b).keys())[0]

    with tempfile.TemporaryDirectory() as td:
        hyb = os.path.join(td, 'hybrid.fasta')
        open(hyb, 'w').write(f'>hybrid\n{sa[:a.splice_a]}{sb[a.splice_b:]}\n')

        rows = []
        for label, reads_fa in (('read', a.reads), ('control', a.controls)):
            if not reads_fa: continue
            ba = best_hit(reads_fa, a.ref_a)
            bb = best_hit(reads_fa, a.ref_b)
            bh = best_hit(reads_fa, hyb)
            for rid in sorted(set(ba) | set(bb) | set(bh)):
                g = lambda d: d.get(rid, (0.0, 0, 0.0))
                verdict = 'RECOMBINANT' if g(bh)[2] > max(g(ba)[2], g(bb)[2]) * 1.15 else 'single-element'
                rows.append((label, rid, *g(ba), *g(bb), *g(bh), verdict))

    hdr = ['set', 'read', f'{name_a}_pid', f'{name_a}_len', f'{name_a}_bits',
           f'{name_b}_pid', f'{name_b}_len', f'{name_b}_bits',
           'hybrid_pid', 'hybrid_len', 'hybrid_bits', 'verdict']
    print('\t'.join(hdr))
    for r in rows:
        print('\t'.join(str(x) for x in r))
    if a.out:
        with open(a.out, 'w') as fh:
            fh.write('\t'.join(hdr) + '\n')
            for r in rows:
                fh.write('\t'.join(str(x) for x in r) + '\n')
        print(f'\nwritten {a.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
