#!/usr/bin/env python3
"""For each recipient/donor Y' pair, map the homology available for a crossover.

A recombination junction needs a stretch where donor and recipient are similar enough for
strand invasion. This aligns the two reference elements and reports every high-identity
block between them, so a pair with no usable homology stands out.

Usage: pair_homology.py --pairs <tsv: expected_element donor_element ...> --lib <elements.fasta>
                        [--min-identity 98] [--min-block 100] [--out <tsv>]
"""
import argparse, csv, os, re, subprocess, sys, tempfile


def lib_seqs(path):
    out, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name: out[name] = ''.join(buf)
            m = re.search(r'/E-(chr\w+?[LR])-(\d+)$', line.strip())
            name = f'E-{m.group(1)}-{m.group(2)}' if m else line[1:].split()[0]
            buf = []
        else: buf.append(line.strip())
    if name: out[name] = ''.join(buf)
    return out


def hsps(a_seq, b_seq):
    with tempfile.TemporaryDirectory() as td:
        q = os.path.join(td, 'q.fa'); open(q, 'w').write(f'>a\n{a_seq}\n')
        r = os.path.join(td, 'r.fa'); open(r, 'w').write(f'>b\n{b_seq}\n')
        db = os.path.join(td, 'db')
        subprocess.run(['makeblastdb', '-in', r, '-dbtype', 'nucl', '-out', db],
                       check=True, capture_output=True)
        out = subprocess.run(['blastn', '-query', q, '-db', db, '-evalue', '1e-5',
                              '-outfmt', '6 pident length qstart qend sstart send bitscore'],
                             check=True, capture_output=True, text=True).stdout
    res = []
    for line in out.splitlines():
        pid, ln, qs, qe, ss, se, bs = line.split('\t')
        res.append(dict(pident=float(pid), length=int(ln), qstart=int(qs), qend=int(qe),
                        sstart=int(ss), send=int(se), bits=float(bs)))
    return res


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--pairs', required=True); p.add_argument('--lib', required=True)
    p.add_argument('--min-identity', type=float, default=98.0)
    p.add_argument('--min-block', type=int, default=100)
    p.add_argument('--out')
    a = p.parse_args()
    lib = lib_seqs(a.lib)

    seen, rows = set(), []
    for c in csv.DictReader(open(a.pairs), delimiter='\t'):
        e = c.get('expected_element') or c.get('true_element')
        d = c.get('donor_element') or c.get('matched_element')
        e = e if e.startswith('E-') else 'E-' + e
        d = d if d.startswith('E-') else 'E-' + d
        if (e, d) in seen: continue
        seen.add((e, d))
        if e not in lib or d not in lib: continue
        hs = hsps(lib[e], lib[d])
        good = [h for h in hs if h['pident'] >= a.min_identity and h['length'] >= a.min_block]
        good.sort(key=lambda h: -h['length'])
        total = sum(h['length'] for h in good)
        best = good[0] if good else None
        rows.append({
            'expected_element': e.replace('E-', ''), 'donor_element': d.replace('E-', ''),
            'expected_len': len(lib[e]), 'donor_len': len(lib[d]),
            'n_blocks': len(good), 'total_homologous_bp': total,
            'longest_block_bp': best['length'] if best else 0,
            'longest_block_identity': f"{best['pident']:.2f}%" if best else '-',
            'longest_block_in_expected': f"{best['qstart']}-{best['qend']}" if best else '-',
            'longest_block_in_donor': f"{min(best['sstart'],best['send'])}-{max(best['sstart'],best['send'])}" if best else '-',
            'usable_for_crossover': 'yes' if best and best['length'] >= 200 else
                                    ('marginal' if best else 'NO'),
        })
    cols = list(rows[0].keys())
    if a.out:
        with open(a.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
            for r in rows: w.writerow(r)
        print(f'written {a.out}', file=sys.stderr)
    print('\t'.join(cols))
    for r in rows: print('\t'.join(str(r[c]) for c in cols))


if __name__ == '__main__':
    main()
