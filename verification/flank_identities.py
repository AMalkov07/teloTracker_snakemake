#!/usr/bin/env python3
"""For each candidate recombinant read, measure how well each HALF of its Y' matches the
expected element versus the donor element -- the evidence that justifies a mid-Y' junction.

The junction is located from a sliding-window scan (identical to
scan_recombinant_junctions.py): the boundary between the last window won by one element and
the first window won by the other. The Y' region is split there, and each half is aligned to
both references, giving a 2x2:

                        vs expected element   vs donor element
    anchor-side half          should be high        should be lower
    telomere-side half        should be lower       should be high

Halves are named by biological side using the read's telo_side, not by read coordinate, so
reads sequenced in either direction are directly comparable.

Usage: flank_identities.py --reads <fasta, "sample|read_id"> --lib-dir <dir>
                           --table <tsv with read_id, chr_end, true/matched element, telo_side,
                                    yp_start, yp_end> --sample <name> [--out <tsv>]
"""
import argparse, csv, os, re, subprocess, sys, tempfile
from collections import Counter

WINDOW, STEP, MARGIN = 300, 150, 1.0


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name: seqs[name] = ''.join(buf)
            name = line[1:].split()[0].strip(); buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = ''.join(buf)
    return seqs


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


def blast_pident(query_seq, ref_seq):
    """Best HSP of query vs ref -> (pident, aln_len, pct_of_query_aligned)."""
    if len(query_seq) < 100: return (0.0, 0, 0.0)
    with tempfile.TemporaryDirectory() as td:
        q = os.path.join(td, 'q.fasta'); open(q, 'w').write(f'>q\n{query_seq}\n')
        r = os.path.join(td, 'r.fasta'); open(r, 'w').write(f'>r\n{ref_seq}\n')
        db = os.path.join(td, 'db')
        subprocess.run(['makeblastdb', '-in', r, '-dbtype', 'nucl', '-out', db],
                       check=True, capture_output=True)
        out = subprocess.run(['blastn', '-query', q, '-db', db, '-evalue', '1e-5',
                              '-outfmt', '6 pident length bitscore'],
                             check=True, capture_output=True, text=True).stdout
    best = (0.0, 0, 0.0)
    for line in out.splitlines():
        pid, ln, bs = line.split('\t')
        if float(bs) > best[2]: best = (float(pid), int(ln), float(bs))
    return (best[0], best[1], round(100.0 * best[1] / len(query_seq), 1))


def window_calls(sub, ref_a, ref_b):
    """[(offset, 'A'|'B'|'-'|'.')] over the read region."""
    positions = list(range(0, max(1, len(sub) - WINDOW), STEP))
    with tempfile.TemporaryDirectory() as td:
        qs = os.path.join(td, 'q.fasta')
        with open(qs, 'w') as fh:
            for p in positions: fh.write(f'>w{p}\n{sub[p:p+WINDOW]}\n')
        res = {}
        for tag, ref in (('A', ref_a), ('B', ref_b)):
            rp = os.path.join(td, f'{tag}.fasta'); open(rp, 'w').write(f'>r\n{ref}\n')
            db = os.path.join(td, f'db{tag}')
            subprocess.run(['makeblastdb', '-in', rp, '-dbtype', 'nucl', '-out', db],
                           check=True, capture_output=True)
            out = subprocess.run(['blastn', '-query', qs, '-db', db, '-evalue', '1e-5',
                                  '-outfmt', '6 qseqid pident length bitscore'],
                                 check=True, capture_output=True, text=True).stdout
            d = {}
            for line in out.splitlines():
                qi, pid, ln, bs = line.split('\t')
                if int(ln) < WINDOW * 0.8: continue
                if qi not in d or float(bs) > d[qi][1]: d[qi] = (float(pid), float(bs))
            res[tag] = d
    calls = []
    for p in positions:
        a = res['A'].get(f'w{p}', (0.0, 0.0))[0]; b = res['B'].get(f'w{p}', (0.0, 0.0))[0]
        if a == 0.0 and b == 0.0: calls.append((p, '.'))
        elif a > b + MARGIN: calls.append((p, 'A'))
        elif b > a + MARGIN: calls.append((p, 'B'))
        else: calls.append((p, '-'))
    return calls


def junction_offset(calls):
    """Midpoint between the last decisive call of the first kind and the first of the other."""
    dec = [(p, c) for p, c in calls if c in 'AB']
    if len(dec) < 4: return None
    first_kind = Counter(c for _, c in dec[:max(1, len(dec)//3)]).most_common(1)[0][0]
    other = 'B' if first_kind == 'A' else 'A'
    last_first = max((p for p, c in dec if c == first_kind), default=None)
    first_other = min((p for p, c in dec if c == other), default=None)
    if last_first is None or first_other is None: return None
    lo = min(p for p, c in dec if c == first_kind)
    # take the transition: last first-kind BEFORE the first other-kind run that follows it
    after = [p for p, c in dec if c == other and p > lo]
    if not after: return None
    fo = min(after)
    lf = max((p for p, c in dec if c == first_kind and p < fo), default=lo)
    return (lf + fo) // 2 + WINDOW // 2


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--reads', required=True); p.add_argument('--lib-dir', required=True)
    p.add_argument('--table', required=True); p.add_argument('--sample', required=True)
    p.add_argument('--out')
    a = p.parse_args()

    reads = read_fasta(a.reads)
    lib = lib_seqs(os.path.join(a.lib_dir, f'elem_{a.sample}.fasta'))
    rows = []
    for c in csv.DictReader(open(a.table), delimiter='\t'):
        rid = c['read_id']
        key = f'{a.sample}|{rid}'
        if key not in reads: continue
        exp_e = c.get('true_element') or ('E-' + c['expected_element'])
        don_e = c.get('matched_element') or ('E-' + c['observed_element'])
        if exp_e not in lib or don_e not in lib: continue
        ys, ye = int(c['yp_start']), int(c['yp_end'])
        pad = 200
        sub = reads[key][max(0, ys - pad):ye + pad]
        calls = window_calls(sub, lib[exp_e], lib[don_e])
        j = junction_offset(calls)
        if j is None or j <= 200 or j >= len(sub) - 200:
            rows.append({'read_id': rid, 'chr_end': c['chr_end'],
                         'expected_element': exp_e.replace('E-', ''),
                         'donor_element': don_e.replace('E-', ''),
                         'junction_in_read': '', 'anchor_half_bp': '', 'telo_half_bp': '',
                         'anchorHalf_vs_expected': '', 'anchorHalf_vs_donor': '',
                         'teloHalf_vs_expected': '', 'teloHalf_vs_donor': '',
                         'note': 'no clean junction'})
            continue
        left, right = sub[:j], sub[j:]
        telo_first = c['telo_side'] == 'beginning'
        telo_half, anchor_half = (left, right) if telo_first else (right, left)
        ae, ad = blast_pident(anchor_half, lib[exp_e]), blast_pident(anchor_half, lib[don_e])
        te, td = blast_pident(telo_half, lib[exp_e]), blast_pident(telo_half, lib[don_e])
        fmt = lambda x: f'{x[0]:.2f}% ({x[1]}bp)'
        rows.append({'read_id': rid, 'chr_end': c['chr_end'],
                     'expected_element': exp_e.replace('E-', ''),
                     'donor_element': don_e.replace('E-', ''),
                     'junction_in_read': ys - pad + j if ys > pad else j,
                     'anchor_half_bp': len(anchor_half), 'telo_half_bp': len(telo_half),
                     'anchorHalf_vs_expected': fmt(ae), 'anchorHalf_vs_donor': fmt(ad),
                     'teloHalf_vs_expected': fmt(te), 'teloHalf_vs_donor': fmt(td),
                     'note': 'anchor half -> expected, telo half -> donor'
                             if (ae[0] > ad[0] and td[0] > te[0]) else 'pattern not clean'})
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
