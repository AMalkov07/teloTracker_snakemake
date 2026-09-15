#!/usr/bin/env python3
"""Decide whether a read that matches no single Y' element well is a RECOMBINANT.

A recombinant Y' -- anchor-proximal part from one element, telomere-distal part from
another -- matches no single library entry end to end, so per-element scoring returns a
mediocre best hit and looks like a matching failure. The signature is instead positional:
slide a window along the read's Y' region and score each window against both candidate
elements. A recombinant shows a CROSSOVER -- one element wins consistently at one end of
the read, the other wins consistently at the other end. A read that is simply noisy, or
whose reference is wrong, shows one element winning throughout (or neither).

Classification per read:
  RECOMBINANT   both elements win a run of >= min_run decisive windows, on opposite sides,
                and the two runs do not interleave
  single        one element wins all decisive windows
  unresolved    too few decisive windows (the two elements are too similar to tell apart)

"Decisive" means the two identities differ by more than --margin (default 1.0 %).

Usage: scan_recombinant_junctions.py --reads <fasta, headers "sample|read_id">
                                     --lib-dir <dir of elem_<sample>.fasta>
                                     --candidates <tsv: sample read_id true_element matched_element ...>
                                     [--window 300] [--step 150] [--margin 1.0] [--min-run 2]
                                     [--out <tsv>]
"""
import argparse, csv, os, re, subprocess, sys, tempfile
from collections import defaultdict


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
    """element id -> sequence, from a make_element_yprime_lib.py library"""
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


def window_identities(read_seq, ref_seq, window, step):
    """[(read_pos, pident or 0.0)] for each window of the read against one reference."""
    with tempfile.TemporaryDirectory() as td:
        qs = os.path.join(td, 'q.fasta')
        positions = list(range(0, max(1, len(read_seq) - window), step))
        with open(qs, 'w') as fh:
            for p in positions:
                fh.write(f'>w{p}\n{read_seq[p:p+window]}\n')
        db = os.path.join(td, 'db')
        open(db + '.fasta', 'w').write(f'>ref\n{ref_seq}\n')
        subprocess.run(['makeblastdb', '-in', db + '.fasta', '-dbtype', 'nucl', '-out', db],
                       check=True, capture_output=True)
        out = subprocess.run(['blastn', '-query', qs, '-db', db, '-evalue', '1e-5',
                              '-outfmt', '6 qseqid pident length bitscore'],
                             check=True, capture_output=True, text=True).stdout
    best = {}
    for line in out.splitlines():
        q, pid, ln, bs = line.split('\t')
        if int(ln) < window * 0.8: continue
        pid, bs = float(pid), float(bs)
        if q not in best or bs > best[q][1]: best[q] = (pid, bs)
    return [(p, best.get(f'w{p}', (0.0, 0.0))[0]) for p in positions]


def classify(wa, wb, margin, min_run):
    """wa/wb: [(pos, pident)] for the true and matched element. Returns (verdict, detail)."""
    calls = []
    for (p, a), (_, b) in zip(wa, wb):
        if a == 0.0 and b == 0.0: calls.append((p, '.'))
        elif a > b + margin: calls.append((p, 'A'))
        elif b > a + margin: calls.append((p, 'B'))
        else: calls.append((p, '-'))
    decisive = [c for c in calls if c[1] in 'AB']
    if len(decisive) < min_run * 2:
        return 'unresolved', ''.join(c[1] for c in calls)
    seq = [c[1] for c in decisive]
    runs = []
    for ch in seq:
        if runs and runs[-1][0] == ch: runs[-1][1] += 1
        else: runs.append([ch, 1])
    big = [r for r in runs if r[1] >= min_run]
    kinds = {r[0] for r in big}
    verdict = 'single'
    if len(kinds) == 2 and len([r for r in big if r[1] >= min_run]) <= 3:
        verdict = 'RECOMBINANT'
    elif len(kinds) == 2:
        verdict = 'mixed'
    return verdict, ''.join(c[1] for c in calls)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--reads', required=True)
    p.add_argument('--lib-dir', required=True)
    p.add_argument('--candidates', required=True)
    p.add_argument('--window', type=int, default=300)
    p.add_argument('--step', type=int, default=150)
    p.add_argument('--margin', type=float, default=1.0)
    p.add_argument('--min-run', type=int, default=2)
    p.add_argument('--out')
    a = p.parse_args()

    reads = read_fasta(a.reads)
    cands = list(csv.DictReader(open(a.candidates), delimiter='\t'))
    # one row per (sample, read); if a read has several bad copies, take the first
    seen, todo = set(), []
    for c in cands:
        k = (c['sample'], c['read_id'])
        if k in seen: continue
        seen.add(k); todo.append(c)

    libs = {}
    rows = []
    for i, c in enumerate(todo, 1):
        s, rid = c['sample'], c['read_id']
        key = f'{s}|{rid}'
        if key not in reads: continue
        if s not in libs: libs[s] = lib_seqs(os.path.join(a.lib_dir, f'elem_{s}.fasta'))
        lib = libs[s]
        ta, tb = c['true_element'], c['matched_element']
        if ta not in lib or tb not in lib: continue
        try:
            ys, ye = int(c['yp_start']), int(c['yp_end'])
        except (ValueError, KeyError):
            ys, ye = 0, len(reads[key])
        sub = reads[key][max(0, ys - 200):ye + 200]
        wa = window_identities(sub, lib[ta], a.window, a.step)
        wb = window_identities(sub, lib[tb], a.window, a.step)
        verdict, pattern = classify(wa, wb, a.margin, a.min_run)
        rows.append({'sample': s, 'read_id': rid, 'chr_end': c['chr_end'],
                     'true_element': ta, 'matched_element': tb,
                     'verdict': verdict, 'pattern': pattern})
        if i % 25 == 0: print(f'  ...{i}/{len(todo)}', file=sys.stderr)

    hdr = ['sample', 'read_id', 'chr_end', 'true_element', 'matched_element', 'verdict', 'pattern']
    if a.out:
        with open(a.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=hdr, delimiter='\t'); w.writeheader()
            for r in rows: w.writerow(r)
    from collections import Counter
    print('\nverdicts:')
    for k, v in Counter(r['verdict'] for r in rows).most_common():
        print(f'  {k:<14} {v}')
    print('\n(pattern: A = true element wins the window, B = matched element wins, '
          '- = too close to call, . = neither matches)')
    if a.out: print(f'\nwritten {a.out}')


if __name__ == '__main__':
    main()
