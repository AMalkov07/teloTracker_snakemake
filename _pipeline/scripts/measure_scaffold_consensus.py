#!/usr/bin/env python3
"""
measure_scaffold_consensus.py -- how similar are the candidate scaffold reads?

For each chr_end, takes the N candidates around the 75th-percentile position that
select_reads would choose from, extracts the region that would actually be
GRAFTED (the telomere-side portion past the anchor -- what becomes the soft-clip
extension), aligns them all-vs-all, and reports how well they agree.

Purpose: decide, from data rather than guesswork, whether a "pick the read that
agrees with >=K others" rule is viable, and what "agrees" should mean
numerically.

Compares the extension region, not the whole read, because two reads can agree
across the anchor and diverge only past it -- which is exactly what a bad
scaffold looks like, and the extension is the part that lands in the reference.

CAVEAT ON THE `similarity` COLUMN
---------------------------------
similarity = pident * min(1, coverage) using the SINGLE BEST HSP. A read that is
fully colinear with the others but carries one large insertion therefore scores
LOW coverage despite being structurally colinear (measured: 7172 chr4L's
SRR33298432.121462 scores 0.60 coverage while summing its 2 HSPs gives 1.00).
This is acceptable for *scaffold selection*, where an indel-bearing read is
just as unwanted as a rearranged one, but it means this column measures
"structural disagreement", NOT "recombination". Classify the flagged reads
separately if you need to know which is which.
"""
import argparse
import os
import subprocess
import sys
import tempfile

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--post-tsv', required=True,
                   help='<base>_post_y_prime_probe.tsv (selection criteria live here)')
    p.add_argument('--all-matches', required=True,
                   help='all_matches_<base>_blasted_<anchor>.tsv (anchor positions)')
    p.add_argument('--reads-fasta', required=True, help='<base>.fasta (needs .fai)')
    p.add_argument('--n', type=int, default=5, help='candidates around the 75th pct')
    p.add_argument('--min-ext', type=int, default=1000,
                   help='ignore candidates whose extension is shorter than this')
    p.add_argument('--output', required=True, help='per-pair TSV')
    p.add_argument('--threads', type=int, default=4)
    return p.parse_args()


def load_fai(path):
    idx = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 5:
            idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return idx


def fetch(fh, rec):
    length, offset, line_bases, line_width = rec
    fh.seek(offset)
    n_lines = (length + line_bases - 1) // line_bases
    return ''.join(fh.read(n_lines * line_width).decode().split())[:length]


def main():
    a = parse_args()
    fai = a.reads_fasta + '.fai'
    if not os.path.exists(fai):
        sys.exit(f'ERROR: missing {fai}')
    idx = load_fai(fai)

    post = pd.read_csv(a.post_tsv, sep='\t')
    am = pd.read_csv(a.all_matches, sep='\t')
    am = am.sort_values('bitscore', ascending=False).drop_duplicates('read_id')
    anchor = am.set_index('read_id')[
        ['match_start_on_read', 'match_end_on_read', 'wanted_section_of_read', 'pident']
    ].to_dict('index')

    # replicate select_reads' candidate filter
    cand = post[(post['repeat_length'] >= 30) & (post['Adapter_After_Telomere'] == True)]

    rows = []
    summary = []
    with open(a.reads_fasta, 'rb') as fh, tempfile.TemporaryDirectory() as tmp:
        for chr_end in sorted(cand['chr_end'].unique(),
                              key=lambda x: (int(''.join(c for c in str(x) if c.isdigit()) or 99), str(x))):
            sub = cand[cand['chr_end'] == chr_end]
            if sub.empty:
                continue
            mode_y = sub['y_prime_probe_count'].mode()
            if len(mode_y) == 0:
                continue
            sub = sub[sub['y_prime_probe_count'] == mode_y[0]]
            if sub.empty:
                continue
            sub = sub.sort_values('repeat_length')
            n_tot = len(sub)
            q = min(int(n_tot * 0.75), n_tot - 1)          # the index select_reads picks
            lo = max(0, q - a.n // 2)
            picks = sub.iloc[lo:lo + a.n]

            recs = []
            for rid in picks['read_id']:
                if rid not in idx or rid not in anchor:
                    continue
                info = anchor[rid]
                seq = fetch(fh, idx[rid])
                if info['wanted_section_of_read'] == 'before_match_start_on_read':
                    ext = seq[:int(info['match_start_on_read'])]
                else:
                    ext = seq[int(info['match_end_on_read']):]
                if len(ext) >= a.min_ext:
                    recs.append((rid, ext, info['pident']))
            if len(recs) < 2:
                summary.append((chr_end, n_tot, len(recs), 0, 0, None, 'too few candidates'))
                continue

            qf = os.path.join(tmp, 'q.fa')
            with open(qf, 'w') as out:
                for rid, ext, _ in recs:
                    out.write(f'>{rid}\n{ext}\n')
            db = os.path.join(tmp, 'db')
            subprocess.run(['makeblastdb', '-in', qf, '-dbtype', 'nucl', '-out', db],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            res = subprocess.run(
                ['blastn', '-query', qf, '-db', db, '-task', 'dc-megablast',
                 '-outfmt', '6 qseqid sseqid pident length', '-num_threads', str(a.threads),
                 '-max_target_seqs', '50'],
                capture_output=True, text=True)

            lens = {rid: len(ext) for rid, ext, _ in recs}
            best = {}
            for line in res.stdout.splitlines():
                f = line.split('\t')
                if len(f) < 4 or f[0] == f[1]:
                    continue
                key = tuple(sorted((f[0], f[1])))
                pid, ln = float(f[2]), int(f[3])
                cov = ln / max(1, min(lens[f[0]], lens[f[1]]))
                score = pid * min(1.0, cov)
                if key not in best or score > best[key][0]:
                    best[key] = (score, pid, ln, cov)

            ids = [r[0] for r in recs]
            for i in range(len(ids)):
                for j in range(i + 1, len(ids)):
                    k = tuple(sorted((ids[i], ids[j])))
                    s, pid, ln, cov = best.get(k, (0.0, 0.0, 0, 0.0))
                    rows.append({'chr_end': chr_end, 'read_a': k[0], 'read_b': k[1],
                                 'similarity': round(s, 2), 'pident': round(pid, 2),
                                 'align_len': ln, 'coverage': round(cov, 3),
                                 'len_a': lens[k[0]], 'len_b': lens[k[1]]})

            def biggest_group(thresh):
                bestn = 1
                for seed in ids:
                    grp = [seed] + [o for o in ids if o != seed and
                                    best.get(tuple(sorted((seed, o))), (0,))[0] >= thresh]
                    ok = [m for m in grp if all(
                        m == o or best.get(tuple(sorted((m, o))), (0,))[0] >= thresh for o in grp)]
                    bestn = max(bestn, len(ok))
                return bestn

            g90 = biggest_group(90)
            g95 = biggest_group(95)
            med = pd.Series([r['similarity'] for r in rows if r['chr_end'] == chr_end]).median()
            summary.append((chr_end, n_tot, len(recs), g90, g95, med, ''))

    df = pd.DataFrame(rows)
    df.to_csv(a.output, sep='\t', index=False)

    print(f"{'chr_end':<8} {'pool':>6} {'cand':>5} {'grp>=90':>8} {'grp>=95':>8} {'med sim':>8}  note")
    print('-' * 66)
    for ce, pool, nc, g90, g95, med, note in summary:
        m = f'{med:.1f}' if med == med and med is not None else '-'
        flag = note if note else ('<== no 3-of-n consensus at 95' if g95 < 3 else '')
        print(f'{ce:<8} {pool:>6} {nc:>5} {g90:>8} {g95:>8} {m:>8}  {flag}')
    if not df.empty:
        print(f'\npairs compared: {len(df)}')
        print(f'similarity  median={df.similarity.median():.1f}  '
              f'p10={df.similarity.quantile(.1):.1f}  min={df.similarity.min():.1f}')
        for t in (99, 97, 95, 90, 80):
            print(f'   pairs >= {t}: {(df.similarity >= t).sum():>4} / {len(df)} '
                  f'({100 * (df.similarity >= t).mean():.0f}%)')
    print(f'\nwrote {a.output}')


if __name__ == '__main__':
    main()
