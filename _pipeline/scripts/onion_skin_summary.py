#!/usr/bin/env python3
"""
onion_skin_summary.py -- per-end summary of how gained Y' arrays were built.

"Onion skin" is the pattern of a read gaining several Y' copies that all come from the SAME
donor end: a donor piece copied in tandem (a circle, rolling-circle style) or the same donor
contributing more than one piece of the gained array. The per-read evidence is already in the
recombination features: analyze_features.py (v2) parses every gained array into donor pieces
with yprime_path.py and writes

    y_prime_path          chr14L[1-3]:ID1,ID1,ID2 > chr14L[3-4]:ID2,ID2
    y_prime_path_circles  chr13L[1-3]:ID2,ID1,ID2x1.7:strong

This script only aggregates those columns per chromosome end; it makes no new calls.

A read counts as a SAME-DONOR REPEAT when its path contains a circle, or when one definite
donor end contributes two or more pieces. Pieces marked tentative ('chr14L?[3]', a single Y'
whose ID exists at several ends) or ambiguous ('chr12R|chr4R') never name a donor.

n_circle* come from y_prime_path_circles, which the path parser fills only for circles whose
donor it can name. A tandem copy on a tentative or ambiguous piece ('chr13L?[1]:ID2,ID2,ID2
(circ x3.0 moderate)') is left out of that column by design; it is counted separately as
n_circle_unassigned so the pattern is not silently lost.

Usage:
    python onion_skin_summary.py <base_name> <recombination_dir> <out_tsv> <read_ids_out>

Writes <out_tsv> (one row per end plus ALL) and <read_ids_out> (the gain-like read IDs, for
plot_yprime_copies.py --reads).
"""

import csv
import glob
import os
import re
import sys
from collections import Counter

GAIN_LIKE = ("Y' Gain", "1st Y' Change", "Y' Recombination")
SEG_RE = re.compile(r'^(?P<donor>[^\[:(]+?)(?P<tent>\?)?(?:\[(?P<pos>[^\]]*)\])?:(?P<ids>[^(]*)(?:\((?P<note>.*)\))?$')

COLS = ['chr_end', 'n_reads', 'n_gain_like', 'n_gain_2plus', 'n_with_path',
        'n_single_donor', 'n_multi_donor', 'n_same_donor_repeat',
        'n_circle', 'n_circle_strong', 'n_circle_moderate', 'n_circle_weak',
        'n_circle_unassigned', 'top_donors', 'top_circles']


def parse_path(path):
    """[(donor or None, is_circle)] per piece; donor is None when tentative or ambiguous."""
    pieces = []
    for seg in (path or '').split(' > '):
        seg = seg.strip()
        if not seg:
            continue
        m = SEG_RE.match(seg)
        if not m:
            pieces.append((None, 'circ' in seg))
            continue
        donor = m.group('donor').strip()
        definite = not m.group('tent') and '|' not in donor
        pieces.append((donor if definite else None, 'circ' in (m.group('note') or '')))
    return pieces


def parse_circles(field):
    """[(unit_string, support)] from y_prime_path_circles."""
    out = []
    for c in (field or '').split(';'):
        c = c.strip()
        if c and ':' in c:
            unit, support = c.rsplit(':', 1)
            out.append((unit, support))
    return out


def summarise(rows):
    s = Counter()
    donors, circles = Counter(), Counter()
    for r in rows:
        s['n_reads'] += 1
        if r.get('y_prime_recombination_status') not in GAIN_LIKE:
            continue
        s['n_gain_like'] += 1
        gained = [x for x in (r.get('y_prime_gained_segment') or '').split(',') if x]
        if len(gained) >= 2:
            s['n_gain_2plus'] += 1
        pieces = parse_path(r.get('y_prime_path'))
        if not pieces:
            continue
        s['n_with_path'] += 1
        s['n_single_donor' if len(pieces) == 1 else 'n_multi_donor'] += 1
        circ = parse_circles(r.get('y_prime_path_circles'))
        definite = Counter(d for d, _ in pieces if d)
        if circ or any(n >= 2 for n in definite.values()):
            s['n_same_donor_repeat'] += 1
        if not circ and any(d is None and is_circ for d, is_circ in pieces):
            s['n_circle_unassigned'] += 1
        if circ:
            s['n_circle'] += 1
            for unit, support in circ:
                if support in ('strong', 'moderate', 'weak'):
                    s[f'n_circle_{support}'] += 1
                circles[unit] += 1
        primary = r.get('y_prime_path_primary_donor') or ''
        if primary and '|' not in primary and not primary.endswith('?'):
            donors[primary] += 1
    s['top_donors'] = ', '.join(f'{d}({n})' for d, n in donors.most_common(3))
    s['top_circles'] = ', '.join(f'{u}({n})' for u, n in circles.most_common(3))
    return s


def end_key(name):
    m = re.match(r'chr(\d+)([LR])', name)
    return (int(m.group(1)), m.group(2)) if m else (999, name)


def main():
    base, recomb_dir, out_tsv, ids_out = sys.argv[1:5]
    per_end, all_rows, gain_ids = {}, [], []
    for f in glob.glob(os.path.join(recomb_dir, f'{base}_chr*_features.tsv')):
        end = os.path.basename(f)[len(base) + 1:-len('_features.tsv')]
        try:
            with open(f) as fh:
                rows = list(csv.DictReader(fh, delimiter='\t'))
        except (OSError, csv.Error):
            rows = []
        if not rows or 'y_prime_recombination_status' not in rows[0]:
            continue                      # skipped end (empty features file)
        per_end[end] = rows
        all_rows.extend(rows)
        gain_ids += [r['read_id'] for r in rows if r.get('y_prime_recombination_status') in GAIN_LIKE]

    if all_rows and 'y_prime_path' not in all_rows[0]:
        print("WARNING: features carry no y_prime_path column (legacy attribution?); "
              "donor and circle columns will be empty")

    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    with open(out_tsv, 'w') as fh:
        fh.write(f'# sample: {base}\n')
        fh.write('# gain-like = y_prime_recombination_status in ' + ', '.join(GAIN_LIKE) + '\n')
        fh.write('# same-donor repeat = a circle in the path, or one definite donor giving >= 2 pieces\n')
        fh.write('\t'.join(COLS) + '\n')
        for end in sorted(per_end, key=end_key) + ['ALL']:
            s = summarise(all_rows if end == 'ALL' else per_end[end])
            fh.write('\t'.join([end] + [str(s.get(c, 0 if c.startswith('n_') else ''))
                                        for c in COLS[1:]]) + '\n')

    os.makedirs(os.path.dirname(os.path.abspath(ids_out)), exist_ok=True)
    with open(ids_out, 'w') as fh:
        fh.write(''.join(f'{i}\n' for i in gain_ids))

    tot = summarise(all_rows)
    print(f'{base}: {tot["n_gain_like"]} gain-like reads over {len(per_end)} ends; '
          f'{tot["n_same_donor_repeat"]} same-donor repeats, {tot["n_circle"]} with a circle '
          f'({tot["n_circle_strong"]} strong)')
    print(f'Written: {out_tsv}')


if __name__ == '__main__':
    main()
