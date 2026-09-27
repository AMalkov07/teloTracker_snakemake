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

This script makes no new calls: it classifies those paths and aggregates them per end.

Per read (gain-like statuses only):
  NAMED REPEAT        one named donor end (the read's own end counts, as 'self') either
                      gives a strong / moderate circle, or gives two or more pieces.
  UNCONFIRMED REPEAT  not a named repeat, but the path holds a circle that cannot be pinned
                      to one donor: the donor is ambiguous ('chr14L|chr4R', several ends explain
                      it equally) or tentative ('chr14L?[3]', one Y' ID found at several ends),
                      or the circle is weak (one verbatim copy of the donor explains it as well).
                      These are runs of the same Y' ID, kept visible but out of the repeat count.

Per event: reads are collapsed into distinct gained arrays. Two reads are the same array when
they come from the same end, carry the same Y' IDs in the same order, and every ITS agrees
within yprime_path.ITS_TOL bp. One event carried by a clone is then counted once, with the
number of reads that support it, instead of once per read.

Usage:
    python onion_skin_summary.py <base_name> <recombination_dir> <out_tsv> <read_ids_out> \
        [--events-out <events_tsv>]

Writes <out_tsv> (one row per end plus ALL), <read_ids_out> (the gain-like read IDs, for
plot_yprime_copies.py --reads) and, with --events-out, one row per distinct gained array.
"""

import argparse
import csv
import glob
import os
import re
from collections import Counter

import yprime_path as ypath

GAIN_LIKE = ("Y' Gain", "1st Y' Change", "Y' Recombination")
SEG_RE = re.compile(r'^(?P<donor>[^\[:(]+?)(?P<tent>\?)?(?:\[(?P<pos>[^\]]*)\])?:(?P<ids>[^(]*)(?:\((?P<note>.*)\))?$')
CIRC_RE = re.compile(r'circ x[\d.]+(?: (?P<support>strong|moderate|weak))?')

NAMED, UNCONFIRMED, NONE = 'named_repeat', 'unconfirmed_repeat', 'none'

COLS = ['chr_end', 'n_reads', 'n_gain_like', 'n_gain_events', 'n_gain_2plus', 'n_with_path',
        'n_single_donor', 'n_multi_donor',
        'n_named_repeat', 'n_named_repeat_events', 'n_unconfirmed_repeat', 'n_unconfirmed_repeat_events',
        'n_circle', 'n_circle_strong', 'n_circle_moderate', 'n_circle_weak',
        'n_circle_unassigned', 'top_donors', 'top_circles']

EVENT_COLS = ['event_id', 'chr_end', 'n_reads', 'status', 'n_copies', 'divergence_idx',
              'array', 'gained', 'y_prime_path', 'tier', 'repeat_donor', 'read_ids']


# ---------------------------------------------------------------------------
# path parsing
# ---------------------------------------------------------------------------

def parse_segments(path):
    """Pieces of a y_prime_path string, anchor-to-telomere order.

    Each piece: {'donor': name or None, 'n_ids': int, 'circle': bool, 'support': str}.
    donor is None when the piece is tentative ('chr14L?[3]') or ambiguous ('chr12R|chr4R');
    'self' is a named donor (the read's own end)."""
    pieces = []
    for seg in (path or '').split(' > '):
        seg = seg.strip()
        if not seg:
            continue
        m = SEG_RE.match(seg)
        if not m:
            pieces.append({'donor': None, 'n_ids': 1, 'circle': 'circ' in seg, 'support': ''})
            continue
        donor = m.group('donor').strip()
        definite = not m.group('tent') and '|' not in donor and donor != '?'
        c = CIRC_RE.search(m.group('note') or '')
        pieces.append({'donor': donor if definite else None,
                       'n_ids': len([x for x in m.group('ids').split(',') if x]),
                       'circle': bool(c), 'support': (c.group('support') or '') if c else ''})
    return pieces


def parse_path(path):
    """[(donor or None, is_circle)] per piece; donor is None when tentative or ambiguous."""
    return [(p['donor'], p['circle']) for p in parse_segments(path)]


def repeat_tier(pieces):
    """(tier, repeat_donor) for one read's pieces; see the module docstring."""
    per_donor = Counter(p['donor'] for p in pieces if p['donor'])
    for p in pieces:
        if p['donor'] and p['circle'] and p['support'] in ('strong', 'moderate'):
            return NAMED, p['donor']
    for donor, n in per_donor.most_common():
        if n >= 2:
            return NAMED, donor
    if any(p['circle'] for p in pieces):
        return UNCONFIRMED, ''
    return NONE, ''


def parse_circles(field):
    """[(unit_string, support)] from y_prime_path_circles."""
    out = []
    for c in (field or '').split(';'):
        c = c.strip()
        if c and ':' in c:
            unit, support = c.rsplit(':', 1)
            out.append((unit, support))
    return out


# ---------------------------------------------------------------------------
# events: distinct gained arrays
# ---------------------------------------------------------------------------

def read_tokens(r):
    """(ID, ITS) tokens of a read, anchor-to-telomere. Falls back to the ID array alone
    (no ITS) when the read carries no hit coordinates."""
    toks = ypath.read_tokens_from_positions(r.get('y_prime_positions', ''), r.get('telo_side', 'end'))
    if toks:
        return toks
    ids = [x for x in (r.get('y_prime_observed_array') or '').split(',') if x]
    return [(i, None) for i in ids]


def same_array(a, b, tol=ypath.ITS_TOL):
    """Same Y' IDs in the same order and every ITS within tol bp (None matches anything)."""
    if len(a) != len(b) or any(x[0] != y[0] for x, y in zip(a, b)):
        return False
    return all(x[1] is None or y[1] is None or abs(x[1] - y[1]) <= tol for x, y in zip(a, b))


def format_array(tokens, div=-1):
    """'ID1(166) ID2(166) | ID2(167) ID2' -- ITS after each copy; '|' marks the divergence."""
    out = []
    for i, (yid, its) in enumerate(tokens):
        out.append(('| ' if i == div else '') + (f'{yid}({its})' if its is not None else yid))
    return ' '.join(out)


def collapse_events(rows):
    """Group gain-like reads of ONE end into distinct gained arrays.

    Reads are taken in read_id order so the grouping (and the representative whose ITS
    values the others are compared to) does not depend on file order."""
    events = []
    for r in sorted(rows, key=lambda x: x.get('read_id', '')):
        toks = read_tokens(r)
        for ev in events:
            if same_array(ev['tokens'], toks):
                ev['members'].append(r)
                break
        else:
            events.append({'tokens': toks, 'members': [r]})
    return events


def describe_event(ev, end, k):
    members = ev['members']
    rep = members[0]
    path = Counter(m.get('y_prime_path', '') for m in members).most_common(1)[0][0]
    tier, donor = repeat_tier(parse_segments(path))
    try:
        div = int(rep.get('y_prime_divergence_idx', -1))
    except (TypeError, ValueError):
        div = -1
    return {'event_id': f'{end}_E{k}', 'chr_end': end, 'n_reads': len(members),
            'status': Counter(m.get('y_prime_recombination_status', '') for m in members).most_common(1)[0][0],
            'n_copies': len(ev['tokens']), 'divergence_idx': div,
            'array': format_array(ev['tokens'], div), 'gained': rep.get('y_prime_gained_segment', ''),
            'y_prime_path': path, 'tier': tier, 'repeat_donor': donor,
            'read_ids': ','.join(m.get('read_id', '') for m in members)}


# ---------------------------------------------------------------------------
# per-end summary
# ---------------------------------------------------------------------------

def summarise(rows, events=None):
    """Counters for one end (or ALL). `events` are the described events of those rows."""
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
        pieces = parse_segments(r.get('y_prime_path'))
        if not pieces:
            continue
        s['n_with_path'] += 1
        s['n_single_donor' if len(pieces) == 1 else 'n_multi_donor'] += 1
        tier, _ = repeat_tier(pieces)
        if tier == NAMED:
            s['n_named_repeat'] += 1
        elif tier == UNCONFIRMED:
            s['n_unconfirmed_repeat'] += 1
        named_circ = [p for p in pieces if p['circle'] and p['donor']]
        if named_circ:
            s['n_circle'] += 1
            for p in named_circ:
                if p['support']:
                    s[f'n_circle_{p["support"]}'] += 1
        elif any(p['circle'] for p in pieces):
            s['n_circle_unassigned'] += 1
        for unit, _support in parse_circles(r.get('y_prime_path_circles')):
            if '|' not in unit.split(':', 1)[0]:
                circles[unit] += 1
        primary = r.get('y_prime_path_primary_donor') or ''
        if primary and '|' not in primary and not primary.endswith('?'):
            donors[primary] += 1
    for ev in events or []:
        s['n_gain_events'] += 1
        if ev['tier'] == NAMED:
            s['n_named_repeat_events'] += 1
        elif ev['tier'] == UNCONFIRMED:
            s['n_unconfirmed_repeat_events'] += 1
    s['top_donors'] = ', '.join(f'{d}({n})' for d, n in donors.most_common(3))
    s['top_circles'] = ', '.join(f'{u}({n})' for u, n in circles.most_common(3))
    return s


def end_key(name):
    m = re.match(r'chr(\d+)([LR])', name)
    return (int(m.group(1)), m.group(2)) if m else (999, name)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('base_name')
    ap.add_argument('recombination_dir')
    ap.add_argument('out_tsv')
    ap.add_argument('read_ids_out')
    ap.add_argument('--events-out', default='', help='one row per distinct gained array')
    args = ap.parse_args()
    base, recomb_dir = args.base_name, args.recombination_dir

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

    events = {}
    for end in sorted(per_end, key=end_key):
        gains = [r for r in per_end[end] if r.get('y_prime_recombination_status') in GAIN_LIKE]
        events[end] = [describe_event(ev, end, k + 1) for k, ev in enumerate(
            sorted(collapse_events(gains), key=lambda e: (-len(e['members']), e['members'][0].get('read_id', ''))))]
    all_events = [e for end in events for e in events[end]]

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, 'w') as fh:
        fh.write(f'# sample: {base}\n')
        fh.write('# gain-like = y_prime_recombination_status in ' + ', '.join(GAIN_LIKE) + '\n')
        fh.write('# named repeat = one named donor gives a strong/moderate circle or >= 2 pieces\n')
        fh.write('# unconfirmed repeat = a circle whose donor is ambiguous/tentative, or a weak circle\n')
        fh.write(f'# events = distinct gained arrays (same end, same Y\' IDs, every ITS within {ypath.ITS_TOL} bp)\n')
        fh.write('\t'.join(COLS) + '\n')
        for end in sorted(per_end, key=end_key) + ['ALL']:
            s = (summarise(all_rows, all_events) if end == 'ALL'
                 else summarise(per_end[end], events.get(end)))
            fh.write('\t'.join([end] + [str(s.get(c, 0 if c.startswith('n_') else ''))
                                        for c in COLS[1:]]) + '\n')

    if args.events_out:
        os.makedirs(os.path.dirname(os.path.abspath(args.events_out)), exist_ok=True)
        with open(args.events_out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=EVENT_COLS, delimiter='\t', lineterminator='\n')
            w.writeheader()
            for e in all_events:
                w.writerow(e)

    os.makedirs(os.path.dirname(os.path.abspath(args.read_ids_out)), exist_ok=True)
    with open(args.read_ids_out, 'w') as fh:
        fh.write(''.join(f'{i}\n' for i in gain_ids))

    tot = summarise(all_rows, all_events)
    print(f'{base}: {tot["n_gain_like"]} gain-like reads = {tot["n_gain_events"]} distinct arrays over '
          f'{len(per_end)} ends; {tot["n_named_repeat"]} named same-donor repeats '
          f'({tot["n_named_repeat_events"]} arrays), {tot["n_unconfirmed_repeat"]} unconfirmed')
    print(f'Written: {args.out_tsv}')
    if args.events_out:
        print(f'Written: {args.events_out}')


if __name__ == '__main__':
    main()
