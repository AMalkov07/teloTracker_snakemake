#!/usr/bin/env python3
"""
onion_skin_timecourse.py -- follow Y' "onion skin" layers across the timepoints of one run.

The per-sample onion-skin step (onion_skin_summary.py) reduces every sample to its distinct
gained Y' arrays. A single read shows one snapshot of an end; this step asks whether an end
KEEPS gaining copies from the same donor as time goes on, by comparing arrays between samples.

Array B EXTENDS array A (same end) when B is A plus more Y' copies at the telomere side: A's
IDs are a prefix of B's and every ITS inside A agrees within yprime_path.ITS_TOL bp. The
reference array of the end is the root, so "reference + k copies" is a first layer. Where
B's added copies come from the same named donor as the layer beneath them, the extension is a
SAME-DONOR LAYER (the read's own end counts as a donor, 'self'). A chain of same-donor layers
whose first appearances do not go backwards in time is the onion-skin signature:

    chr13L:  REF ID1(166) ID2  ->  +1 ID2 (day 0)  ->  +2 ID2 (day 3)  ->  +5 ID2 (day 5)

Separate samples cannot prove that one lineage grew; coexisting subclones give the same
picture. The time order of first appearance is reported so the two readings can be weighed.

Samples must be given in time order.

Inputs per sample: <results>/<sample>/_pipeline/recombination_events/<sample>_onion_skin_events.tsv
and ..._onion_skin_summary.tsv (read counts per end), plus the day-0 BED and Y' library that
the recombination step used (to rebuild the reference arrays).

Outputs (in --out-dir):
  <prefix>_onion_arrays.tsv      every distinct gained array, pooled over samples, with reads per sample
  <prefix>_onion_layers.tsv      each direct extension (A -> B), its added copies, donor and junction ITS
  <prefix>_onion_ladders.tsv     per end and sample: reads that are the reference plus k copies
  <prefix>_onion_timecourse.tsv  per end: the longest same-donor chain and whether it is time-ordered
  <prefix>_<end>_onion_ladder.png  one figure per end with a ladder in at least two samples
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

import yprime_path as ypath

TOL = ypath.ITS_TOL
MAX_K_BIN = 4           # ladder figure bins: +1, +2, +3, +4 or more


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------

TOKEN_RE = re.compile(r'^(?P<id>[^()]+?)(?:\((?P<its>-?\d+)\))?$')


def parse_array(text):
    """'ID1(166) | ID2(167) ID2' -> [('ID1', 166), ('ID2', 167), ('ID2', None)]."""
    toks = []
    for t in (text or '').replace('|', ' ').split():
        m = TOKEN_RE.match(t)
        if m:
            toks.append((m.group('id'), int(m.group('its')) if m.group('its') is not None else None))
    return toks


def read_tsv(path):
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith('#')]
    return list(csv.DictReader(lines, delimiter='\t'))


def load_reference_tokens(bed, y_prime_lib, id_level='family'):
    """{chr_end: [(ID, ITS), ...]} exactly as the recombination step builds them."""
    import analyze_features as af
    af.Y_PRIME_ID_LEVEL = id_level
    return af.build_reference_tokens(bed, af.build_y_prime_info(y_prime_lib))


def find_sample_file(results, sample, name):
    for d in (os.path.join(results, sample, '_pipeline', 'recombination_events'),
              os.path.join(results, sample, 'recombination_events')):
        p = os.path.join(d, f'{sample}_{name}')
        if os.path.isfile(p):
            return p
    return ''


# ---------------------------------------------------------------------------
# arrays and extensions
# ---------------------------------------------------------------------------

def its_agree(a, b):
    return a is None or b is None or abs(a - b) <= TOL


def extends(a, b):
    """Index where b's added copies start if b = a + copies at the telomere side, else -1.
    The ITS after a's last copy is b's junction (new sequence), so it is not compared."""
    n = len(a)
    if len(b) <= n or any(x[0] != y[0] for x, y in zip(a, b[:n])):
        return -1
    if not all(its_agree(x[1], y[1]) for x, y in zip(a[:n - 1], b[:n - 1])):
        return -1
    return n


def segment_spans(path, div):
    """[(start, end_exclusive, donor or None)] of each path piece on the array's token index."""
    import onion_skin_summary as oss
    spans, i = [], max(div, 0)
    for p in oss.parse_segments(path):
        spans.append((i, i + p['n_ids'], p['donor']))
        i += p['n_ids']
    return spans


def donors_over(spans, lo, hi):
    """Named donors of the pieces overlapping token range [lo, hi); None marks an unnamed piece."""
    return [d for s, e, d in spans if s < hi and e > lo]


def pool_arrays(samples, events_by_sample):
    """Merge each sample's events into arrays shared across samples (same end, same IDs, ITS
    within TOL). Returns a list of arrays in first-seen order."""
    arrays = []
    for si, sample in enumerate(samples):
        for ev in events_by_sample.get(sample, []):
            toks = parse_array(ev['array'])
            if not toks:
                continue
            for a in arrays:
                if a['end'] == ev['chr_end'] and len(a['tokens']) == len(toks) and all(
                        x[0] == y[0] and its_agree(x[1], y[1]) for x, y in zip(a['tokens'], toks)):
                    break
            else:
                a = {'end': ev['chr_end'], 'tokens': toks, 'reads': Counter(), 'first': si,
                     'path': ev['y_prime_path'], 'div': int(ev.get('divergence_idx') or -1),
                     'tier': ev.get('tier', ''), 'repeat_donor': ev.get('repeat_donor', '')}
                arrays.append(a)
            a['reads'][sample] += int(ev['n_reads'])
    per_end = defaultdict(int)
    for a in arrays:
        per_end[a['end']] += 1
        a['id'] = f"{a['end']}_A{per_end[a['end']]}"
        a['spans'] = segment_spans(a['path'], a['div'])
    return arrays


MAX_ITS = 400           # bp; a longer (or negative) gap is a missed / partial Y', not a junction ITS


def minimal_period(ids):
    for p in range(1, len(ids) + 1):
        if all(ids[i] == ids[i - p] for i in range(p, len(ids))):
            return p
    return len(ids)


def beneath_unit(node, ref_tokens_end):
    """(unit IDs, known ITS values) of the layer directly beneath any copies added to `node`.

    For the reference it is the terminal Y' copy and the ITS in front of it. For a gained
    array it is the last piece of its path: the copies of that piece and the ITS inside and
    in front of it."""
    toks = node['tokens']
    if not toks:
        return [], []
    if node['id'] == 'REF' or not node['spans']:
        s, e = len(toks) - 1, len(toks)
    else:
        s, e = node['spans'][-1][0], node['spans'][-1][1]
        s, e = max(0, min(s, len(toks) - 1)), min(e, len(toks))
    ids = [t[0] for t in toks[s:e]]
    its = [t[1] for t in toks[max(s - 1, 0):e - 1] if t[1] is not None and 0 <= t[1] <= MAX_ITS]
    return ids, its


def continues_unit(unit, added):
    """True when `added` carries on the repeating pattern of `unit` (phase included)."""
    if not unit or not added:
        return False
    p = minimal_period(unit)
    return all(a == unit[len(unit) - p + (j % p)] for j, a in enumerate(added))


def classify_layer(a, b, k):
    """'repeat'  the added copies continue the unit beneath, and every new ITS (the junction
                 onto the array beneath and the ITS between added copies) matches an ITS
                 already in that unit; with no ITS known beneath, the new ITS must agree
                 with each other;
       'repeat_its_mismatch'  same Y' unit, but a new ITS differs or is a gap;
       'different'  the added copies are another Y' unit.
    Returns (type, junction ITS)."""
    unit, its_known = beneath_unit(a, None)
    added = [t[0] for t in b['tokens'][k:]]
    if not continues_unit(unit, added):
        return 'different', None
    junction = b['tokens'][k - 1][1] if k > 0 else None
    new_its = [t[1] for t in b['tokens'][max(k - 1, 0):len(b['tokens']) - 1] if t[1] is not None]
    if any(not 0 <= x <= MAX_ITS for x in new_its):
        return 'repeat_its_mismatch', junction
    ref_its = its_known or new_its[:1]
    if all(any(abs(x - r) <= TOL for r in ref_its) for x in new_its):
        return 'repeat', junction
    return 'repeat_its_mismatch', junction


def build_layers(arrays, ref_tokens):
    """Direct extensions A -> B within each end (A may be the reference, id 'REF').

    Direct = no third array C with A -> C -> B, so a ladder +1 -> +2 -> +3 is reported step by
    step. Each extension is classified by classify_layer (does it repeat the Y' unit beneath?).
    added_donor is the donor the path parser gave B's added copies, for reference only: the
    parser explains each read on its own and can prefer one ambiguous piece over two pieces
    from the same end, so it is not used to decide whether a layer repeats the one beneath."""
    by_end = defaultdict(list)
    for a in arrays:
        by_end[a['end']].append(a)
    layers = []
    for end, arrs in by_end.items():
        ref = {'id': 'REF', 'end': end, 'tokens': ref_tokens.get(end, []), 'first': -1,
               'spans': [], 'reads': Counter()}
        nodes = [ref] + arrs
        ext = {}
        for a in nodes:
            for b in arrs:
                if a is b:
                    continue
                k = extends(a['tokens'], b['tokens']) if a['tokens'] else (0 if a is ref else -1)
                if k >= 0:
                    ext[(a['id'], b['id'])] = k
        ids = {n['id']: n for n in nodes}
        for (aid, bid), k in ext.items():
            if any((aid, c) in ext and (c, bid) in ext for c in ids if c not in (aid, bid)):
                continue                                   # not direct
            a, b = ids[aid], ids[bid]
            added = b['tokens'][k:]
            kind, junction = classify_layer(a, b, k)
            add_donors = donors_over(b['spans'], k, len(b['tokens']))
            layers.append({
                'chr_end': end, 'from_array': aid, 'to_array': bid,
                'from_n_copies': len(a['tokens']), 'n_added': len(added),
                'added': ' '.join(t[0] for t in added),
                'junction_its': junction if junction is not None else '',
                'layer_type': kind,
                'added_donor': '|'.join(sorted({d or '?' for d in add_donors})),
                'from_first_seen': a['first'], 'to_first_seen': b['first'],
            })
    return layers


def chain_time_ordered(chain):
    firsts = [l['to_first_seen'] for l in chain]
    return all(x <= y for x, y in zip(firsts, firsts[1:]))


def repeat_chains(end, layers):
    """(longest time-ordered chain, length of the longest chain in any order) of 'repeat'
    layers starting at the reference."""
    nxt = defaultdict(list)
    for l in layers:
        if l['chr_end'] == end and l['layer_type'] == 'repeat':
            nxt[l['from_array']].append(l)
    best_ordered, longest = [], 0

    def walk(node, chain, seen):
        nonlocal best_ordered, longest
        longest = max(longest, len(chain))
        if chain_time_ordered(chain) and len(chain) > len(best_ordered):
            best_ordered = list(chain)
        for l in nxt.get(node, []):
            if l['to_array'] not in seen:
                walk(l['to_array'], chain + [l], seen | {l['to_array']})
    walk('REF', [], {'REF'})
    return best_ordered, longest


def describe_chain(chain, samples):
    if not chain:
        return ''
    parts = ['REF']
    for l in chain:
        added = l['added'].split()
        what = f'{len(added)}x {added[0]}' if len(set(added)) == 1 else l['added']
        parts.append(f"+{what} ({samples[l['to_first_seen']]}, ITS {l['junction_its'] or '-'})")
    return ' -> '.join(parts)


# ---------------------------------------------------------------------------
# ladders: reads that are the reference plus k copies
# ---------------------------------------------------------------------------

def ladders(arrays, ref_tokens, samples, n_reads):
    """Per (end, sample): reads whose array is the reference plus k extra copies (by k), and
    how many of those are REPEATS of the reference's terminal Y' (classify_layer 'repeat')."""
    rows = []
    for end in sorted({a['end'] for a in arrays}, key=end_key):
        ref = {'id': 'REF', 'tokens': ref_tokens.get(end, []), 'spans': []}
        ext = []
        for a in arrays:
            if a['end'] != end:
                continue
            k0 = extends(ref['tokens'], a['tokens']) if ref['tokens'] else 0
            if k0 < 0:
                continue
            kind = classify_layer(ref, a, k0)[0] if ref['tokens'] else 'different'
            ext.append((len(a['tokens']) - k0, kind == 'repeat', a))
        for s in samples:
            total = n_reads.get((s, end), 0)
            by_k, rep_k = Counter(), Counter()
            for k, rep, a in ext:
                r = a['reads'].get(s, 0)
                by_k[k] += r
                if rep:
                    rep_k[k] += r
            n_ext, n_rep = sum(by_k.values()), sum(rep_k.values())
            rows.append({'chr_end': end, 'sample': s, 'n_reads': total, 'n_extended': n_ext,
                         'pct_extended': round(100 * n_ext / total, 2) if total else '',
                         'n_repeat_extended': n_rep,
                         'pct_repeat_extended': round(100 * n_rep / total, 2) if total else '',
                         'repeat_by_extra_copies': ';'.join(f'+{k}:{v}' for k, v in sorted(rep_k.items()) if v),
                         'all_by_extra_copies': ';'.join(f'+{k}:{v}' for k, v in sorted(by_k.items()) if v)})
    return rows


def end_key(name):
    m = re.match(r'chr(\d+)([LR])', name)
    return (int(m.group(1)), m.group(2)) if m else (999, name)


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------

def plot_ladder(end, rows, samples, out_png):
    """Stacked bars per sample: % of the end's reads that are the reference plus k repeats of
    its terminal Y' (k binned 1, 2, 3, 4+)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels = ['+1', '+2', '+3', f'+{MAX_K_BIN} or more']
    alphas = [0.35, 0.6, 0.8, 1.0]
    fig, ax = plt.subplots(figsize=(7.5, 0.55 * len(samples) + 2.1))
    for i, r in enumerate(rows):
        bins = Counter()
        for part in filter(None, r['repeat_by_extra_copies'].split(';')):
            k, v = part[1:].split(':')
            bins[min(int(k), MAX_K_BIN)] += int(v)
        left = 0.0
        for b in range(1, MAX_K_BIN + 1):
            w = 100 * bins[b] / r['n_reads'] if r['n_reads'] else 0
            ax.barh(i, w, left=left, color='#C23F39', alpha=alphas[b - 1], height=0.62,
                    label=labels[b - 1] if i == 0 else None)
            left += w
        ax.text(left + 0.3, i, f"{left:.1f}%  ({r['n_repeat_extended']} of {r['n_reads']} reads)",
                va='center', fontsize=8.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r['sample'] for r in rows], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel(f"% of {end} reads = reference array + extra copies of its terminal Y'")
    top = max((float(r['pct_repeat_extended'] or 0) for r in rows), default=0)
    ax.set_xlim(0, max(5.0, top * 1.6))
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(title='extra copies', fontsize=8, title_fontsize=8, frameon=False,
              loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=4)
    ax.set_title(f"{end}: repeats of the terminal Y' added over time", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def write_tsv(path, rows, cols):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--results', default='results', help='results directory holding one dir per sample')
    ap.add_argument('--samples', nargs='+', required=True, help='sample names, in time order')
    ap.add_argument('--day0-bed', required=True)
    ap.add_argument('--y-prime-lib', required=True)
    ap.add_argument('--y-prime-id-level', default='family', choices=['family', 'variant'])
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--prefix', default='timecourse')
    ap.add_argument('--no-plots', action='store_true')
    args = ap.parse_args()

    samples, events, n_reads = [], {}, {}
    for s in args.samples:
        ev = find_sample_file(args.results, s, 'onion_skin_events.tsv')
        sm = find_sample_file(args.results, s, 'onion_skin_summary.tsv')
        if not ev or not sm:
            print(f'WARNING: {s}: no onion-skin events/summary (run the recombination step first); skipped')
            continue
        samples.append(s)
        events[s] = read_tsv(ev)
        for r in read_tsv(sm):
            if r['chr_end'] != 'ALL':
                n_reads[(s, r['chr_end'])] = int(r['n_reads'])
    if len(samples) < 2:
        print('Fewer than two samples with onion-skin output: nothing to compare across time.')
        return 0

    ref_tokens = load_reference_tokens(args.day0_bed, args.y_prime_lib, args.y_prime_id_level)
    arrays = pool_arrays(samples, events)
    by_id = {a['id']: a for a in arrays}
    layers = build_layers(arrays, ref_tokens)
    ladder_rows = ladders(arrays, ref_tokens, samples, n_reads)

    os.makedirs(args.out_dir, exist_ok=True)
    out = lambda name: os.path.join(args.out_dir, f'{args.prefix}_{name}')

    arr_rows = []
    for a in arrays:
        row = {'array_id': a['id'], 'chr_end': a['end'], 'n_copies': len(a['tokens']),
               'array': ' '.join(f'{i}({t})' if t is not None else i for i, t in a['tokens']),
               'y_prime_path': a['path'], 'tier': a['tier'], 'repeat_donor': a['repeat_donor'],
               'first_seen': samples[a['first']]}
        row.update({f'reads_{s}': a['reads'].get(s, 0) for s in samples})
        arr_rows.append(row)
    write_tsv(out('onion_arrays.tsv'), arr_rows,
              ['array_id', 'chr_end', 'n_copies', 'array', 'y_prime_path', 'tier', 'repeat_donor',
               'first_seen'] + [f'reads_{s}' for s in samples])

    for l in layers:
        l['from_first_seen'] = samples[l['from_first_seen']] if l['from_first_seen'] >= 0 else 'reference'
        l['to_first_seen_name'] = samples[l['to_first_seen']]
    write_tsv(out('onion_layers.tsv'),
              sorted(layers, key=lambda l: (end_key(l['chr_end']), l['from_n_copies'], l['to_array'])),
              ['chr_end', 'from_array', 'to_array', 'from_n_copies', 'n_added', 'added', 'junction_its',
               'layer_type', 'added_donor', 'from_first_seen', 'to_first_seen_name'])
    write_tsv(out('onion_ladders.tsv'), ladder_rows,
              ['chr_end', 'sample', 'n_reads', 'n_extended', 'pct_extended', 'n_repeat_extended',
               'pct_repeat_extended', 'repeat_by_extra_copies', 'all_by_extra_copies'])

    tc_rows = []
    for end in sorted({a['end'] for a in arrays}, key=end_key):
        chain, longest = repeat_chains(end, layers)
        lad = [r for r in ladder_rows if r['chr_end'] == end]
        n_rep_samples = sum(1 for r in lad if r['n_repeat_extended'])
        max_k = max((int(p.split(':')[0][1:]) for r in lad
                     for p in filter(None, r['repeat_by_extra_copies'].split(';'))), default=0)
        tc_rows.append({
            'chr_end': end,
            'n_arrays': sum(1 for a in arrays if a['end'] == end),
            'n_layers': sum(1 for l in layers if l['chr_end'] == end),
            'n_repeat_layers': sum(1 for l in layers if l['chr_end'] == end and l['layer_type'] == 'repeat'),
            'repeat_chain_len': len(chain),
            'longest_chain_any_order': longest,
            'chain': describe_chain(chain, samples),
            'pct_repeat_extended_by_sample': ' | '.join(str(r['pct_repeat_extended']) for r in lad),
            'max_repeat_copies': max_k,
        })
        if not args.no_plots and n_rep_samples >= 2 and max_k >= 2:
            try:
                plot_ladder(end, lad, samples, out(f'{end}_onion_ladder.png'))
            except Exception as e:                           # a figure must never sink the tables
                print(f'WARNING: {end}: ladder figure failed: {e}')
    tc_rows.sort(key=lambda r: (-r['repeat_chain_len'], -r['max_repeat_copies'], end_key(r['chr_end'])))
    with open(out('onion_timecourse.tsv'), 'w') as fh:
        fh.write('# samples in time order: ' + ', '.join(samples) + '\n')
        fh.write("# repeat layer = copies added at the telomere side that continue the Y' unit beneath, "
                 f"with a junction ITS matching that unit's ITS (+-{TOL} bp)\n")
        fh.write('# chain = the longest run of repeat layers from the reference whose first appearances '
                 'never go back in time\n')
    with open(out('onion_timecourse.tsv'), 'a', newline='') as fh:
        cols = ['chr_end', 'n_arrays', 'n_layers', 'n_repeat_layers', 'repeat_chain_len',
                'longest_chain_any_order', 'chain', 'pct_repeat_extended_by_sample', 'max_repeat_copies']
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in tc_rows:
            w.writerow(r)

    print(f'{len(samples)} samples, {len(arrays)} distinct gained arrays, {len(layers)} extensions '
          f"({sum(1 for l in layers if l['layer_type'] == 'repeat')} repeat layers)")
    for r in tc_rows[:5]:
        if r['repeat_chain_len'] >= 2:
            print(f"  {r['chr_end']}: {r['chain']}")
    print(f'Written: {args.out_dir}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
