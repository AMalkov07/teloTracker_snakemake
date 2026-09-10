#!/usr/bin/env python3
"""
verify_recombination.py -- Part B checks on recombination outputs.

Works on a snapshot directory laid out as
  <snapshot>/<sample>/recombination/<sample>_<chr_end>_features.tsv
  <snapshot>/<sample>/recombination/<sample>_recombination_summary.tsv
  <snapshot>/<sample>/<sample>_post_telo_trimming.tsv
  <snapshot>/<sample>/<sample>_post_y_prime_probe.tsv
and the day-0 reference's simp.bed + extracted Y' library.

Subcommands
  null              day-0 self-run: non-Loss recombination should be ~0 per end;
                    Y' Loss rate is recorded as the per-end baseline
  replicates        two samples of the same library/strain: per-end deltas
  positive-control  a known event (e.g. 7172 chr11L -> chr11R): fraction of reads
                    attributed to the expected source
  loss-vs-truncation  classify every "Y' Loss" read as real_contraction /
                    rm_miss / unconfirmed_end / other using telomere + probe TSVs
  truth-set         build & score the alternating-array truth set (reads at other
                    ends that carry the unique Y' fingerprint of --fingerprint-end)

Every subcommand writes <out>/<name>.tsv and appends to <out>/partB_report.md.
"""

import argparse
import glob
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from verify_day0_reference import parse_lib, end_sort_key, load_bed

CONFIRMED_MIN_REPEAT = 30          # read_summary.tsv "qualifying" definition
LOSS = "Y' Loss"
GAINLIKE = ("Y' Gain", "Y' Recombination", "1st Y' Change")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_features(snapshot, sample):
    frames = []
    for f in sorted(glob.glob(os.path.join(snapshot, sample, 'recombination', f'{sample}_chr*_features.tsv'))):
        if os.path.getsize(f) == 0:
            continue
        try:
            df = pd.read_csv(f, sep='\t')
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df['recombination_source'] = df['recombination_source'].fillna('').astype(str)
    df['y_prime_observed_array'] = df['y_prime_observed_array'].fillna('').astype(str)
    return df


def load_telo(snapshot, sample):
    p = os.path.join(snapshot, sample, f'{sample}_post_telo_trimming.tsv')
    if not os.path.exists(p):
        return pd.DataFrame()
    t = pd.read_csv(p, sep='\t')
    if t.columns[0].startswith('Unnamed'):
        t = t.drop(columns=t.columns[0])
    t['repeat_length'] = pd.to_numeric(t['repeat_length'], errors='coerce')
    t['end_confirmed'] = (t['Adapter_After_Telomere'].astype(str) == 'True') & (t['repeat_length'] >= CONFIRMED_MIN_REPEAT)
    return t[['read_id', 'repeat_length', 'Adapter_After_Telomere', 'end_confirmed']].drop_duplicates('read_id')


def load_probe(snapshot, sample):
    p = os.path.join(snapshot, sample, f'{sample}_post_y_prime_probe.tsv')
    if not os.path.exists(p):
        return pd.DataFrame()
    t = pd.read_csv(p, sep='\t')
    return t[['read_id', 'y_prime_probe_count', 'reference_y_primes']].drop_duplicates('read_id')


def load_summary(snapshot, sample):
    p = os.path.join(snapshot, sample, 'recombination', f'{sample}_recombination_summary.tsv')
    return pd.read_csv(p, sep='\t') if os.path.exists(p) else pd.DataFrame()


def reference_arrays(day0_lib):
    """{chr_end: [ID at pos 1, ID at pos 2, ...]} from the Y' library headers."""
    elem, _, _ = parse_lib(day0_lib)
    arrays = {}
    for (ce, pos), rec in elem.items():
        arrays.setdefault(ce, {})[pos] = rec['id']
    return {ce: [d[k] for k in sorted(d)] for ce, d in arrays.items()}


def ref_counts_from_bed(day0_bed):
    c = Counter()
    for f in load_bed(day0_bed):
        if f['ftype'] == 'y_prime':
            c[f['chr_end']] += 1
    return c


def armless(src):
    return bool(src) and src not in ('ambiguous',) and not re.match(r'^chr\d+[LR]$', src)


def md_table(df, max_rows=80):
    if df is None or df.empty:
        return '_(none)_\n'
    d = df.head(max_rows)
    cols = list(d.columns)
    out = ['| ' + ' | '.join(map(str, cols)) + ' |', '|' + '---|' * len(cols)]
    for _, r in d.iterrows():
        out.append('| ' + ' | '.join('' if pd.isna(v) else str(v) for v in r.values) + ' |')
    if len(df) > max_rows:
        out.append(f'\n_({len(df) - max_rows} more rows in the TSV)_')
    return '\n'.join(out) + '\n'


def append_report(out, title, lines):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'partB_report.md'), 'a') as fh:
        fh.write(f'\n## {title}\n\n' + '\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# B1 null
# ---------------------------------------------------------------------------

def cmd_null(a):
    df = load_features(a.snapshot, a.sample)
    if df.empty or 'y_prime_recombination_status' not in df.columns:
        print(f'cmd_null: no features for {a.sample} in {a.snapshot} -- skipped'); return
    telo = load_telo(a.snapshot, a.sample)
    df = df.merge(telo, on='read_id', how='left')
    rows = []
    for ce, g in df.groupby('chr_end'):
        n = len(g)
        det = g[g['recombination_detected'] == True]
        loss = g[g['y_prime_recombination_status'] == LOSS]
        nonloss = det[det['y_prime_recombination_status'] != LOSS]
        rows.append({
            'chr_end': ce, 'n_reads': n,
            'n_recomb': len(det), 'pct_recomb': round(100 * len(det) / n, 1),
            'n_recomb_nonloss': len(nonloss), 'pct_recomb_nonloss': round(100 * len(nonloss) / n, 1),
            'n_loss': len(loss), 'pct_loss': round(100 * len(loss) / n, 1),
            'n_loss_confirmed_end': int(loss['end_confirmed'].fillna(False).sum()),
            'n_gain': int((g['y_prime_recombination_status'] == "Y' Gain").sum()),
            'n_first_change': int((g['y_prime_recombination_status'] == "1st Y' Change").sum()),
            'n_yp_recomb': int((g['y_prime_recombination_status'] == "Y' Recombination").sum()),
            'n_spacer_switch': int(g['spacer_recombination'].isin(['switch_detected', 'full_switch']).sum()),
            'n_x_switch': int(g['x_element_recombination'].isin(['switch_detected', 'full_switch']).sum()),
            'n_ambiguous': int((det['recombination_source'] == 'ambiguous').sum()),
            'n_armless_source': int(det['recombination_source'].map(armless).sum()),
            'nonloss_sources': ','.join(f'{k}:{v}' for k, v in Counter(nonloss['recombination_source']).most_common(4)),
            'verdict': 'PASS' if (100 * len(nonloss) / n) <= a.max_nonloss_pct else 'FAIL',
        })
    res = pd.DataFrame(rows).sort_values('chr_end', key=lambda s: s.map(end_sort_key))
    res.to_csv(os.path.join(a.out, f'B1_null_{a.sample}.tsv'), sep='\t', index=False)
    fails = res[res['verdict'] == 'FAIL']
    append_report(a.out, f'B1 null (day-0 self-run): {a.sample}', [
        f'Criterion: non-Loss recombination <= {a.max_nonloss_pct}% per end. '
        f"Y' Loss is reported as the per-end baseline (subclonal copy-number heterogeneity), not as a failure.",
        f'**{len(fails)} of {len(res)} ends FAIL**; total reads {res["n_reads"].sum()}; '
        f'arm-less sources: {res["n_armless_source"].sum()}; ambiguous: {res["n_ambiguous"].sum()}.\n',
        '### Failing ends', md_table(fails),
        '### Ends with Y\' Loss baseline > 0', md_table(res[res['n_loss'] > 0][['chr_end', 'n_reads', 'n_loss', 'pct_loss', 'n_loss_confirmed_end']]),
    ])
    print(f'B1 {a.sample}: {len(fails)}/{len(res)} ends fail non-Loss<={a.max_nonloss_pct}%')


# ---------------------------------------------------------------------------
# B2 replicates
# ---------------------------------------------------------------------------

def per_end_rates(df):
    rows = {}
    for ce, g in df.groupby('chr_end'):
        n = len(g)
        det = g[g['recombination_detected'] == True]
        src = det['recombination_source']
        src = src[(src != '') & (src != 'ambiguous')]
        rows[ce] = {'n': n, 'pct_recomb': 100 * len(det) / n,
                    'pct_gain': 100 * (g['y_prime_recombination_status'] == "Y' Gain").sum() / n,
                    'pct_loss': 100 * (g['y_prime_recombination_status'] == LOSS).sum() / n,
                    'top_source': src.mode().iloc[0] if not src.empty else ''}
    return rows


def cmd_replicates(a):
    A, B = load_features(a.snapshot, a.sample_a), load_features(a.snapshot, a.sample_b)
    if A.empty or B.empty:
        print(f'cmd_replicates: missing features for {a.sample_a} or {a.sample_b} -- skipped'); return
    ra, rb = per_end_rates(A), per_end_rates(B)
    rows = []
    for ce in sorted(set(ra) | set(rb), key=end_sort_key):
        x, y = ra.get(ce), rb.get(ce)
        if not x or not y:
            rows.append({'chr_end': ce, 'verdict': 'missing_in_one'})
            continue
        d = y['pct_recomb'] - x['pct_recomb']
        rows.append({'chr_end': ce, 'n_a': x['n'], 'n_b': y['n'],
                     'pct_recomb_a': round(x['pct_recomb'], 1), 'pct_recomb_b': round(y['pct_recomb'], 1), 'delta': round(d, 1),
                     'pct_gain_a': round(x['pct_gain'], 1), 'pct_gain_b': round(y['pct_gain'], 1),
                     'pct_loss_a': round(x['pct_loss'], 1), 'pct_loss_b': round(y['pct_loss'], 1),
                     'top_source_a': x['top_source'], 'top_source_b': y['top_source'],
                     'source_agrees': (x['top_source'] == y['top_source']) or not (x['top_source'] and y['top_source']),
                     'verdict': 'FLAG' if abs(d) > a.max_delta else 'ok'})
    res = pd.DataFrame(rows)
    name = f'B2_replicates_{a.sample_a}__vs__{a.sample_b}'
    res.to_csv(os.path.join(a.out, name + '.tsv'), sep='\t', index=False)
    flagged = res[res['verdict'] == 'FLAG']
    append_report(a.out, f'B2 replicate concordance: {a.sample_a} vs {a.sample_b}', [
        f'Criterion: |delta pct_recombination| <= {a.max_delta} points per end.',
        f'**{len(flagged)} of {len(res)} ends flagged**.\n', md_table(flagged)])
    print(f'B2 {a.sample_a} vs {a.sample_b}: {len(flagged)}/{len(res)} ends flagged')


# ---------------------------------------------------------------------------
# B3 positive control
# ---------------------------------------------------------------------------

def cmd_positive(a):
    rows = []
    for s in a.samples.split(','):
        df = load_features(a.snapshot, s)
        if df.empty:
            continue
        g = df[df['chr_end'] == a.chr_end]
        n = len(g)
        hit = g[g['recombination_source'] == a.expected_source]
        rows.append({'sample': s, 'chr_end': a.chr_end, 'n_reads': n,
                     'n_expected_source': len(hit), 'pct_expected_source': round(100 * len(hit) / max(n, 1), 1),
                     'mean_confidence_expected': round(hit['overall_confidence'].mean(), 3) if len(hit) else float('nan'),
                     'pct_complex': round(100 * (g['is_complex_event'] == True).sum() / max(n, 1), 1),
                     'sources': ','.join(f'{k}:{v}' for k, v in Counter(g['recombination_source']).most_common(5)),
                     'verdict': 'PASS' if (100 * len(hit) / max(n, 1)) >= a.min_pct else 'FAIL'})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(a.out, f'B3_positive_control_{a.chr_end}_{a.expected_source}.tsv'), sep='\t', index=False)
    append_report(a.out, f'B3 positive control: {a.chr_end} -> {a.expected_source}', [
        f'Criterion: >= {a.min_pct}% of {a.chr_end} reads attributed to {a.expected_source}.\n', md_table(res)])
    print('B3:', ' '.join(f"{r['sample']}={r['pct_expected_source']}%({r['verdict']})" for r in rows))


# ---------------------------------------------------------------------------
# B4 loss vs truncation
# ---------------------------------------------------------------------------

def classify_loss(row):
    if not bool(row['end_confirmed']):
        return 'unconfirmed_end'
    pc, rm, ref = row['y_prime_probe_count'], row['y_prime_count_on_read'], row['ref_count']
    if pd.notna(pc):
        if pc >= ref:
            return 'rm_miss'            # probe sees the full array -> RepeatMasker missed copies
        if pc > rm:
            return 'rm_miss'
    if row['tail_bp'] > 500:
        return 'other'                  # non-Y' sequence between the Y' region and the telomere
    return 'real_contraction'


def cmd_loss(a):
    df = load_features(a.snapshot, a.sample)
    if df.empty or 'y_prime_recombination_status' not in df.columns:
        print(f'cmd_loss: no features for {a.sample} in {a.snapshot} -- skipped'); return
    telo, probe = load_telo(a.snapshot, a.sample), load_probe(a.snapshot, a.sample)
    refc = ref_counts_from_bed(a.day0_bed)
    loss = df[df['y_prime_recombination_status'] == LOSS].copy()
    # v2 features carry their own copies of these (C8); use the raw TSVs here
    loss = loss.drop(columns=[c for c in ('y_prime_probe_count', 'telomere_end_confirmed',
                                          'telomere_repeat_length', 'y_prime_tail_bp') if c in loss.columns])
    loss = loss.merge(telo, on='read_id', how='left').merge(probe, on='read_id', how='left')
    loss['end_confirmed'] = loss['end_confirmed'].fillna(False)
    loss['ref_count'] = loss['chr_end'].map(refc).fillna(0).astype(int)
    loss['tail_bp'] = loss.apply(lambda r: (r['y_prime_start'] if r['telo_side'] == 'beginning'
                                            else r['read_length'] - r['y_prime_end']) if r['y_prime_count_on_read'] > 0 else -1, axis=1)
    loss['loss_class'] = loss.apply(classify_loss, axis=1)
    keep = ['read_id', 'chr_end', 'read_length', 'telo_side', 'y_prime_count_on_read', 'ref_count', 'y_prime_probe_count',
            'y_prime_observed_array', 'end_confirmed', 'repeat_length', 'tail_bp', 'recombination_source', 'overall_confidence', 'loss_class']
    loss[keep].to_csv(os.path.join(a.out, f'B4_loss_reads_{a.sample}.tsv'), sep='\t', index=False)
    rows = []
    for ce, g in df.groupby('chr_end'):
        L = loss[loss['chr_end'] == ce]
        c = Counter(L['loss_class'])
        nl = len(L)
        rows.append({'chr_end': ce, 'n_reads': len(g), 'ref_yprimes': refc.get(ce, 0), 'n_loss': nl,
                     'pct_loss': round(100 * nl / len(g), 1),
                     'real_contraction': c.get('real_contraction', 0), 'rm_miss': c.get('rm_miss', 0),
                     'unconfirmed_end': c.get('unconfirmed_end', 0), 'other': c.get('other', 0),
                     'pct_real': round(100 * c.get('real_contraction', 0) / max(nl, 1), 1),
                     'pct_rm_miss': round(100 * c.get('rm_miss', 0) / max(nl, 1), 1),
                     'rm_miss_flag': 'FLAG' if nl >= 5 and 100 * c.get('rm_miss', 0) / nl > a.max_rm_miss_pct else ''})
    res = pd.DataFrame(rows).sort_values('chr_end', key=lambda s: s.map(end_sort_key))
    res.to_csv(os.path.join(a.out, f'B4_loss_summary_{a.sample}.tsv'), sep='\t', index=False)
    append_report(a.out, f"B4 Y' Loss classification: {a.sample}", [
        'real_contraction = telomere end confirmed (adapter after telomere, repeat >= 30 bp), probe count agrees with RepeatMasker, '
        'no long non-Y\' tail; rm_miss = probe count > RepeatMasker count (library/RepeatMasker miss); '
        'unconfirmed_end = read may be truncated; other = >500 bp non-Y\' sequence between Y\' and telomere.\n',
        md_table(res[res['n_loss'] > 0])])
    print(f"B4 {a.sample}: {len(loss)} Loss reads; classes {dict(Counter(loss['loss_class']))}; "
          f"rm_miss flags: {int((res['rm_miss_flag'] == 'FLAG').sum())}")


# ---------------------------------------------------------------------------
# B5 truth set: the alternating fingerprint
# ---------------------------------------------------------------------------

def longest_alternating_window(ids, pair):
    """Longest run inside `ids` that alternates between the two IDs of `pair`
    (either phase). Returns (length, start)."""
    a, b = pair
    best, best_i = 0, -1
    i = 0
    n = len(ids)
    while i < n:
        if ids[i] not in (a, b):
            i += 1
            continue
        j = i
        while j + 1 < n and ids[j + 1] in (a, b) and ids[j + 1] != ids[j]:
            j += 1
        if j - i + 1 > best:
            best, best_i = j - i + 1, i
        i = j + 1
    return best, best_i


def cmd_truth(a):
    arrays = reference_arrays(a.day0_lib)
    fp = arrays.get(a.fingerprint_end, [])
    uniq = list(dict.fromkeys(fp))
    if len(uniq) != 2 or len(fp) < 3:
        sys.exit(f'{a.fingerprint_end} reference array {fp} is not a 2-ID alternating pattern')
    pair = (uniq[0], uniq[1])
    # which other ends could produce an alternating window of length k?
    others = {ce: longest_alternating_window(arr, pair)[0] for ce, arr in arrays.items() if ce != a.fingerprint_end}
    max_other = max(others.values()) if others else 0
    min_len = max(a.min_window, max_other + 1)
    # ends whose own array holds BOTH fingerprint IDs cannot give an unambiguous
    # truth read (their own array can explain a short alternating window)
    excluded_ends = {ce for ce, arr in arrays.items() if ce != a.fingerprint_end and set(pair) <= set(arr)}
    rows = []
    for s in a.samples.split(','):
        df = load_features(a.snapshot, s)
        if df.empty:
            continue
        g = df[(df['chr_end'] != a.fingerprint_end) & (df['y_prime_recombination_status'].isin(GAINLIKE))
               & (~df['chr_end'].isin(excluded_ends))]
        for _, r in g.iterrows():
            ids = [x for x in r['y_prime_observed_array'].split(',') if x]
            L, i = longest_alternating_window(ids, pair)
            if L >= min_len:
                rows.append({'sample': s, 'read_id': r['read_id'], 'chr_end': r['chr_end'],
                             'y_prime_status': r['y_prime_recombination_status'],
                             'observed_array': r['y_prime_observed_array'], 'window_len': L, 'window_start': i,
                             'expected_source': a.fingerprint_end,
                             'recombination_source': r['recombination_source'],
                             'overall_confidence': r['overall_confidence'],
                             'compatible_ends': r.get('y_prime_compatible_ends', ''),
                             'fingerprint_source': r.get('y_prime_fingerprint_source', ''),
                             'mechanism': r.get('recombination_mechanism', ''),
                             'hit': r['recombination_source'] == a.fingerprint_end})
    ts = pd.DataFrame(rows)
    ts.to_csv(os.path.join(a.out, f'B5_truth_set_{a.fingerprint_end}.tsv'), sep='\t', index=False)
    if a.truth_out:
        os.makedirs(os.path.dirname(a.truth_out) or '.', exist_ok=True)
        ts[['sample', 'read_id', 'chr_end', 'observed_array', 'window_len', 'expected_source']].to_csv(a.truth_out, sep='\t', index=False)
    per_sample = ts.groupby('sample').agg(n=('read_id', 'size'), n_hit=('hit', 'sum')).reset_index() if not ts.empty else pd.DataFrame()
    if not per_sample.empty:
        per_sample['pct_hit'] = (100 * per_sample['n_hit'] / per_sample['n']).round(1)
    # primary criterion: a full period (window >= len(fp)) -- unambiguously the fingerprint
    strict = ts[ts['window_len'] >= len(fp)] if not ts.empty else ts
    frac = 100 * strict['hit'].mean() if not strict.empty else float('nan')
    frac_all = 100 * ts['hit'].mean() if not ts.empty else float('nan')
    conf_hit = strict[strict['hit']]['overall_confidence'].mean() if not strict.empty and strict['hit'].any() else float('nan')
    verdict = 'PASS' if (not strict.empty and frac >= a.min_pct) else 'FAIL'
    by_len = ts.groupby('window_len').agg(n=('hit', 'size'), n_hit=('hit', 'sum')).reset_index() if not ts.empty else pd.DataFrame()
    append_report(a.out, f'B5 truth set: alternating {pair[0]}/{pair[1]} fingerprint of {a.fingerprint_end}', [
        f'Reference array at {a.fingerprint_end}: `{",".join(fp)}`. Longest alternating window at any other end: {max_other}; '
        f'ends whose own array holds both IDs are excluded ({", ".join(sorted(excluded_ends)) or "none"}).',
        f'Truth set: reads at other ends with a gain-like Y\' status carrying an alternating window >= {min_len}. '
        f'**Strict set (window >= {len(fp)}, a full period): {len(strict)} reads, {frac:.1f}% attributed to {a.fingerprint_end} '
        f'(mean confidence of hits {conf_hit:.2f}) -> {verdict}** (criterion >= {a.min_pct}%). '
        f'All windows >= {min_len}: {len(ts)} reads, {frac_all:.1f}%.\n',
        '### By window length', md_table(by_len),
        '### Per sample', md_table(per_sample),
        '### Attributed sources', md_table(pd.DataFrame(Counter(ts['recombination_source']).most_common(), columns=['source', 'n']) if not ts.empty else None),
        '### Reads', md_table(ts[['sample', 'read_id', 'chr_end', 'observed_array', 'window_len', 'recombination_source', 'overall_confidence', 'mechanism']], 80)])
    print(f'B5: strict {len(strict)} reads {frac:.1f}% -> {a.fingerprint_end} ({verdict}); all >= {min_len}: {len(ts)} reads {frac_all:.1f}%')


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--snapshot', required=True)
    ap.add_argument('--out', required=True)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('null'); p.add_argument('--sample', required=True); p.add_argument('--max-nonloss-pct', type=float, default=1.0); p.set_defaults(fn=cmd_null)
    p = sub.add_parser('replicates'); p.add_argument('--sample-a', required=True); p.add_argument('--sample-b', required=True); p.add_argument('--max-delta', type=float, default=5.0); p.set_defaults(fn=cmd_replicates)
    p = sub.add_parser('positive-control'); p.add_argument('--samples', required=True); p.add_argument('--chr-end', required=True); p.add_argument('--expected-source', required=True); p.add_argument('--min-pct', type=float, default=90.0); p.set_defaults(fn=cmd_positive)
    p = sub.add_parser('loss-vs-truncation'); p.add_argument('--sample', required=True); p.add_argument('--day0-bed', required=True); p.add_argument('--max-rm-miss-pct', type=float, default=20.0); p.set_defaults(fn=cmd_loss)
    p = sub.add_parser('truth-set'); p.add_argument('--samples', required=True); p.add_argument('--day0-lib', required=True); p.add_argument('--fingerprint-end', default='chr13L'); p.add_argument('--min-window', type=int, default=3); p.add_argument('--min-pct', type=float, default=80.0); p.add_argument('--truth-out', default=''); p.set_defaults(fn=cmd_truth)

    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    a.fn(a)


if __name__ == '__main__':
    main()
