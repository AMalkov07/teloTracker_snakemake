#!/usr/bin/env python3
"""
confidence_compare.py -- old overall_confidence vs the v3 split scores, on the same reads.

Both are computed in one pass (the old score is still written, deprecated), so every read has
both and the comparison is like-for-like. Checks:

  0. Regression: the old score and the donor calls are identical to the reference snapshot, so
     the new code changed nothing else.
  1. Shape: how many reads sit on a hardcoded value, how many distinct values.
  2. Per-end summaries: how much of each per-end mean is just the recombination rate (R^2).
  3. Donor calls vs truth (AUC): 7172 chr11L -> chr11R positive control; chr13L strict set.
  4. Real events vs day-0 calls (AUC).
  5. Worked examples, one per kind of read, with before / after and the evidence behind them.

Usage:
  python verification/confidence_compare.py <new_snapshot> <reference_snapshot> <truth_set_tsv> <out_dir>
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confidence_audit import auc, md, summaries, STRICT_WINDOW  # noqa: E402

NO_YP = "no Y' change"


def load(snapshot):
    frames = []
    for d in sorted(glob.glob(os.path.join(snapshot, '*'))):
        s = os.path.basename(d)
        for f in glob.glob(os.path.join(d, 'recombination', f'{s}_chr*_features.tsv')):
            try:
                df = pd.read_csv(f, sep='\t', dtype=str, keep_default_na=False)
            except pd.errors.EmptyDataError:
                continue
            if df.empty:
                continue
            df['sample'] = s
            frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df['recomb'] = df['recombination_detected'].str.lower() == 'true'
    for c in ('overall_confidence', 'recombination_confidence', 'donor_confidence'):
        if c in df:
            df[c + '_n'] = pd.to_numeric(df[c], errors='coerce')
    df['grp'] = df['sample'].map(lambda s: 'day0' if '_day0_' in s else ('survivor' if 'survivor' in s else 'later'))
    return df


def r2_vs_rate(t, col):
    t = t.dropna(subset=[col, 'pct_recombination'])
    p = t['pct_recombination'] / 100.0
    y = t[col]
    X = np.column_stack([1 - p, p])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    ss_res = ((y - X @ coef) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    return 1 - ss_res / ss_tot if ss_tot else float('nan'), len(t)


def pick(df, mask, sort_col=None, ascending=False):
    sub = df[mask]
    if sub.empty:
        return None
    if sort_col:
        sub = sub.sort_values(sort_col, ascending=ascending)
    return sub.iloc[0]


def main():
    new_snap, ref_snap, truth_tsv, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    new, ref = load(new_snap), load(ref_snap)
    rep = ['# Recombination confidence: old vs new', '',
           f'Same reads, one pass: `{new_snap}` ({new["sample"].nunique()} samples, {len(new):,} reads). '
           f'`overall_confidence` is the old score (still written, deprecated); '
           '`recombination_confidence` and `donor_confidence` are the new ones.', '']

    # ---- 0. regression -------------------------------------------------------------------
    key = ['sample', 'chr_end', 'read_id']
    j = new[key + ['overall_confidence', 'recombination_source']].merge(
        ref[key + ['overall_confidence', 'recombination_source']], on=key, suffixes=('_new', '_ref'))
    same_conf = (j['overall_confidence_new'] == j['overall_confidence_ref']).mean()
    same_src = (j['recombination_source_new'] == j['recombination_source_ref']).mean()
    rep += ['## 0. Nothing else changed', '',
            f'{len(j):,} reads matched to `{ref_snap}`: old score identical for **{100*same_conf:.3f}%**, '
            f'donor call identical for **{100*same_src:.3f}%**.', '']

    rec = new[new['recomb']]
    unc = new[~new['recomb']]

    # ---- 1. shape --------------------------------------------------------------------------
    def shape(s, consts):
        s = s.dropna()
        return (len(s), s.nunique(), 100 * s.round(4).isin(consts).mean(), s.mean(), s.median())
    rows = [
        ['old overall_confidence', 'recombinant', *shape(rec['overall_confidence_n'], [0.3, 0.33])],
        ['new recombination_confidence', 'recombinant', *shape(rec['recombination_confidence_n'], [])],
        ['new donor_confidence', 'recombinant', *shape(rec['donor_confidence_n'], [])],
        ['old overall_confidence', 'unchanged', *shape(unc['overall_confidence_n'], [0.95, 0.9])],
        ['new recombination_confidence', 'unchanged', *shape(unc['recombination_confidence_n'], [0.95, 0.9])],
    ]
    t1 = pd.DataFrame(rows, columns=['score', 'reads', 'n', 'distinct values', '% on a hardcoded value',
                                     'mean', 'median'])
    rep += ['## 1. Shape', '',
            'Hardcoded values: 0.30 / 0.33 (old recombinant floor), 0.95 / 0.90 (old unchanged constants). '
            'The new unchanged-read score still takes two values by design: 0.95 when the read reaches '
            'the telomere, 0.70 when it stops early -- which, unlike before, is information.', '',
            md(t1, 3), '']

    # ---- 2. per-end summaries --------------------------------------------------------------
    t = summaries(new_snap)
    for c in ('mean_recombination_confidence', 'mean_donor_confidence'):
        t[c] = pd.to_numeric(t.get(c), errors='coerce')
    r_old, n_old = r2_vs_rate(t, 'mean_confidence')
    r_rc, _ = r2_vs_rate(t, 'mean_recombination_confidence')
    r_dc, _ = r2_vs_rate(t, 'mean_donor_confidence')
    rep += ['## 2. Per-end summary columns', '',
            f'R^2 of each per-end mean against the recombination rate, over {n_old:,} (sample, end) pairs. '
            'High = the column mostly restates the rate printed next to it.', '',
            md(pd.DataFrame([['mean_confidence (old)', r_old],
                             ['mean_recombination_confidence (new)', r_rc],
                             ['mean_donor_confidence (new)', r_dc]], columns=['column', 'R^2']), 3), '']

    # ---- 3. donor vs truth -----------------------------------------------------------------
    pc = rec[rec['sample'].isin(['7172_day4_with_selection', '7172_day6_with_selection',
                                 '7172_day9_with_selection']) & (rec['chr_end'] == 'chr11L')]
    ts = pd.read_csv(truth_tsv, sep='\t', dtype=str, keep_default_na=False)
    ts = ts[pd.to_numeric(ts['window_len'], errors='coerce') >= STRICT_WINDOW]
    tsj = ts[['sample', 'read_id']].merge(rec, on=['sample', 'read_id'])
    tsj['hit'] = tsj['recombination_source'] == 'chr13L'
    rows = []
    for name, d, ok in (('positive control 7172 chr11L -> chr11R', pc, pc['recombination_source'] == 'chr11R'),
                        ('chr13L truth set, strict', tsj, tsj['hit'])):
        rows.append([name, int(ok.sum()), int((~ok).sum()),
                     auc(d[ok]['overall_confidence_n'], d[~ok]['overall_confidence_n']),
                     auc(d[ok]['donor_confidence_n'], d[~ok]['donor_confidence_n']),
                     d[ok]['donor_confidence_n'].mean(), d[~ok]['donor_confidence_n'].mean()])
    t3 = pd.DataFrame(rows, columns=['truth set', 'correct', 'wrong', 'AUC old', 'AUC new (donor)',
                                     'new: mean correct', 'new: mean wrong'])
    rep += ['## 3. Does the donor score separate right from wrong donors?', '', md(t3, 3), '',
            'The chr13L strict set has very few wrong calls, so its AUC is noisy either way.', '']

    # ---- 4. real vs day-0 ------------------------------------------------------------------
    d0, later = rec[rec['grp'] == 'day0'], rec[rec['grp'] == 'later']
    t4 = pd.DataFrame([
        ['old overall_confidence', auc(later['overall_confidence_n'], d0['overall_confidence_n']),
         d0['overall_confidence_n'].median(), later['overall_confidence_n'].median()],
        ['new recombination_confidence', auc(later['recombination_confidence_n'], d0['recombination_confidence_n']),
         d0['recombination_confidence_n'].median(), later['recombination_confidence_n'].median()],
    ], columns=['score', 'AUC later vs day-0', 'median day-0', 'median later'])
    rep += ['## 4. Does the call score separate real events from day-0 calls?', '',
            'Day-0 calls include genuine standing variation, so no score can reach 1.0 here.', '',
            md(t4, 3), '']

    # ---- 5. worked examples ----------------------------------------------------------------
    st = rec['y_prime_recombination_status']
    tc = rec['telomere_end_confirmed'].str.lower() == 'true'
    sw = lambda c: rec[c].isin(['switch_detected', 'full_switch'])
    ex = [
        ("Y'-only donor transfer that sat on the old 0.30 floor",
         (rec['overall_confidence_n'].round(2) == 0.3) & (st == "Y' Gain") & (rec['recombination_mechanism'] == 'donor_transfer'),
         'donor_confidence_n', False),
        ("Spacer switch + Y' gain, where a weak spacer switch hid the Y' evidence",
         sw('spacer_recombination') & (st == "Y' Gain") & (rec['overall_confidence_n'] < 0.1),
         'donor_confidence_n', False),
        ('Positive control: chr11L -> chr11R', rec.index.isin(pc[pc['recombination_source'] == 'chr11R'].index),
         'donor_confidence_n', False),
        ("Y' Loss on a read that never reaches the telomere (day 0)",
         (st == "Y' Loss") & ~tc & (rec['grp'] == 'day0'), 'recombination_confidence_n', True),
        ('Donor left ambiguous', rec['recombination_source'] == 'ambiguous', 'recombination_confidence_n', False),
        ('Tandem amplification of the read\'s own end', rec['recombination_mechanism'] == 'tandem_amplification_same_end',
         'donor_confidence_n', False),
    ]
    rep += ['## 5. Worked examples: before and after', '']
    for title, mask, col, asc in ex:
        r = pick(rec, mask, col, asc)
        if r is None:
            continue
        rep += [f'### {title}', '',
                f'`{r["sample"]}` {r["chr_end"]} read `{r["read_id"]}`', '',
                f'- Y\' array: `{r.get("y_prime_observed_array", "")}` ({r["y_prime_recombination_status"] or NO_YP}; '
                f'telomere reached: {r.get("telomere_end_confirmed", "")})',
                f'- spacer: {r.get("spacer_recombination", "")} (conf {r.get("spacer_confidence", "")}); '
                f'X element: {r.get("x_element_recombination", "")} (conf {r.get("x_element_confidence", "")})',
                f'- Y\' fingerprint: `{r.get("y_prime_fingerprint_source", "") or "-"}` '
                f'(length {r.get("y_prime_fingerprint_len", "")}); path: `{r.get("y_prime_path", "") or "-"}`',
                f'- donor named: **{r["recombination_source"]}** ({r.get("recombination_mechanism", "")})', '',
                '| | old | new |', '|---|---|---|',
                f'| score | `overall_confidence` = **{r["overall_confidence"]}** | '
                f'`recombination_confidence` = **{r["recombination_confidence"]}**, '
                f'`donor_confidence` = **{r["donor_confidence"]}** |',
                f'| why | `{r.get("confidence_factors", "")}` | `{r.get("confidence_basis", "")}` |', '']
    # an unchanged read that stops before the telomere
    r = pick(new, (~new['recomb']) & (new['telomere_end_confirmed'].str.lower() == 'false'))
    if r is not None:
        rep += ['### Unchanged read that stops before the telomere', '',
                f'`{r["sample"]}` {r["chr_end"]} read `{r["read_id"]}` -- Y\' array `{r.get("y_prime_observed_array", "")}` '
                'matches the reference as far as the read goes.', '',
                '| | old | new |', '|---|---|---|',
                f'| score | `overall_confidence` = **{r["overall_confidence"]}** | '
                f'`recombination_confidence` = **{r["recombination_confidence"]}** |',
                f'| why | `{r.get("confidence_factors", "")}` | `{r.get("confidence_basis", "")}` |', '']

    path = os.path.join(out, 'confidence_old_vs_new.md')
    with open(path, 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
    print('\n'.join(rep))
    print(f'\nWritten: {path}')


if __name__ == '__main__':
    main()
