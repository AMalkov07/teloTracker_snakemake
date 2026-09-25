#!/usr/bin/env python3
"""
confidence_audit.py -- what does the recombination confidence score actually measure?

Read-only audit of the CURRENT overall_confidence / mean_confidence, run before redesigning
them. Four questions, each answered from data that already carries ground truth:

  1. Distribution by event class. Which values does the score take for each kind of read, and
     how many reads sit on a hardcoded value (0.95 no-recombination, 0.90 self-source-only,
     0.30 / 0.33 Y'-only floor)?
  2. Is the per-end mean_confidence just a restatement of the recombination rate? Fit
     mean ~ a*(1-p) + b*p across every (sample, end) and report how much it explains.
  3. Does confidence separate CORRECT from WRONG donor calls? Two truth sets:
       - positive control: 7172 chr11L reads, true donor chr11R (verified by direct BLAST)
       - chr13L truth set: reads at other ends that gained chr13L's unique alternating
         array (B5 table from Part B verification), true donor chr13L
     Reported as AUC (0.5 = no separation, 1.0 = perfect).
  4. Does confidence separate DAY-0 calls (at best standing variation) from calls in later
     timepoints? If it measured "did recombination happen", day-0 calls should score lower.

Usage:
  python verification/confidence_audit.py <snapshot_dir> <truth_set_tsv> <out_dir>
  e.g. verification/snapshot_v2c_path verification/reports/v2c_path/B5_truth_set_chr13L.tsv \
       verification/reports/confidence_audit
"""
import glob
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd

NO_CHANGE = {'no_change', '', 'nan', 'none', 'None'}
YP_NO_CHANGE = {'No Change', '', 'nan'}


def load(snapshot):
    frames = []
    for d in sorted(glob.glob(os.path.join(snapshot, '*'))):
        s = os.path.basename(d)
        for f in glob.glob(os.path.join(d, 'recombination', f'{s}_chr*_features.tsv')):
            try:
                df = pd.read_csv(f, sep='\t', dtype=str, keep_default_na=False)
            except pd.errors.EmptyDataError:
                continue
            if df.empty or 'overall_confidence' not in df:
                continue
            df['sample'] = s
            frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df['conf'] = pd.to_numeric(df['overall_confidence'], errors='coerce')
    df['recomb'] = df['recombination_detected'].astype(str).str.lower() == 'true'
    return df


def summaries(snapshot):
    rows = []
    for d in sorted(glob.glob(os.path.join(snapshot, '*'))):
        s = os.path.basename(d)
        f = os.path.join(d, 'recombination', f'{s}_recombination_summary.tsv')
        if not os.path.isfile(f):
            continue
        t = pd.read_csv(f, sep='\t', dtype=str, keep_default_na=False)
        t = t[t['status'] == 'analyzed'].copy()      # columns read BY NAME, never position
        t['sample'] = s
        rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    for c in ('total_reads', 'pct_recombination', 'mean_confidence'):
        t[c] = pd.to_numeric(t[c], errors='coerce')
    return t


def event_class(r):
    if not r['recomb']:
        return 'no recombination'
    sp = str(r.get('spacer_recombination', '')) not in NO_CHANGE
    xe = str(r.get('x_element_recombination', '')) not in NO_CHANGE
    yp = str(r.get('y_prime_recombination_status', '')) not in YP_NO_CHANGE
    axes = [n for n, on in (('spacer', sp), ('X', xe), ("Y'", yp)) if on]
    name = '+'.join(axes) if axes else 'none of the three axes'
    if str(r.get('is_complex_event', '')).lower() == 'true':
        name += ' (complex)'
    return name


def auc(pos, neg):
    """Mann-Whitney AUC: P(score of a positive > score of a negative), ties = 1/2."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if not len(pos) or not len(neg):
        return float('nan')
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank(method='average').to_numpy()
    rp = ranks[:len(pos)].sum()
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def md(df, floatfmt=3):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join(f'{v:.{floatfmt}f}' if isinstance(v, float) else str(v) for v in r) + ' |')
    return '\n'.join(out)


def main():
    snapshot, truth_tsv, out = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    df = load(snapshot)
    rep = ['# Recombination confidence audit', '',
           f'Snapshot: `{snapshot}` -- {df["sample"].nunique()} samples, {len(df):,} reads. '
           'Read-only: this measures the current `overall_confidence`; nothing is changed.', '']

    # ---- 1. distribution by event class ------------------------------------------------
    df['event_class'] = df.apply(event_class, axis=1)
    g = df.groupby('event_class')['conf']
    t1 = pd.DataFrame({'reads': g.size(), 'mean': g.mean(), 'median': g.median(),
                       'distinct values': g.nunique(),
                       '% at 0.95': g.apply(lambda s: 100 * (s.round(4) == 0.95).mean()),
                       '% at 0.30/0.33': g.apply(lambda s: 100 * s.round(4).isin([0.3, 0.33]).mean())})
    t1 = t1.sort_values('reads', ascending=False).reset_index()
    rec = df[df['recomb']]
    floor = 100 * rec['conf'].round(4).isin([0.3, 0.33]).mean()
    rep += ['## 1. What values the score takes, by kind of read', '',
            f'Of {len(rec):,} recombinant reads, **{floor:.1f}% sit exactly on the 0.30 / 0.33 floor** '
            '(base 0.3, x1.1 when the Y\' array agrees). Every non-recombinant read is a constant.', '',
            md(t1, 3), '']
    top_factors = Counter(rec['confidence_factors']).most_common(8)
    rep += ['Most common `confidence_factors` among recombinant reads:', '']
    rep += [f'- `{k}`: {v:,}' for k, v in top_factors] + ['']

    # ---- 2. per-end mean vs recombination rate -----------------------------------------
    t = summaries(snapshot).dropna(subset=['pct_recombination', 'mean_confidence'])
    p = t['pct_recombination'] / 100.0
    X = np.column_stack([1 - p, p])
    coef, *_ = np.linalg.lstsq(X, t['mean_confidence'], rcond=None)
    pred = X @ coef
    ss_res = ((t['mean_confidence'] - pred) ** 2).sum()
    ss_tot = ((t['mean_confidence'] - t['mean_confidence'].mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot
    r = np.corrcoef(p, t['mean_confidence'])[0, 1]
    rep += ['## 2. Is the per-end mean just the recombination rate restated?', '',
            f'Across {len(t):,} (sample, end) pairs, fitting `mean_confidence = a*(1-p) + b*p` '
            f'(p = fraction recombinant) gives a = {coef[0]:.3f}, b = {coef[1]:.3f}, '
            f'**R^2 = {r2:.3f}** (Pearson r = {r:.3f}). a ~ 0.95 is the constant given to every '
            'non-recombinant read. An R^2 this high means the per-end mean carries almost no '
            'information beyond the recombination rate already reported next to it.', '']

    # ---- 3. correct vs wrong donor -----------------------------------------------------
    pc = df[df['sample'].isin(['7172_day4_with_selection', '7172_day6_with_selection',
                               '7172_day9_with_selection']) & (df['chr_end'] == 'chr11L') & df['recomb']]
    pc_ok = pc[pc['recombination_source'] == 'chr11R']['conf']
    pc_bad = pc[pc['recombination_source'] != 'chr11R']['conf']
    ts = pd.read_csv(truth_tsv, sep='\t', dtype=str, keep_default_na=False)
    ts['conf'] = pd.to_numeric(ts['overall_confidence'], errors='coerce')
    ts_ok = ts[ts['hit'].str.lower() == 'true']['conf']
    ts_bad = ts[ts['hit'].str.lower() != 'true']['conf']
    ok, bad = pd.concat([pc_ok, ts_ok]), pd.concat([pc_bad, ts_bad])
    t3 = pd.DataFrame([
        ['positive control (7172 chr11L -> chr11R)', len(pc_ok), len(pc_bad),
         pc_ok.mean(), pc_bad.mean(), auc(pc_ok, pc_bad)],
        ['chr13L truth set (7302)', len(ts_ok), len(ts_bad), ts_ok.mean(), ts_bad.mean(), auc(ts_ok, ts_bad)],
        ['both pooled', len(ok), len(bad), ok.mean(), bad.mean(), auc(ok, bad)],
    ], columns=['truth set', 'correct donor', 'wrong donor', 'mean conf (correct)',
                'mean conf (wrong)', 'AUC'])
    rep += ['## 3. Does confidence separate correct from wrong donor calls?', '',
            'AUC = probability a correctly attributed read scores higher than a wrongly attributed '
            'one. 0.5 = no information, 1.0 = perfect.', '', md(t3, 3), '']
    wrong_sources = Counter(pd.concat([pc[pc['recombination_source'] != 'chr11R']['recombination_source'],
                                       ts[ts['hit'].str.lower() != 'true']['recombination_source']]))
    rep += ['Wrong calls went to: ' + ', '.join(f'`{k or "(blank)"}` {v}' for k, v in wrong_sources.most_common(6)), '']

    # ---- 4. day-0 vs later timepoints --------------------------------------------------
    day0 = rec[rec['sample'].str.contains('_day0_')]['conf']
    later = rec[~rec['sample'].str.contains('_day0_') & ~rec['sample'].str.contains('survivor')]['conf']
    surv = rec[rec['sample'].str.contains('survivor')]['conf']
    t4 = pd.DataFrame([['day-0 self-runs', len(day0), day0.mean(), day0.median()],
                       ['later timepoints', len(later), later.mean(), later.median()],
                       ['survivors', len(surv), surv.mean(), surv.median()]],
                      columns=['recombinant reads in', 'n', 'mean conf', 'median conf'])
    rep += ['## 4. Does confidence separate day-0 calls from real events?', '',
            'At day 0 there should be no recombination, so day-0 recombinant calls are at best '
            'standing variation and at worst false. If the score measured "did recombination '
            'happen", they should score clearly lower than calls at later timepoints.', '',
            md(t4, 3), '',
            f'**AUC (later timepoint vs day-0) = {auc(later, day0):.3f}.**', '']

    path = os.path.join(out, 'confidence_audit.md')
    with open(path, 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
    t1.to_csv(os.path.join(out, 'confidence_by_event_class.tsv'), sep='\t', index=False)
    print('\n'.join(rep))
    print(f'\nWritten: {path}')


if __name__ == '__main__':
    main()
