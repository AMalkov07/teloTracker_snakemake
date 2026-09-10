#!/usr/bin/env python3
"""
diff_runs.py -- per-read comparison of two recombination snapshots (e.g. the
Argon v1 run vs the v2 reprocessing, or v2-reprocessed vs v2-rerun on Argon).

Usage: diff_runs.py <snapshot_a> <snapshot_b> <out_dir> [sample ...]
Writes <out_dir>/diff_<sample>.tsv (per-read rows whose status/source/detected
changed), <out_dir>/diff_summary.tsv (per sample/end counts) and
<out_dir>/diff_report.md.
"""
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
import pandas as pd
from verify_recombination import load_features, end_sort_key, md_table

KEY = ['read_id', 'chr_end']
COLS = ['y_prime_recombination_status', 'recombination_source', 'recombination_detected', 'overall_confidence',
        'is_complex_event', 'recombination_mechanism']


def main():
    a_dir, b_dir, out = sys.argv[1:4]
    samples = sys.argv[4:] or sorted(d for d in os.listdir(b_dir) if os.path.isdir(os.path.join(b_dir, d, 'recombination')))
    os.makedirs(out, exist_ok=True)
    summary, report = [], ['# Per-read diff: A = %s, B = %s\n' % (a_dir, b_dir)]
    for s in samples:
        A, B = load_features(a_dir, s), load_features(b_dir, s)
        if A.empty or B.empty:
            continue
        for df in (A, B):
            for c in COLS:
                if c not in df.columns:
                    df[c] = ''
        m = A[KEY + COLS].merge(B[KEY + COLS], on=KEY, how='inner', suffixes=('_a', '_b'))
        m['status_changed'] = m['y_prime_recombination_status_a'] != m['y_prime_recombination_status_b']
        m['source_changed'] = m['recombination_source_a'].astype(str) != m['recombination_source_b'].astype(str)
        m['detected_changed'] = m['recombination_detected_a'].astype(str) != m['recombination_detected_b'].astype(str)
        m['dconf'] = (pd.to_numeric(m['overall_confidence_b'], errors='coerce') - pd.to_numeric(m['overall_confidence_a'], errors='coerce')).round(3)
        changed = m[m['status_changed'] | m['source_changed'] | m['detected_changed']]
        changed.to_csv(os.path.join(out, f'diff_{s}.tsv'), sep='\t', index=False)
        transitions = Counter(zip(changed['recombination_source_a'].astype(str), changed['recombination_source_b'].astype(str)))
        for ce, g in m.groupby('chr_end'):
            summary.append({'sample': s, 'chr_end': ce, 'n': len(g),
                            'n_status_changed': int(g['status_changed'].sum()),
                            'n_source_changed': int(g['source_changed'].sum()),
                            'n_detected_changed': int(g['detected_changed'].sum()),
                            'n_detected_a': int((g['recombination_detected_a'].astype(str) == 'True').sum()),
                            'n_detected_b': int((g['recombination_detected_b'].astype(str) == 'True').sum()),
                            'n_ambiguous_a': int((g['recombination_source_a'].astype(str) == 'ambiguous').sum()),
                            'n_ambiguous_b': int((g['recombination_source_b'].astype(str) == 'ambiguous').sum()),
                            'mean_dconf': round(g['dconf'].mean(), 3)})
        report += [f'## {s}', f'reads compared: {len(m)}; status changed: {int(m["status_changed"].sum())}; '
                   f'source changed: {int(m["source_changed"].sum())}; detected changed: {int(m["detected_changed"].sum())}; '
                   f'ambiguous {int((m["recombination_source_a"].astype(str)=="ambiguous").sum())} -> '
                   f'{int((m["recombination_source_b"].astype(str)=="ambiguous").sum())}',
                   '', 'Top source transitions (A -> B):',
                   md_table(pd.DataFrame([{'from': k[0], 'to': k[1], 'n': v} for k, v in transitions.most_common(15)])), '']
        print(f'{s}: {len(m)} reads, status changed {int(m["status_changed"].sum())}, source changed {int(m["source_changed"].sum())}, '
              f'detected changed {int(m["detected_changed"].sum())}')
    pd.DataFrame(summary).to_csv(os.path.join(out, 'diff_summary.tsv'), sep='\t', index=False)
    with open(os.path.join(out, 'diff_report.md'), 'w') as fh:
        fh.write('\n'.join(report))


if __name__ == '__main__':
    main()
