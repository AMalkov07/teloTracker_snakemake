#!/usr/bin/env python3
"""
circle_report.py -- candidate circular / rolling-circle Y' events from v2 features
(ID-fingerprint level; the (ID, ITS) path report is the finer view).

A read is a circle candidate when the gained part of its Y' array repeats a unit
(>= 2 Y', >= 2 copies), is a run of one ID (>= 3 copies, weak), is explained by
another end only as a ROTATION, or repeats its own end's array.

Usage: circle_report.py <snapshot_dir> <out_dir> <sample> [sample ...]
"""
import os, sys
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
import pandas as pd
from verify_recombination import load_features, end_sort_key, md_table
from analyze_features import _periodic_unit

def main():
    snap, out = sys.argv[1], sys.argv[2]; samples = sys.argv[3:]
    os.makedirs(out, exist_ok=True); rep = [f'# Circle candidates -- snapshot {snap}\n']
    for s in samples:
        df = load_features(snap, s)
        if df.empty or 'y_prime_fingerprint_kind' not in df.columns:
            rep.append(f'## {s}\n_(no v2 columns)_\n'); continue
        g = df[df['y_prime_recombination_status'].isin(["Y' Gain", "1st Y' Change", "Y' Recombination"])].copy()
        g['kind'] = g['y_prime_fingerprint_kind'].fillna('')
        g['donor'] = g['y_prime_fingerprint_source'].fillna('')
        g['self'] = g['y_prime_self_match'].astype(str) == 'True'
        def repeat_info(seg):
            ids = [x for x in str(seg or '').split(',') if x]
            if len(ids) < 2: return ('', 0)
            if len(set(ids)) == 1: return (ids[0], len(ids))
            u = _periodic_unit(ids)
            if u and len(set(u)) > 1: return (','.join(u), len(ids) / len(u))
            best = ('', 0)
            for L in range(len(ids), 3, -1):
                for i in range(0, len(ids) - L + 1):
                    w = ids[i:i + L]; u = _periodic_unit(w)
                    if u and len(set(w)) > 1 and len(w) / len(u) >= 2 and L > best[1] * max(len(best[0].split(',')), 1):
                        best = (','.join(u), len(w) / len(u))
                if best[0]: break
            return best
        ri = g['y_prime_gained_segment'].map(repeat_info)
        g['unit'] = ri.map(lambda t: t[0]); g['copies'] = ri.map(lambda t: round(t[1], 1))
        g['unit_len'] = g['unit'].map(lambda u: len(u.split(',')) if u else 0)
        cand = g[((g['unit_len'] >= 2) & (g['copies'] >= 2)) | (g['kind'] == 'rotation') |
                 ((g['unit_len'] == 1) & (g['copies'] >= 3))].copy()
        def cls(r):
            if r['unit_len'] == 1: return 'homopolymer_run'
            if r['self']: return 'tandem_amplification_own_array'
            if r['kind'] == 'rotation': return 'circle_reinserted_other_phase'
            return 'repeated_unit_from_donor' if r['donor'] else 'repeated_unit_donor_ambiguous'
        cand['circle_class'] = cand.apply(cls, axis=1)
        cols = ['read_id', 'chr_end', 'y_prime_recombination_status', 'y_prime_observed_array', 'y_prime_gained_segment',
                'y_prime_fingerprint_window', 'unit', 'copies', 'kind', 'donor', 'y_prime_array_matches', 'recombination_source',
                'overall_confidence', 'recombination_mechanism', 'telomere_end_confirmed', 'circle_class']
        cand[cols].to_csv(os.path.join(out, f'circles_{s}.tsv'), sep='\t', index=False)
        by = cand.groupby(['chr_end', 'circle_class', 'donor', 'unit']).agg(n_reads=('read_id', 'size'), max_copies=('copies', 'max')).reset_index()
        by = by.sort_values(['n_reads'], ascending=False)
        per_end = df.groupby('chr_end').size()
        by['pct_of_end'] = by.apply(lambda r: round(100 * r['n_reads'] / per_end.get(r['chr_end'], 1), 1), axis=1)
        rep += [f'## {s}', f'{len(g)} gain-like reads; **{len(cand)} circle candidates** '
                f'({dict(Counter(cand["circle_class"]))}); ends with >= 3 candidate reads: '
                f'{", ".join(sorted((cand.groupby("chr_end").size()[lambda x: x >= 3]).index, key=end_sort_key)) or "none"}\n',
                '### Recurrent (end, class, donor, repeat unit) with >= 2 reads', md_table(by[by['n_reads'] >= 2], 60),
                '### Donor ends of circle candidates', md_table(pd.DataFrame(Counter(cand[cand['donor'] != '']['donor']).most_common(), columns=['donor', 'n'])), '']
        print(f'{s}: {len(cand)} circle candidates of {len(g)} gain-like reads; classes {dict(Counter(cand["circle_class"]))}')
    open(os.path.join(out, 'circle_report.md'), 'w').write('\n'.join(rep))
main()
