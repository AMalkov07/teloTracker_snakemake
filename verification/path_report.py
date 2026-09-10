#!/usr/bin/env python3
"""path_report.py <snapshot> <out_dir> <sample>... -- what the (ID, ITS) path adds.
Per sample: how many gain-like reads get a unique donor from the ID fingerprint alone vs
from the path; the primary donors; and circle segments (unit x repeats) per end/donor."""
import os, sys
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_pipeline', 'scripts'))
import pandas as pd
from verify_recombination import load_features, end_sort_key, md_table
GAIN = ("Y' Gain", "1st Y' Change", "Y' Recombination")

def main():
    snap, out, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
    os.makedirs(out, exist_ok=True); rep = [f'# Y\' path report -- {snap}\n']
    for s in samples:
        df = load_features(snap, s)
        if df.empty or 'y_prime_path' not in df.columns:
            rep.append(f'## {s}\n_(no path columns)_\n'); continue
        g = df[df['y_prime_recombination_status'].isin(GAIN)].copy()
        g['fp'] = g['y_prime_fingerprint_source'].fillna('').astype(str)
        g['pd'] = g['y_prime_path_primary_donor'].fillna('').astype(str)
        g['self_fp'] = g['y_prime_self_match'].astype(str) == 'True'
        g['nseg'] = pd.to_numeric(g['y_prime_path_n_segments'], errors='coerce').fillna(0).astype(int)
        g['ngained'] = g['y_prime_gained_segment'].fillna('').map(lambda x: len([i for i in x.split(',') if i]))
        multi = g[g['ngained'] >= 2]
        def cat(r):
            if r['self_fp'] or r['pd'] == 'self': return 'self (tandem amplification)'
            if r['fp']: return 'unique donor by IDs'
            if r['pd'] and '|' not in r['pd']: return 'unique donor by IDs+ITS'
            if r['pd'] and '|' in r['pd']: return 'ambiguous (2+ donors)'
            return 'unresolved'
        multi = multi.assign(cat=multi.apply(cat, axis=1))
        cats = multi['cat'].value_counts()
        donors = Counter(multi[multi['cat'].str.startswith('unique')].apply(lambda r: r['fp'] or r['pd'], axis=1))
        # circles from the path: donor[piece]:unit x reps:support
        rows = []
        for _, r in g.iterrows():
            for c in str(r['y_prime_path_circles'] or '').split(';'):
                if c.count(':') != 2: continue
                donor, unit_reps, support = c.split(':')
                unit, reps = unit_reps.rsplit('x', 1)
                rows.append({'chr_end': r['chr_end'], 'donor': donor, 'unit': unit, 'repeats': float(reps),
                             'support': support, 'read_id': r['read_id'], 'path': r['y_prime_path'], 'source': r['recombination_source']})
        circ = pd.DataFrame(rows)
        circ.to_csv(os.path.join(out, f'path_circles_{s}.tsv'), sep='\t', index=False)
        multi[['read_id', 'chr_end', 'y_prime_observed_array', 'y_prime_gained_segment', 'y_prime_path', 'cat', 'recombination_source', 'overall_confidence']].to_csv(os.path.join(out, f'path_reads_{s}.tsv'), sep='\t', index=False)
        by = (circ.groupby(['donor', 'unit']).agg(n_reads=('read_id', 'size'), n_ends=('chr_end', 'nunique'), max_repeats=('repeats', 'max'),
                                                 n_strong=('support', lambda x: int((x == 'strong').sum())), n_weak=('support', lambda x: int((x == 'weak').sum())),
                                                 ends=('chr_end', lambda x: ','.join(sorted(set(x), key=end_sort_key)))).reset_index().sort_values('n_reads', ascending=False)) if not circ.empty else pd.DataFrame()
        comp = g[g['nseg'] >= 2]
        rep += [f'## {s}', f'{len(g)} gain-like reads, {len(multi)} with >= 2 gained Y\'.',
                '### Donor resolution of multi-Y\' gains', md_table(cats.rename_axis('category').reset_index(name='n_reads')),
                '### Donors (unique calls)', md_table(pd.DataFrame(donors.most_common(), columns=['donor', 'n_reads'])),
                f'### Circles from the path ({len(circ)} reads): donor x repeat unit', md_table(by, 30),
                f'### Composite paths (>= 2 segments): {len(comp)} reads; examples', md_table(comp[['chr_end', 'y_prime_path']].head(12)), '']
        print(f"{s}: {len(multi)} multi-Y' gains -> {dict(cats)}; circles {len(circ)} reads; top donors {donors.most_common(3)}")
    open(os.path.join(out, 'path_report.md'), 'w').write('\n'.join(rep))
main()
