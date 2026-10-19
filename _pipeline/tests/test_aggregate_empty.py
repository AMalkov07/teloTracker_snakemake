"""Unit tests for the per-end recombination summary (aggregate_recombination.summarize).
Run: python _pipeline/tests/test_aggregate_empty.py

An end whose features file is empty (analyze_features found no reads to report) must appear in the
summary as 'no_events', not vanish from it.
"""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import pandas as pd
import aggregate_recombination as agg
from recombination_utils import write_results_tsv


def summarize(rows_by_end, sidecars=()):
    with tempfile.TemporaryDirectory() as d:
        for end, rows in rows_by_end.items():
            write_results_tsv(rows, os.path.join(d, f'S_{end}_features.tsv'))
        for end in sidecars:
            open(os.path.join(d, f'S_{end}_features.tsv.skipped'), 'w').write('insufficient_coverage\n')
        out = os.path.join(d, 'summary.tsv')
        agg.summarize(d, 'S', out)
        return pd.read_csv(out, sep='\t').set_index('chr_end')


def test_empty_features_file_is_no_events():
    row = {'read_id': 'r1', 'recombination_detected': False, 'y_prime_recombination_status': 'No Change'}
    s = summarize({'chr1L': [], 'chr2R': [row]})
    assert set(s.index) == {'chr1L', 'chr2R'}
    assert s.loc['chr1L', 'status'] == 'no_events' and s.loc['chr1L', 'total_reads'] == 0
    assert s.loc['chr2R', 'total_reads'] == 1


def test_skipped_sidecar_still_reported_as_skipped():
    s = summarize({'chr3L': []}, sidecars=['chr3L'])
    assert s.loc['chr3L', 'status'] == 'skipped'


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
