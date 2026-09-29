"""Unit tests for the X-element call (call_x_element). Run: python _pipeline/tests/test_x_element_call.py

A switch to another end's X needs a hit covering half of that X and a clear identity lead
over the read's hit to its own X.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import analyze_features as af


def hit(cluster, end, pident, length, subject_len=720, qstart=10000):
    return {'pident': pident, 'length': length, 'qstart': qstart, 'qend': qstart + length,
            'bitscore': 1.8 * length * pident / 100, 'subject_len': subject_len,
            'cluster_id': cluster, 'members': [end], 'rep_chr_end': end}


def call(hits):
    return af.call_x_element(hits, 'chr10R', 'ID4')


def test_own_x_is_no_change():
    r = call([hit('ID4', 'chr10R', 99.5, 715), hit('ID19', 'chr4L', 95.0, 700)])
    assert r['x_element_recombination'] == 'no_change' and r['x_element_source'] == 'chr10R'


def test_clear_full_length_switch_is_kept():
    r = call([hit('ID19', 'chr4L', 99.0, 712), hit('ID4', 'chr10R', 95.5, 700)])
    assert r['x_element_recombination'] == 'full_switch' and r['x_element_source'] == 'chr4L'


def test_short_fragment_is_no_data():
    # 70 bp of a 713 bp X right before the telomere, best matching another end
    r = call([hit('ID19', 'chr4L', 89.0, 70, subject_len=713)])
    assert r['x_element_recombination'] == 'no_data'
    assert r['x_element_cluster_id'] == 'ID19'


def test_short_library_x_can_still_switch():
    # chr5L's X is 110 bp: a 105 bp hit covers it
    r = call([hit('ID21', 'chr5L', 99.0, 105, subject_len=110)])
    assert r['x_element_recombination'] == 'full_switch' and r['x_element_source'] == 'chr5L'


def test_near_tie_with_own_x_is_no_change():
    r = call([hit('ID19', 'chr4L', 97.2, 740), hit('ID4', 'chr10R', 96.6, 700)])
    assert r['x_element_recombination'] == 'no_change' and r['x_element_source'] == 'chr10R'


def test_own_x_higher_identity_but_lower_bitscore_is_no_change():
    r = call([hit('ID19', 'chr4L', 87.3, 760), hit('ID4', 'chr10R', 98.5, 600)])
    assert r['x_element_recombination'] == 'no_change'


def test_no_own_hit_keeps_the_switch():
    r = call([hit('ID19', 'chr4L', 93.0, 700)])
    assert r['x_element_recombination'] == 'full_switch'


def test_no_hits_is_no_data():
    assert call([])['x_element_recombination'] == 'no_data'


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
