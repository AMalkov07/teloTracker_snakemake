"""Unit tests for the onion-skin summary. Run: python _pipeline/tests/test_onion_skin.py

Paths below are real y_prime_path strings from the 7302 day-5 run.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import onion_skin_summary as oss

G = "Y' Gain"


def row(path, circles='', gained='ID2,ID1', status=G, donor=''):
    return {'y_prime_recombination_status': status, 'y_prime_path': path,
            'y_prime_path_circles': circles, 'y_prime_gained_segment': gained,
            'y_prime_path_primary_donor': donor}


def test_parse_path_donors_and_circles():
    p = oss.parse_path('chr13L[1-3]:ID2,ID1,ID2,ID2,ID1(circ x1.7 strong)')
    assert p == [('chr13L', True)]
    p = oss.parse_path('chr14L?[3]:ID2 > chr12R[2-5]:ID1,ID1,ID1,ID1')
    assert p == [(None, False), ('chr12R', False)]          # tentative piece names no donor
    p = oss.parse_path('chr13L[1-3]:ID2,ID1,ID2 > chr12R|chr4R:ID1,ID1,ID1')
    assert p == [('chr13L', False), (None, False)]          # ambiguous piece names no donor
    p = oss.parse_path('chr13L[2-3]:ID2,ID1,ID2,ID1(circ x2.0 moderate | alt chr13L[1-3] + chr13L[2])')
    assert p == [('chr13L', True)]


def test_same_donor_repeat_by_circle_or_repeated_donor():
    s = oss.summarise([
        row('chr13L[1-3]:ID2,ID1,ID2,ID2,ID1(circ x1.7 strong)', 'chr13L[1-3]:ID2,ID1,ID2x1.7:strong'),
        row('chr14L[1-3]:ID1,ID1,ID2 > chr14L[3-4]:ID2,ID2'),          # same donor twice, no circle
        row('chr13L[1-4]:ID2,ID1,ID2,ID1'),                             # one piece, not a repeat
        row('chr14L?[3]:ID2 > chr14L?[3]:ID2'),                         # tentative: never a donor
    ])
    assert s['n_same_donor_repeat'] == 2
    assert s['n_circle'] == 1 and s['n_circle_strong'] == 1
    assert s['n_single_donor'] == 2 and s['n_multi_donor'] == 2


def test_unassigned_circle_counted_separately():
    s = oss.summarise([row('chr13L?[1]:ID2,ID2,ID2(circ x3.0 moderate)')])
    assert s['n_circle'] == 0 and s['n_circle_unassigned'] == 1


def test_only_gain_like_statuses_count():
    s = oss.summarise([row('chr13L[1-2]:ID2,ID1', status='No Change'),
                       row('chr13L[1-2]:ID2,ID1', status="Y' Loss"),
                       row('chr13L[1-2]:ID2,ID1', status="1st Y' Change")])
    assert s['n_reads'] == 3 and s['n_gain_like'] == 1


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
