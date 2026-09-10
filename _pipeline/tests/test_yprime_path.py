"""Unit tests for the (ID, ITS) path parser. Run: python _pipeline/tests/test_yprime_path.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import yprime_path as yp

REF = {'chr4R': [('ID1', 10), ('ID1', 10), ('ID1', 10), ('ID1', 171), ('ID1', 172), ('ID1', 10), ('ID1', None)],
       'chr12R': [('ID5', 160), ('ID1', 170), ('ID1', 172), ('ID1', 172), ('ID1', 172), ('ID1', None)],
       'chr13L': [('ID2', 163), ('ID1', 152), ('ID2', 163), ('ID1', None)],
       'chr14L': [('ID1', 169), ('ID1', 172), ('ID2', 175), ('ID2', 175), ('ID2', None)],
       'chr2L': [('ID4', None)]}


def P(tokens, self_end='chr2L'):
    return yp.parse_path(tokens, REF, self_end)


def test_read_tokens_orientation():
    assert yp.read_tokens_from_positions('ID2:100-5600;ID1:5763-12400', 'end') == [('ID2', 163), ('ID1', None)]
    # telomere at the beginning: hits are in decreasing coordinate along anchor->telomere
    assert yp.read_tokens_from_positions('ID1:100-6700;ID2:6863-12300', 'beginning') == [('ID2', 163), ('ID1', None)]


def test_its_separates_identical_id_runs():
    p = P([('ID1', 10), ('ID1', 10), ('ID1', None)])
    assert p['primary_donor'] == 'chr4R' and p['n_segments'] == 1
    p = P([('ID1', 170), ('ID1', 172), ('ID1', None)])
    assert p['segments'][0]['tag'] == 'ambiguous' and set(p['segments'][0]['donor'].split('|')) == {'chr12R', 'chr4R'}
    p = P([('ID1', 170), ('ID1', 172), ('ID1', 172), ('ID1', None)])
    assert p['primary_donor'] == 'chr12R'                     # a third 172-bp ITS rules out chr4R (10 bp there)
    p = P([('ID1', None)])
    assert p['primary_donor'] == '' and p['segments'][0]['tag'] == 'ambiguous'


def test_alternating_and_circle():
    p = P([('ID2', 163), ('ID1', 152), ('ID2', 163), ('ID1', None)])
    assert p['primary_donor'] == 'chr13L' and p['its_verified'] == 3
    p = P([('ID2', 163), ('ID1', 40), ('ID2', 163), ('ID1', 37), ('ID2', 163), ('ID1', None)])
    assert p['primary_donor'] == 'chr13L' and p['circles'] and p['circles'][0]['unit'] == ['ID2', 'ID1']
    assert abs(p['circles'][0]['repeats'] - 3.0) < 1e-9
    assert yp.format_circles(p).startswith('chr13L[1-2]:ID2,ID1x3.0:')
    assert p['circles'][0]['circle_support'] == 'strong'       # 163 inside every unit: not a whole-array copy


def test_composite_path():
    p = P([('ID2', 163), ('ID1', 152), ('ID2', 163), ('ID1', 10), ('ID1', 10), ('ID1', None)])
    assert [s['donor'] for s in p['segments']] == ['chr13L', 'chr4R']
    assert yp.format_path(p) == 'chr13L[1-4]:ID2,ID1,ID2,ID1 > chr4R[1-2]:ID1,ID1'


def test_self_is_labelled():
    p = P([('ID2', 163), ('ID1', 152), ('ID2', None)], self_end='chr13L')
    assert p['primary_donor'] == 'self' and p['segments'][0]['tag'] == 'self'


def test_circle_junction_rules():
    # consistent junction ITS equal to the donor's flanking ITS -> chr4R, but a one-copy circle is tentative
    p = P([('ID1', 10), ('ID1', 10), ('ID1', 10), ('ID1', 10), ('ID1', 10), ('ID1', None)])
    assert p['circles'] and p['segments'][0]['donor'] == 'chr4R' and p['segments'][0]['tag'] == 'tentative'
    assert p['primary_donor'] == ''                            # a one-copy circle never names the donor outright
    # inconsistent junctions break the repeat
    p = P([('ID2', 163), ('ID1', 40), ('ID2', 163), ('ID1', 120), ('ID2', 163), ('ID1', None)])
    assert p['n_segments'] >= 2
    # a run of one ID with a 172-bp junction: chr4R, chr12R and chr14L all carry ID1 with a
    # ~172-bp ITS, so a single-copy circle stays ambiguous among them (and chr5R etc. drop out)
    p = P([('ID1', 172), ('ID1', 172), ('ID1', 172), ('ID1', 172), ('ID1', 172), ('ID1', 172), ('ID1', None)])
    assert set(p['segments'][0]['donor'].split('|')) == {'chr12R', 'chr14L', 'chr4R'}, yp.format_path(p)


def test_circle_any_phase():
    # circle of chr13L copies 1-2 (ID2 -163- ID1) inserted starting on its ID1, junction ~152
    p = P([('ID1', 152), ('ID2', 163), ('ID1', 152), ('ID2', 163), ('ID1', None)])
    assert p['primary_donor'] == 'chr13L' and p['circles'], yp.format_path(p)
    assert p['circles'][0]['unit'] == ['ID1', 'ID2'] and p['n_segments'] == 1
    # same but with a junction that matches nothing in the donor (new sequence): still one circle segment
    p = P([('ID1', 40), ('ID2', 163), ('ID1', 40), ('ID2', 163), ('ID1', None)])
    assert p['n_segments'] == 1 and p['circles'] and p['primary_donor'] == 'chr13L', yp.format_path(p)


def test_circle_support_levels():
    # 5 alternating copies: no single linear copy of the 4-copy chr13L array can give 5 -> strong
    p = P([('ID2', 163), ('ID1', 152), ('ID2', 163), ('ID1', 152), ('ID2', None)])
    assert p['circles'] and p['circles'][0]['circle_support'] == 'strong', yp.format_path(p)
    # 4 alternating copies in the donor's own phase = one verbatim copy: the circle reading is weak
    p = P([('ID2', 163), ('ID1', 152), ('ID2', 163), ('ID1', None)])
    assert (not p['circles']) or p['circles'][0]['circle_support'] == 'weak', yp.format_path(p)
    # 4 alternating copies starting on ID1 = copies 2-4 then copy 1: not one linear piece -> moderate, alt shown
    p = P([('ID1', 152), ('ID2', 163), ('ID1', 152), ('ID2', None)])
    assert p['circles'] and p['circles'][0]['circle_support'] == 'moderate' and 'chr13L[2-4]' in p['circles'][0]['alt'], yp.format_path(p)
    # 6 alternating copies with 163 inside every unit and a junction unlike any chr13L ITS: only a circle explains it
    p = P([('ID2', 163), ('ID1', 130), ('ID2', 163), ('ID1', 130), ('ID2', 163), ('ID1', None)])   # 130-bp junctions: unlike any chr13L ITS
    assert p['circles'] and p['circles'][0]['circle_support'] == 'strong', yp.format_path(p)


def test_its_tolerance():
    assert P([('ID1', 15), ('ID1', 6), ('ID1', None)])['primary_donor'] == 'chr4R'      # within tolerance
    assert P([('ID1', 30), ('ID1', 30), ('ID1', None)])['primary_donor'] != 'chr4R'    # outside


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
