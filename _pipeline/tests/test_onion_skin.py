"""Unit tests for the onion-skin summary. Run: python _pipeline/tests/test_onion_skin.py

Paths below are real y_prime_path strings from the 7302 day-5 and 7372 runs.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import onion_skin_summary as oss

G = "Y' Gain"


def row(path, circles='', gained='ID2,ID1', status=G, donor='', positions='', read_id='r'):
    return {'y_prime_recombination_status': status, 'y_prime_path': path,
            'y_prime_path_circles': circles, 'y_prime_gained_segment': gained,
            'y_prime_path_primary_donor': donor, 'y_prime_positions': positions,
            'telo_side': 'end', 'read_id': read_id, 'y_prime_divergence_idx': '1'}


def positions(ids_its):
    """[('ID1', 166), ('ID2', None)] -> 'ID1:0-6000;ID2:6166-12166' (telo_side 'end')."""
    out, x = [], 0
    for yid, its in ids_its:
        out.append(f'{yid}:{x}-{x + 6000}')
        x += 6000 + (its or 0)
    return ';'.join(out)


def test_parse_path_donors_and_circles():
    p = oss.parse_path('chr13L[1-3]:ID2,ID1,ID2,ID2,ID1(circ x1.7 strong)')
    assert p == [('chr13L', True)]
    p = oss.parse_path('chr14L?[3]:ID2 > chr12R[2-5]:ID1,ID1,ID1,ID1')
    assert p == [(None, False), ('chr12R', False)]          # tentative piece names no donor
    p = oss.parse_path('chr13L[1-3]:ID2,ID1,ID2 > chr12R|chr4R:ID1,ID1,ID1')
    assert p == [('chr13L', False), (None, False)]          # ambiguous piece names no donor
    p = oss.parse_path('chr13L[2-3]:ID2,ID1,ID2,ID1(circ x2.0 moderate | alt chr13L[1-3] + chr13L[2])')
    assert p == [('chr13L', True)]


def test_segments_carry_support_and_length():
    s = oss.parse_segments('self[1]:ID2 > self[3]:ID1,ID1,ID1,ID1(circ x4.0 weak | alt chr14L[3-6]) > self[1]:ID2')
    assert [(x['donor'], x['n_ids'], x['support']) for x in s] == [('self', 1, ''), ('self', 4, 'weak'), ('self', 1, '')]


def test_named_repeat_by_circle_or_repeated_donor():
    s = oss.summarise([
        row('chr13L[1-3]:ID2,ID1,ID2,ID2,ID1(circ x1.7 strong)', 'chr13L[1-3]:ID2,ID1,ID2x1.7:strong'),
        row('chr14L[1-3]:ID1,ID1,ID2 > chr14L[3-4]:ID2,ID2'),          # same donor twice, no circle
        row('chr13L[1-4]:ID2,ID1,ID2,ID1'),                             # one piece, not a repeat
        row('chr14L?[3]:ID2 > chr14L?[3]:ID2'),                         # tentative: never a donor
    ])
    assert s['n_named_repeat'] == 2 and s['n_unconfirmed_repeat'] == 0
    assert s['n_circle'] == 1 and s['n_circle_strong'] == 1
    assert s['n_single_donor'] == 2 and s['n_multi_donor'] == 2


def test_ambiguous_donor_circle_is_unconfirmed_not_named():
    # 7372 day 3 chr2R: three ID2 whose ITS fit three ends -- the old code counted it as a
    # same-donor repeat and as a circle
    s = oss.summarise([row('chr14L|chr4R:ID2,ID2,ID2(circ x1.5)', 'chr14L|chr4R:ID2,ID2x1.5:na')])
    assert s['n_named_repeat'] == 0 and s['n_unconfirmed_repeat'] == 1
    assert s['n_circle'] == 0 and s['n_circle_unassigned'] == 1


def test_weak_circle_alone_is_unconfirmed():
    s = oss.summarise([row('self[3]:ID1,ID1,ID1,ID1(circ x4.0 weak | alt chr14L[3-6])')])
    assert s['n_named_repeat'] == 0 and s['n_unconfirmed_repeat'] == 1
    assert s['n_circle'] == 1 and s['n_circle_weak'] == 1


def test_tentative_circle_counted_as_unassigned():
    s = oss.summarise([row('chr13L?[1]:ID2,ID2,ID2(circ x3.0 moderate)')])
    assert s['n_circle'] == 0 and s['n_circle_unassigned'] == 1 and s['n_unconfirmed_repeat'] == 1


def test_only_gain_like_statuses_count():
    s = oss.summarise([row('chr13L[1-2]:ID2,ID1', status='No Change'),
                       row('chr13L[1-2]:ID2,ID1', status="Y' Loss"),
                       row('chr13L[1-2]:ID2,ID1', status="1st Y' Change")])
    assert s['n_reads'] == 3 and s['n_gain_like'] == 1


def test_clonal_reads_collapse_into_one_event():
    # 7372 day 5 chr13L: three reads of the same five-copy gain, ITS within 4 bp
    path = 'self[2]:ID2,ID2,ID2,ID2,ID2(circ x5.0 strong)'
    reads = [row(path, positions=positions(t), read_id=f'r{i}') for i, t in enumerate([
        [('ID1', 166), ('ID2', 166), ('ID2', 167), ('ID2', 166), ('ID2', 167), ('ID2', 166), ('ID2', None)],
        [('ID1', 166), ('ID2', 166), ('ID2', 167), ('ID2', 166), ('ID2', 166), ('ID2', 166), ('ID2', None)],
        [('ID1', 164), ('ID2', 166), ('ID2', 170), ('ID2', 167), ('ID2', 166), ('ID2', 166), ('ID2', None)]])]
    other = row('self[2]:ID2,ID2,ID2,ID2(circ x4.0 strong)', read_id='r9',
                positions=positions([('ID1', 166), ('ID2', 166), ('ID2', 167), ('ID2', 166), ('ID2', 166), ('ID2', None)]))
    events = oss.collapse_events(reads + [other])
    assert sorted(len(e['members']) for e in events) == [1, 3]
    desc = [oss.describe_event(e, 'chr13L', k + 1) for k, e in enumerate(events)]
    s = oss.summarise(reads + [other], desc)
    assert s['n_gain_like'] == 4 and s['n_gain_events'] == 2
    assert s['n_named_repeat'] == 4 and s['n_named_repeat_events'] == 2


def test_different_its_is_a_different_event():
    a = [('ID1', 166), ('ID2', 166), ('ID2', None)]
    b = [('ID1', 166), ('ID2', 10), ('ID2', None)]          # a 10 bp junction: another event
    assert oss.same_array(a, a) and not oss.same_array(a, b)
    assert oss.format_array(a, 2) == 'ID1(166) ID2(166) | ID2'


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
