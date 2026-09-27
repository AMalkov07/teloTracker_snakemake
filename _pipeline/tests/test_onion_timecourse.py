"""Unit tests for the onion-skin timecourse step. Run: python _pipeline/tests/test_onion_timecourse.py

The arrays mirror the 7372 chr13L ladder: reference ID1(166) ID2, then more and more copies of
the terminal ID2, joined by chr13L's own 166 bp ITS.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import onion_skin_timecourse as tc

REF = {'chr13L': [('ID1', 166), ('ID2', None)]}
SAMPLES = ['d0', 'd3', 'd5']


def ev(end, array, path, n=1, div=2):
    return {'chr_end': end, 'array': array, 'y_prime_path': path, 'n_reads': str(n),
            'divergence_idx': str(div), 'tier': '', 'repeat_donor': ''}


def ladder_events():
    return {
        'd0': [ev('chr13L', 'ID1(166) ID2(166) | ID2', 'self[2]:ID2', 2)],
        'd3': [ev('chr13L', 'ID1(166) ID2(166) | ID2', 'self[2]:ID2', 14),
               ev('chr13L', 'ID1(166) ID2(166) | ID2(167) ID2', 'chr14L|chr4R:ID2,ID2', 5)],
        'd5': [ev('chr13L', 'ID1(166) ID2(165) | ID2(166) ID2(167) ID2', 'self[2]:ID2,ID2,ID2(circ x3.0 strong)', 3),
               ev('chr13L', 'ID1(166) ID2(166) | ID2(10) ID2', 'chr4R[1-2]:ID2,ID2', 1)],
    }


def test_parse_array():
    assert tc.parse_array('ID1(166) ID2(166) | ID2') == [('ID1', 166), ('ID2', 166), ('ID2', None)]


def test_extends_needs_prefix_and_matching_its():
    a = [('ID1', 166), ('ID2', None)]
    assert tc.extends(a, [('ID1', 166), ('ID2', 166), ('ID2', None)]) == 2
    assert tc.extends(a, [('ID1', 190), ('ID2', 166), ('ID2', None)]) == -1   # ITS inside a differs
    assert tc.extends(a, [('ID2', 166), ('ID2', None)]) == -1                   # not a prefix
    assert tc.extends(a, a) == -1                                               # nothing added


def test_continues_unit_with_phase():
    assert tc.continues_unit(['ID2'], ['ID2', 'ID2'])
    assert tc.continues_unit(['ID2', 'ID1', 'ID2', 'ID1'], ['ID2', 'ID1'])
    assert tc.continues_unit(['ID1', 'ID2', 'ID1'], ['ID2'])      # period 2, next is ID2
    assert not tc.continues_unit(['ID2'], ['ID1'])


def test_pool_merges_same_array_across_samples():
    arrays = tc.pool_arrays(SAMPLES, ladder_events())
    first = [a for a in arrays if len(a['tokens']) == 3][0]
    assert first['reads'] == {'d0': 2, 'd3': 14} and first['first'] == 0
    assert len(arrays) == 4


def test_ladder_chain_is_time_ordered_repeats():
    arrays = tc.pool_arrays(SAMPLES, ladder_events())
    layers = tc.build_layers(arrays, REF)
    kinds = sorted((l['from_n_copies'], l['n_added'], l['layer_type']) for l in layers)
    # REF->+1 repeat; +1->+2 repeat (the parser called +2 'chr14L|chr4R', which must not matter);
    # +2->+3 repeat; +1->(10 bp junction) is the same unit but a new junction ITS
    assert (2, 1, 'repeat') in kinds and (3, 1, 'repeat') in kinds and (4, 1, 'repeat') in kinds
    assert (3, 1, 'repeat_its_mismatch') in kinds
    chain, longest = tc.repeat_chains('chr13L', layers)
    assert len(chain) == 3 and longest == 3
    assert tc.chain_time_ordered(chain)
    assert [SAMPLES[l['to_first_seen']] for l in chain] == ['d0', 'd3', 'd5']


def test_different_unit_is_not_a_repeat():
    events = {'d0': [], 'd3': [ev('chr13L', 'ID1(166) ID2(166) | ID5', 'chr14R[1]:ID5')], 'd5': []}
    layers = tc.build_layers(tc.pool_arrays(SAMPLES, events), REF)
    assert [l['layer_type'] for l in layers] == ['different']


def test_ladders_count_reference_plus_k_repeats():
    arrays = tc.pool_arrays(SAMPLES, ladder_events())
    n_reads = {('d0', 'chr13L'): 48, ('d3', 'chr13L'): 244, ('d5', 'chr13L'): 136}
    rows = {r['sample']: r for r in tc.ladders(arrays, REF, SAMPLES, n_reads)}
    assert rows['d3']['n_repeat_extended'] == 19 and rows['d3']['repeat_by_extra_copies'] == '+1:14;+2:5'
    assert rows['d5']['n_extended'] == 4 and rows['d5']['n_repeat_extended'] == 3   # the 10 bp one is not a repeat
    assert rows['d0']['pct_repeat_extended'] == round(100 * 2 / 48, 2)


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
