"""
Unit tests for the v2 source attribution in analyze_features.py.

Run with   python -m pytest _pipeline/tests/test_attribution.py -q
or plain   python _pipeline/tests/test_attribution.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'scripts'))

import analyze_features as af  # noqa: E402

# The 7302 day-0 Y' library as built on Argon (extracted_yprimes_7302_day0_with_selection.fasta)
LIB_7302 = [
    'Y_Prime_chr12R2,3,4,5;chr4R1,2,3,4,6,7#Long/Tandem/ID1_Gray',
    'Y_Prime_chr14L3,4,5#Short/Tandem/ID2_Red',
    'Y_Prime_chr13L1,3#Short/Tandem/ID2_Red',
    'Y_Prime_chr13L2,4#Long/Tandem/ID1_Gray',
    'Y_Prime_chr10L1#Long/Solo/ID3_Green',
    'Y_Prime_chr12L1#Short/Solo/ID2_Red',
    'Y_Prime_chr12R1#Long/Solo/ID5_Purple',
    'Y_Prime_chr12R6#Long/Solo/ID1_Gray',
    'Y_Prime_chr14L1#Long/Solo/ID1_Gray',
    'Y_Prime_chr14L2#Long/Solo/ID1_Gray',
    'Y_Prime_chr14R1#Long/Solo/ID6_Blue',
    'Y_Prime_chr15R1#Long/Solo/ID1_Gray',
    'Y_Prime_chr16L1#Long/Solo/ID1_Gray',
    'Y_Prime_chr16R1#Short/Solo/ID7_Yellow',
    'Y_Prime_chr2L1#Short/Solo/ID4_Orange',
    'Y_Prime_chr4R5#Long/Solo/ID1_Gray',
    'Y_Prime_chr5L1#Long/Solo/ID8_Cyan',
    'Y_Prime_chr5R1#Long/Solo/ID1_Gray',
    'Y_Prime_chr6L1#Short/Solo/ID4_Orange',
    'Y_Prime_chr7R1#Long/Solo/ID1_Gray',
    'Y_Prime_chr8L1#Short/Solo/ID2_Red',
    'Y_Prime_chr8R1#Short/Solo/ID2_Red',
    'Y_Prime_chr9L1#Long/Solo/ID3_Green',
]


def _lib():
    return {h.split()[0]: af.parse_y_prime_header(h) for h in LIB_7302}


def setup_module(module=None):
    af.ATTRIBUTION_MODE = 'v2'


# ---------------------------------------------------------------------------
# C0 helpers
# ---------------------------------------------------------------------------

def test_parse_origin_locations_multi_end():
    locs = af.parse_origin_locations('Y_Prime_chr12R2,3,4,5;chr4R1,2,3,4,6,7')
    assert ('chr4R', 7) in locs and ('chr12R', 2) in locs and len(locs) == 10


def test_reference_arrays():
    arrays = af.build_reference_arrays(_lib())
    assert arrays['chr13L'] == ['ID2', 'ID1', 'ID2', 'ID1']
    assert arrays['chr4R'] == ['ID1'] * 7
    assert arrays['chr14L'] == ['ID1', 'ID1', 'ID2', 'ID2', 'ID2']
    assert arrays['chr12R'] == ['ID5', 'ID1', 'ID1', 'ID1', 'ID1', 'ID1']


def test_norm_end_and_same_end():
    assert af._norm_end('chr10L') == (10, 'L')
    assert af._norm_end('chr4_extended') == (4, None)
    assert af._same_end('chr4', 'chr4R')
    assert not af._same_end('chr1', 'chr10L')          # the old substring bug
    assert not af._same_end('chr4L', 'chr4R')


# ---------------------------------------------------------------------------
# C1 fingerprint matching
# ---------------------------------------------------------------------------

def test_alternating_fingerprint_is_unique_to_chr13L():
    arrays = af.build_reference_arrays(_lib())
    m = af.match_gained_array(['ID2', 'ID1', 'ID2', 'ID1'], arrays)
    assert m['chr13L'] == 'contiguous' and all(k != 'contiguous' for ce, k in m.items() if ce != 'chr13L')
    for gained in (['ID2', 'ID1', 'ID2', 'ID1'], ['ID1', 'ID2', 'ID1'], ['ID2', 'ID1'] * 3):
        fp = af.find_fingerprint(gained, arrays, 'chr2L')
        assert fp['best'] == ['chr13L'] and fp['source'] == 'chr13L', (gained, fp)


def test_two_element_gain_is_not_specific():
    arrays = af.build_reference_arrays(_lib())
    m = af.match_gained_array(['ID1', 'ID1'], arrays)
    assert {'chr4R', 'chr12R', 'chr14L'} <= set(m)
    assert set(af.match_gained_array(['ID1'], arrays)) >= {'chr4R', 'chr12R', 'chr14L', 'chr5R'}   # single common ID: many ends
    assert af.find_fingerprint(['ID1'], arrays, 'chr2L')['source'] == ''                            # -> not a fingerprint
    assert af.find_fingerprint(['ID8'], arrays, 'chr2L')['source'] == 'chr5L'                       # unique variant -> fingerprint


def test_single_unique_yprime_names_donor():
    lib = _lib()
    arrays = af.build_reference_arrays(lib)
    y = af.compare_y_prime_arrays(['ID4', 'ID8'], _ref(['ID4']), lib, arrays, 'chr2L')   # ID8 exists only at chr5L
    rec = af.reconcile_features(_no_spacer('chr2L'), _no_x('chr2L'), y, [], 'chr2L')
    assert rec['recombination_source'] == 'chr5L' and rec['recombination_mechanism'] == 'donor_transfer'
    y = af.compare_y_prime_arrays(['ID4', 'ID1'], _ref(['ID4']), lib, arrays, 'chr2L')   # ID1 is everywhere
    rec = af.reconcile_features(_no_spacer('chr2L'), _no_x('chr2L'), y, [], 'chr2L')
    assert rec['recombination_source'] == 'ambiguous'


def test_rotation_and_periodic():
    arrays = {'donor': ['A', 'B', 'C']}
    assert af.match_gained_array(['B', 'C', 'A'], arrays) == {'donor': 'rotation'}
    assert af.match_gained_array(['A', 'B', 'A', 'B'], {'d': ['A', 'B']}) == {'d': 'periodic'}
    assert af.match_gained_array(['A', 'B', 'A', 'B', 'A'], {'d': ['A', 'B']}) == {'d': 'periodic'}


def test_homopolymer_run():
    arrays = af.build_reference_arrays(_lib())
    m = af.match_gained_array(['ID1'] * 5, arrays)
    assert m['chr4R'] == 'contiguous' and m['chr12R'] == 'contiguous' and m['chr14L'] == 'periodic'
    assert af.match_gained_array(['ID1'] * 9, arrays)['chr4R'] == 'periodic'   # longer than any run


def test_find_fingerprint_prefers_contiguous_and_trims_ragged_tail():
    arrays = af.build_reference_arrays(_lib())
    # ID1,ID2,ID1 is contiguous in chr13L and only a rotation in chr14L -> chr13L wins
    fp = af.find_fingerprint(['ID1', 'ID2', 'ID1'], arrays, 'chr2L')
    assert fp['best'] == ['chr13L'] and fp['source'] == 'chr13L' and fp['best_kind'] == 'contiguous'
    # ragged read: long alternating stretch followed by an irregular tail
    fp = af.find_fingerprint(['ID2', 'ID1', 'ID2', 'ID1', 'ID2', 'ID1', 'ID2', 'ID1', 'ID1', 'ID2'], arrays, 'chr2L')
    assert fp['source'] == 'chr13L' and fp['matched_len'] == 8 and fp['best_kind'] == 'periodic', fp
    # self end explains the gain -> no donor
    fp = af.find_fingerprint(['ID2', 'ID1', 'ID2'], arrays, 'chr13L')
    assert fp['self_match'] and fp['source'] == ''
    # a long homopolymer run after the alternation must not out-rank it
    fp = af.find_fingerprint(['ID2', 'ID1', 'ID2', 'ID1'] + ['ID2'] * 7, arrays, 'chr11L')
    assert fp['source'] == 'chr13L' and fp['window'][:4] == ['ID2', 'ID1', 'ID2', 'ID1'] and fp['specificity'] == 1.0, fp
    # a pure run still matches (low specificity)
    fp = af.find_fingerprint(['ID1'] * 4, arrays, 'chr2L')
    assert fp['matched_len'] == 4 and set(fp['best']) == {'chr4R', 'chr12R'}


# ---------------------------------------------------------------------------
# C3 compatible ends
# ---------------------------------------------------------------------------

def test_compatible_ends_sees_every_location():
    af.ATTRIBUTION_MODE = 'v2'
    ends = af.find_compatible_ends(['ID1'], _lib())
    assert 'chr4R' in ends and 'chr14L' in ends and 'chr12R' not in ends
    # legacy read only the first chr_end / first position of a multi-location
    # header, so an end named second in a shared record was invisible
    mini = {h.split()[0]: af.parse_y_prime_header(h) for h in
            ['Y_Prime_chr12R2,3;chr4R1,2#Long/Tandem/ID1_Gray', 'Y_Prime_chr12R1#Long/Solo/ID5_Purple']}
    af.ATTRIBUTION_MODE = 'legacy'
    legacy = af.find_compatible_ends(['ID1'], mini)
    af.ATTRIBUTION_MODE = 'v2'
    assert 'chr4R' not in legacy
    assert af.find_compatible_ends(['ID1'], mini) == ['chr4R']


# ---------------------------------------------------------------------------
# compare_y_prime_arrays + reconcile (C2, C5, C7)
# ---------------------------------------------------------------------------

def _ref(ids):
    return [{'feature_name': f'x_Y_Prime_{i+1}', 'id': v, 'start': 0, 'end': 0} for i, v in enumerate(ids)]


def _no_spacer(chr_end):
    return {'spacer_source': chr_end, 'spacer_recombination': 'no_change', 'spacer_confidence': 0.95}


def _no_x(chr_end):
    return {'x_element_source': chr_end, 'x_element_recombination': 'no_change', 'x_element_confidence': 0.95}


def test_yprime_only_gain_of_fingerprint_attributed_to_chr13L():
    lib = _lib()
    arrays = af.build_reference_arrays(lib)
    # chr2L (ref = ID4) gains the chr13L array
    y = af.compare_y_prime_arrays(['ID4', 'ID2', 'ID1', 'ID2', 'ID1'], _ref(['ID4']), lib, arrays, 'chr2L')
    assert y['y_prime_recombination_status'] == "Y' Gain"
    assert y['y_prime_gained_segment'] == 'ID2,ID1,ID2,ID1'
    assert y['y_prime_fingerprint_source'] == 'chr13L'
    rec = af.reconcile_features(_no_spacer('chr2L'), _no_x('chr2L'), y, ['chr5_extended'], 'chr2L')
    assert rec['recombination_source'] == 'chr13L'
    assert rec['recombination_mechanism'] == 'donor_transfer'
    assert rec['overall_confidence'] >= 0.5
    assert rec['source_resolution'] == 'arm'


def test_armless_supplementary_cannot_beat_arm_resolved():
    lib = _lib()
    arrays = af.build_reference_arrays(lib)
    y = af.compare_y_prime_arrays(['ID4', 'ID2', 'ID1', 'ID2'], _ref(['ID4']), lib, arrays, 'chr2L')
    sp = {'spacer_source': 'chr8L', 'spacer_recombination': 'switch_detected', 'spacer_confidence': 0.8}
    rec = af.reconcile_features(sp, _no_x('chr2L'), y, ['chr5'], 'chr2L')
    # the spacer names the proximal donor; the Y' fingerprint (chr13L) becomes the
    # secondary y_prime_donor and the event is complex; chr5 (arm-less) never wins
    assert rec['recombination_source'] == 'chr8L'
    assert rec['y_prime_donor'] == 'chr13L'
    assert rec['is_complex_event'] is True
    assert rec['source_tie'] is False


def test_structural_donor_keeps_source_when_yprime_points_elsewhere():
    # 7172 chr11L -> chr11R BIR: x-element says chr11R, Y' array fingerprints chr14L
    lib = _lib()
    arrays = af.build_reference_arrays(lib)
    y = af.compare_y_prime_arrays(['ID1', 'ID1', 'ID2', 'ID2'], _ref([]), lib, arrays, 'chr11L')
    x = {'x_element_source': 'chr11R', 'x_element_recombination': 'full_switch', 'x_element_confidence': 0.9}
    rec = af.reconcile_features(_no_spacer('chr11L'), x, y, ['chr11'], 'chr11L')
    assert rec['recombination_source'] == 'chr11R'
    assert rec['y_prime_donor'] == 'chr14L'
    assert rec['recombination_mechanism'] == 'subtelomere_switch'


def test_armless_only_evidence_keeps_chromosome():
    lib = _lib()
    y = af.compare_y_prime_arrays(['ID4', 'ID1'], _ref(['ID4']), lib, af.build_reference_arrays(lib), 'chr2L')   # ID1: 6 ends, no Y' vote
    rec = af.reconcile_features(_no_spacer('chr2L'), _no_x('chr2L'), y, ['chr9'], 'chr2L')
    assert rec['recombination_source'] == 'chr9'
    assert rec['source_resolution'] == 'chromosome'
    assert rec['is_complex_event'] is False


def test_same_chromosome_armless_and_arm_agree_not_complex():
    lib = _lib()
    y = af.compare_y_prime_arrays(['ID4'], _ref(['ID4']), lib, af.build_reference_arrays(lib), 'chr2L')
    y['y_prime_recombination_status'] = 'No Change'
    sp = {'spacer_source': 'chr4R', 'spacer_recombination': 'switch_detected', 'spacer_confidence': 0.8}
    rec = af.reconcile_features(sp, _no_x('chr2L'), y, ['chr4'], 'chr2L')
    assert rec['recombination_source'] == 'chr4R'
    assert rec['is_complex_event'] is False         # legacy flagged 'chr4' vs 'chr4R' as complex
    assert rec['recombination_mechanism'] == 'subtelomere_switch'


def test_tandem_amplification_same_end():
    lib = _lib()
    arrays = af.build_reference_arrays(lib)
    y = af.compare_y_prime_arrays(['ID2', 'ID1', 'ID2', 'ID1', 'ID2', 'ID1'], _ref(['ID2', 'ID1', 'ID2', 'ID1']), lib, arrays, 'chr13L')
    assert y['y_prime_recombination_status'] == "Y' Gain"
    rec = af.reconcile_features(_no_spacer('chr13L'), _no_x('chr13L'), y, [], 'chr13L')
    assert rec['recombination_mechanism'] == 'tandem_amplification_same_end'
    assert rec['recombination_source'] == 'chr13L'         # own array explains it -> source = self
    assert rec['source_resolution'] == 'self'
    # chr12R gaining more ID1 copies is also self-explained, even though chr4R could supply ID1 runs
    y = af.compare_y_prime_arrays(['ID5'] + ['ID1'] * 8, _ref(['ID5'] + ['ID1'] * 5), lib, arrays, 'chr12R')
    rec = af.reconcile_features(_no_spacer('chr12R'), _no_x('chr12R'), y, ['chr4'], 'chr12R')
    assert rec['recombination_mechanism'] == 'tandem_amplification_same_end'
    assert rec['recombination_source'] == 'chr12R'         # a lone supplementary hit is not a donor call
    assert rec['source_resolution'] == 'self'


def test_loss_confirmation():
    lib = _lib()
    y = af.compare_y_prime_arrays(['ID2', 'ID1'], _ref(['ID2', 'ID1', 'ID2', 'ID1']), lib, af.build_reference_arrays(lib), 'chr13L')
    assert y['y_prime_recombination_status'] == "Y' Loss"
    ok = af.reconcile_features(_no_spacer('chr13L'), _no_x('chr13L'), y, [], 'chr13L', {'confirmed': True})
    bad = af.reconcile_features(_no_spacer('chr13L'), _no_x('chr13L'), y, [], 'chr13L', {'confirmed': False})
    assert ok['recombination_mechanism'] == 'array_contraction'
    assert bad['recombination_mechanism'] == 'array_contraction_unconfirmed'
    assert bad['overall_confidence'] < ok['overall_confidence']
    assert 'loss_unconfirmed_end' in bad['qc_flags_extra']


def test_spacer_switch_confidence_reflects_per_chunk_gap():
    # 10 chunks: 4 match chr2L at 99%, then 6 match chr8L at 99% while chr2L only reaches 92%
    hits = []
    for i in range(4):
        hits.append((i * 250, {'source': 'chr2L', 'pident': 99.0, 'bitscore': 400, 'expected_pident': 99.0}))
    for i in range(4, 10):
        hits.append((i * 250, {'source': 'chr8L', 'pident': 99.0, 'bitscore': 400, 'expected_pident': 92.0}))
    af.ATTRIBUTION_MODE = 'legacy'
    leg = af.analyze_chunks('r', hits, 'chr2L', 'spacer')
    af.ATTRIBUTION_MODE = 'v2'
    new = af.analyze_chunks('r', hits, 'chr2L', 'spacer')
    assert leg['spacer_recombination'] == new['spacer_recombination'] == 'switch_detected'
    assert leg['spacer_source'] == 'chr8L' and new['spacer_source'] == 'chr8L'
    assert leg['spacer_confidence'] < 0.05                # the broken separation-based score
    assert abs(new['spacer_confidence'] - 0.9 * 0.7) < 1e-6   # 7-point gap -> 0.63
    assert new['spacer_plurality_source'] == 'chr8L'


def test_low_confidence_spacer_switch_does_not_override_fingerprint():
    lib = _lib()
    arrays = af.build_reference_arrays(lib)
    y = af.compare_y_prime_arrays(['ID4', 'ID2', 'ID1', 'ID2', 'ID1'], _ref(['ID4']), lib, arrays, 'chr2L')
    sp = {'spacer_source': 'chr12R', 'spacer_recombination': 'switch_detected', 'spacer_confidence': 0.002}
    rec = af.reconcile_features(sp, _no_x('chr2L'), y, [], 'chr2L')
    assert rec['recombination_source'] == 'chr13L'
    assert rec['y_prime_donor'] == 'chr13L'


def test_spacer_walk_restricted_to_spacer_interval():
    # anchor 0-5000, spacer 5000-10000, x 10000-10700, Y' 11000-40000 (telo_side=end)
    start, end = af.spacer_interval('end', 40000, 0, 5000, 10000, 10700, 11000, 40000)
    assert (start, end) == (5000, 10000)
    hits = [(p, {'source': 'chr12R' if p >= 11000 else 'chr2L', 'pident': 99.0, 'bitscore': 400,
                 'expected_pident': 0.0 if p >= 11000 else 99.0}) for p in range(0, 40000, 250)]
    kept = af.chunks_in_interval(hits, start, end)
    assert all(5000 <= p < 10000 for p, _ in kept) and len(kept) == 20
    af.ATTRIBUTION_MODE = 'v2'
    assert af.analyze_chunks('r', kept, 'chr2L', 'spacer')['spacer_recombination'] == 'no_change'
    assert af.analyze_chunks('r', hits, 'chr2L', 'spacer')['spacer_recombination'] != 'no_change'   # the legacy artefact
    # mirrored read (telomere at the beginning)
    assert af.spacer_interval('beginning', 40000, 35000, 40000, 29300, 30000, 0, 29000) == (30000, 35000)


def test_legacy_mode_unchanged_shape():
    af.ATTRIBUTION_MODE = 'legacy'
    try:
        lib = _lib()
        y = af.compare_y_prime_arrays(['ID4', 'ID2', 'ID1', 'ID2', 'ID1'], _ref(['ID4']), lib, af.build_reference_arrays(lib), 'chr2L')
        assert y['y_prime_fingerprint_source'] == ''      # fingerprint disabled in legacy
        rec = af.reconcile_features(_no_spacer('chr2L'), _no_x('chr2L'), y, ['chr5_extended'], 'chr2L')
        assert rec['recombination_source'] == 'chr5'         # the legacy behaviour the truth set exposes
        assert rec['recombination_mechanism'] == ''
    finally:
        af.ATTRIBUTION_MODE = 'v2'


if __name__ == '__main__':
    setup_module()
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn()
                print(f'PASS {name}')
            except AssertionError as e:
                failed += 1
                print(f'FAIL {name}: {e}')
            except Exception as e:  # noqa: BLE001
                failed += 1
                print(f'ERROR {name}: {type(e).__name__}: {e}')
    sys.exit(1 if failed else 0)
