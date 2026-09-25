"""Unit tests for the v3 confidence scores. Run: python _pipeline/tests/test_confidence.py

recombination_confidence: is the read's call right?  donor_confidence: is the named donor right?
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import analyze_features as af


def rc(**kw):
    args = dict(recombinant=True, spacer_switch=False, spacer_conf=0, x_switch=False, x_conf=0,
                y_status='No Change', n_gained=0, downstream_consistent=True, telo_confirmed=True)
    args.update(kw)
    return af.recombination_confidence(**args)[0]


def test_unchanged_read_depends_on_reaching_the_telomere():
    assert rc(recombinant=False, telo_confirmed=True) == af.NO_CHANGE_CONFIRMED
    assert rc(recombinant=False, telo_confirmed=False) == af.NO_CHANGE_UNCONFIRMED
    assert rc(recombinant=False, telo_confirmed=None) == af.NO_CHANGE_CONFIRMED


def test_more_gained_copies_more_confidence():
    one, two, four = (rc(y_status="Y' Gain", n_gained=n) for n in (1, 2, 4))
    assert 0 < one < two < four < 1


def test_loss_on_unconfirmed_end_scores_lower():
    assert rc(y_status="Y' Loss", telo_confirmed=False) < rc(y_status="Y' Loss", telo_confirmed=True)


def test_id_change_weaker_than_multi_copy_gain_and_mixed_array_weaker_still():
    clean = rc(y_status="1st Y' Change", downstream_consistent=True)
    mixed = rc(y_status="1st Y' Change", downstream_consistent=False)
    assert mixed < clean < rc(y_status="Y' Gain", n_gained=3)


def test_independent_evidence_accumulates():
    y_only = rc(y_status="Y' Gain", n_gained=1)
    both = rc(y_status="Y' Gain", n_gained=1, spacer_switch=True, spacer_conf=0.45)
    assert both > y_only and both <= 1.0


def test_no_floor():
    # the old score never went below 0.30 for a recombinant read; the new one reflects weak evidence
    assert rc(spacer_switch=True, spacer_conf=0.05) < 0.1


def dc(best, pool, evidence):
    return af.donor_confidence(best, pool, evidence)[0]


def test_ambiguous_donor_is_zero():
    assert dc('ambiguous', {'chr13L': 1.0, 'chr14L': 1.0}, [('y_prime', 'chr13L', 0.5)]) == 0.0
    assert dc('', {}, []) == 0.0


def test_unique_long_fingerprint_beats_short_one():
    long_ = dc('chr13L', {'chr13L': 1.0}, [('y_prime', 'chr13L', 1.0)])
    short = dc('chr13L', {'chr13L': 0.6}, [('y_prime', 'chr13L', 0.67)])
    assert long_ > short > 0


def test_close_vote_lowers_confidence():
    clear = dc('chr11R', {'chr11R': 1.5}, [('spacer', 'chr11R', 0.8)])
    close = dc('chr11R', {'chr11R': 1.5, 'chr11L': 1.4}, [('spacer', 'chr11R', 0.8)])
    assert close < clear


def test_dissenting_evidence_lowers_confidence():
    # a weak spacer switch names chr7R while the Y' array fingerprints chr13L: the old code ignored
    # the Y' evidence entirely once any spacer / X switch fired
    agree = dc('chr7R', {'chr7R': 1.0}, [('spacer', 'chr7R', 0.6)])
    dissent = dc('chr7R', {'chr7R': 1.0}, [('spacer', 'chr7R', 0.6), ('y_prime', 'chr13L', 1.0)])
    assert dissent < agree


def test_agreeing_axes_raise_confidence():
    one = dc('chr11R', {'chr11R': 1.0}, [('x_element', 'chr11R', 0.5)])
    two = dc('chr11R', {'chr11R': 1.5}, [('x_element', 'chr11R', 0.5), ('spacer', 'chr11R', 0.5)])
    assert two > one


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
