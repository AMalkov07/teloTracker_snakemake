"""Unit tests for telomere-end confirmation and the scaffold candidate pool.
Run: python _pipeline/tests/test_end_confirmation.py

Numbers follow 7372 day 0 chr14L: 17 complete reads with a 5-copy Y' array whose telomere
starts ~40.75 kb past the anchor, of which only one has an adapter called after the telomere.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import pandas as pd
import analyze_features as af

FIVE = ('ID2', 'ID2', 'ID1', 'ID1', 'ID1')
SIX = FIVE + ('ID1',)


def read(rid, array, start, adapter=False, repeat=120):
    return {'read_id': rid, 'array': array, 'telomere_start': start, 'repeat_length': repeat, 'adapter': adapter}


def test_telomere_start_offset_both_orientations():
    # telomere at the read start: the anchor starts 40,900 bp in, 150 of them telomere repeat
    assert af.telomere_start_offset('beginning', 45940, 40900, 45940, 150) == 40750
    # telomere at the read end: read 45,960 bp, anchor ends at 5,040
    assert af.telomere_start_offset('end', 45960, 0, 5040, 170) == 40750
    assert af.telomere_start_offset('end', 45960, -1, -1, 170) is None


def test_reads_ending_together_confirm_each_other():
    reads = [read(f'r{i}', FIVE, 40750 + d) for i, d in enumerate((-20, -5, 0, 8, 30))]
    reads.append(read('adapter', FIVE, 40760, adapter=True))
    ev = af.concordant_end_confirmation(reads)
    assert ev['adapter'] == 'adapter'
    assert all(ev[f'r{i}'] == 'concordant_ends' for i in range(5))


def test_too_few_or_scattered_reads_stay_unconfirmed():
    ev = af.concordant_end_confirmation([read('a', FIVE, 40750), read('b', FIVE, 40760)])
    assert ev == {'a': '', 'b': ''}                        # two reads are not enough
    ev = af.concordant_end_confirmation([read('a', FIVE, 30000), read('b', FIVE, 35000), read('c', FIVE, 40750)])
    assert set(ev.values()) == {''}                        # three reads, three different ends


def test_different_arrays_do_not_pool():
    reads = [read('a', FIVE, 40750), read('b', FIVE, 40755), read('c', SIX, 40752)]
    assert set(af.concordant_end_confirmation(reads).values()) == {''}


def test_short_repeat_does_not_count():
    reads = [read(f'r{i}', FIVE, 40750, repeat=20) for i in range(4)]
    assert set(af.concordant_end_confirmation(reads).values()) == {''}


def test_spacer_break_guard():
    # 3 reads stop after 2 copies where 40 reads continue into a third: the signature of reads
    # broken inside the same ITS. Without an adapter-confirmed member they stay unconfirmed;
    # with one, the cluster is a real (terminal-loss) end.
    stop = [read(f's{i}', ('ID2', 'ID2'), 11500 + i) for i in range(3)]
    cont = [read(f'c{i}', FIVE, 40750 + (i % 7)) for i in range(40)]
    ev = af.concordant_end_confirmation(stop + cont)
    assert all(ev[f's{i}'] == '' for i in range(3))
    assert all(ev[f'c{i}'] == 'concordant_ends' for i in range(40))
    stop[0]['adapter'] = True
    ev = af.concordant_end_confirmation(stop + cont)
    assert ev['s0'] == 'adapter' and ev['s1'] == ev['s2'] == 'concordant_ends'


def test_apply_end_confirmation_updates_telo_info():
    info = {f'r{i}': {'confirmed': False, 'evidence': '', 'repeat_length': 150.0, 'probe_count': 5}
            for i in range(3)}
    reads = [(f'r{i}', 'beginning', 45940 + i * 10, 40900 + i * 10, 45940 + i * 10, 'ID2,ID2,ID1,ID1,ID1')
             for i in range(3)]
    assert af.apply_end_confirmation(info, reads) == 3
    assert all(v['confirmed'] and v['evidence'] == 'concordant_ends' for v in info.values())
    assert af.apply_end_confirmation(None, reads) == 0


def pool_df(n_adapter, n_plain, probe=5):
    rows = [{'read_id': f'a{i}', 'chr_end': '14L', 'repeat_length': 150 + i, 'Adapter_After_Telomere': True,
             'y_prime_probe_count': probe} for i in range(n_adapter)]
    rows += [{'read_id': f'p{i}', 'chr_end': '14L', 'repeat_length': 60 + i, 'Adapter_After_Telomere': False,
              'y_prime_probe_count': probe} for i in range(n_plain)]
    rows += [{'read_id': 'short', 'chr_end': '14L', 'repeat_length': 10, 'Adapter_After_Telomere': False,
              'y_prime_probe_count': probe}]
    return pd.DataFrame(rows)


def test_scaffold_pool_prefers_adapter_reads():
    su = __import__('subtelomere_reference_pipeline_utils')
    pool, rule = su.scaffold_candidate_pool(pool_df(4, 15), '14L', min_agree=3)
    assert rule == 'adapter' and len(pool) == 4


def test_scaffold_pool_widens_when_adapter_reads_cannot_vote():
    su = __import__('subtelomere_reference_pipeline_utils')
    pool, rule = su.scaffold_candidate_pool(pool_df(2, 16), '14L', min_agree=3)
    assert rule == 'telomere_repeat' and len(pool) == 18      # the 10 bp-repeat read is excluded
    assert list(pool['repeat_length']) == sorted(pool['repeat_length'])
    pool, rule = su.scaffold_candidate_pool(pool_df(0, 0), '14L', min_agree=3)
    assert len(pool) == 0


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
