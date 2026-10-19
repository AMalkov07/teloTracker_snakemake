"""Unit tests for the profiling options: the qualifying-read rule, the stand-in for the adapter check,
splitting contigs into per-end halves, and how run_pipeline.py writes and checks these settings.
Run: python _pipeline/tests/test_profile_options.py
"""
import gzip, os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(HERE, '..', 'scripts'), os.path.join(HERE, '..', '..')]
import pandas as pd
import telomere_filters as tf
import no_adapter_check as nac
import split_contig_ends as sce
import run_pipeline as rp


def exits(fn, *args):
    try:
        fn(*args)
    except SystemExit:
        return True
    return False


def test_qualifying_needs_the_adapter_only_when_it_was_checked():
    df = pd.DataFrame({'repeat_length': [250, 250, 29, None, 250],
                       'Adapter_After_Telomere': [True, False, True, True, float('nan')]})
    assert list(tf.qualifying_mask(df, require_adapter=True)) == [True, False, False, False, False]
    assert list(tf.qualifying_mask(df, require_adapter=False)) == [True, True, False, False, True]
    text = pd.DataFrame({'repeat_length': ['250', '250'], 'Adapter_After_Telomere': ['True', 'False']})
    assert list(tf.qualifying_mask(text)) == [True, False]


def test_stand_in_adapter_table_keeps_every_anchored_read_unchecked():
    top = pd.DataFrame({'read_id': ['r1', 'r2'], 'repeat_type': ['AC', 'TG'],
                        'anchor_name': ['chr1L_anchor', 'chr9R_anchor'], 'pident': [99.0, 98.0]})
    t = nac.adapter_table(top)
    assert list(t.columns) == nac.COLUMNS and list(t['Repeat_Type']) == ['AC', 'TG']
    assert (t['Adapter_After_Telomere'] == '').all()
    # blank = not checked: never counted as an adapter, and the read still qualifies when the
    # adapter is not required
    t['repeat_length'] = 100
    assert not tf.qualifying_mask(t, require_adapter=True).any()
    assert tf.qualifying_mask(t, require_adapter=False).all()


def test_contigs_are_cut_into_one_half_per_end():
    assert sce.halves('chr1', 'AACCGGTT') == [('chr1__start', 'AACC'), ('chr1__end', 'GGTT')]
    assert sce.halves('c', 'ACG') == [('c__start', 'A'), ('c__end', 'CG')]
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, 'asm.fa.gz')
        with gzip.open(src, 'wt') as fh:
            fh.write('>chrI some description\nAAAA\nCCCC\n>chrII\nGGTT\n')
        out = os.path.join(d, 'out.fasta')
        sys.argv = ['split_contig_ends.py', src, out]
        sce.main()
        assert open(out).read() == '>chrI__start\nAAAA\n>chrI__end\nCCCC\n>chrII__start\nGG\n>chrII__end\nTT\n'


def test_driver_writes_profiling_options():
    cfg = rp.normalize_cfg({'base_name': 's', 'adapter_check': False, 'min_qscore': 0})
    text = rp.snakemake_config_text(cfg, recombination=False)
    assert 'input_type: "reads"' in text and 'adapter_check: false' in text
    assert 'min_qscore: 0' in text and 'min_read_length: 2000' in text
    text = rp.snakemake_config_text(rp.normalize_cfg({'base_name': 's', 'input_type': 'contigs'}), False)
    assert 'input_type: "contigs"' in text and 'adapter_check: false' in text   # never for contigs
    text = rp.snakemake_config_text(rp.normalize_cfg({'base_name': 's'}), False)
    assert 'adapter_check: true' in text


def test_driver_rejects_bad_profiling_options():
    ok = rp.normalize_cfg({'base_name': 's', 'input_type': 'contigs'})
    assert not exits(rp.validate_cfg, ok, ['profile'])
    assert exits(rp.validate_cfg, dict(ok, strain='x'), ['recombination'])   # contigs: profile only
    assert exits(rp.validate_cfg, rp.normalize_cfg({'base_name': 's', 'input_type': 'scaffolds'}), ['profile'])
    assert exits(rp.validate_cfg, rp.normalize_cfg({'base_name': 's', 'adapter_check': 'no'}), ['profile'])
    assert exits(rp.validate_cfg, rp.normalize_cfg({'base_name': 's', 'min_qscore': -1}), ['profile'])


def test_sample_inputs_follow_the_input_type():
    with tempfile.TemporaryDirectory() as root:
        os.makedirs(os.path.join(root, 'samples_dorado_basecalled'))
        open(os.path.join(root, 'samples_dorado_basecalled', 'asm.fna.gz'), 'w').close()
        contigs = rp.normalize_cfg({'base_name': 'asm', 'input_type': 'contigs'})
        assert rp.missing_sample_inputs(contigs, ['asm'], root) == []
        reads = rp.normalize_cfg({'base_name': 'asm'})
        assert len(rp.missing_sample_inputs(reads, ['asm'], root)) == 1


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
