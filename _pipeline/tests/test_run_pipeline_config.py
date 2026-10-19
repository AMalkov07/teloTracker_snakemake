"""Unit tests for run_pipeline.py's config handling: normalisation, the Snakemake config it writes,
validation of the Y' grouping, and shell-variable patching.
Run: python _pipeline/tests/test_run_pipeline_config.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import run_pipeline as rp


def normalized(**cfg):
    return rp.normalize_cfg(dict(cfg))


def exits(fn, *args):
    try:
        fn(*args)
    except SystemExit:
        return True
    return False


def test_numeric_strain_and_bam_suffix_are_normalised():
    cfg = normalized(base_name='sample_day0.bam', strain=7302, timepoint_samples=['d3.bam', 4])
    assert cfg['strain'] == '7302' and cfg['base_name'] == 'sample_day0'
    assert cfg['timepoint_samples'] == ['d3', '4']


def test_snakemake_config_uses_normalised_name_and_writes_attribution_keys():
    text = rp.snakemake_config_text(normalized(base_name='s.bam', strain=7302))
    assert 'base_name: "s"' in text and 'strain: "7302"' in text
    assert 'attribution_mode: "v2"' in text and 'y_prime_id_level: "family"' in text
    assert 'yprime_identity_threshold' not in text          # read by label_regions only
    text = rp.snakemake_config_text(normalized(base_name='s', strain='x', y_prime_id_level='group'))
    assert 'y_prime_id_level: "group"' in text


def test_yprime_grouping_resolves_to_the_older_keys():
    cfg = normalized(yprime_grouping=99.5)
    assert cfg['yprime_stop_mode'] == 'threshold' and cfg['yprime_identity_threshold'] == 99.5
    cfg = normalized(yprime_grouping='silhouette', yprime_identity_threshold=97.0)
    assert cfg['yprime_stop_mode'] == 'silhouette' and cfg['yprime_identity_threshold'] == 97.0
    assert exits(lambda: normalized(yprime_grouping='strict'))


def test_threshold_given_as_a_fraction_is_rejected():
    cfg = normalized(base_name='s', strain='x', yprime_grouping=0.995)
    assert exits(rp.validate_cfg, cfg, ['label_regions'])
    cfg = normalized(base_name='s', strain='x', yprime_grouping=99.5)
    assert not exits(rp.validate_cfg, cfg, ['label_regions'])


def test_recombination_requires_strain():
    assert exits(rp.validate_cfg, normalized(base_name='s'), ['recombination'])


def test_substitute_var_forms():
    script = ('A="old"  # keep me\n'
              "  export B='old value'\n"
              'C=old\r\n'
              '#A="commented out"\n')
    out = rp.substitute_var(script, 'A', 'new', 't.sh')
    out = rp.substitute_var(out, 'B', 'x y', 't.sh')
    out = rp.substitute_var(out, 'C', r'C:\dir "q"', 't.sh')
    assert 'A="new"  # keep me\n' in out
    assert '  export B="x y"\n' in out
    assert 'C="C:\\\\dir \\"q\\""\r\n' in out
    assert '#A="commented out"' in out


def test_substitute_var_missing_variable_is_an_error():
    assert exits(rp.substitute_var, 'X=1\n', 'Y', '2', 't.sh')


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
