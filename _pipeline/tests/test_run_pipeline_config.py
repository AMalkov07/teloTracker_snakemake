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


def test_profile_config_needs_no_strain_or_day0_keys():
    text = rp.snakemake_config_text(normalized(base_name='s'), recombination=False)
    assert 'recombination: false' in text and 'base_name: "s"' in text
    for key in ('strain:', 'day0_ref', 'day0_bed', 'y_prime_lib', 'spacer_lib_dir', 'x_element_lib',
                'attribution_mode'):
        assert key not in text, key
    assert 'anchors: "_pipeline/references/telomerase_shutoff_anchors.fasta"' in text
    assert not exits(rp.validate_cfg, normalized(base_name='s'), ['profile'])
    assert exits(rp.validate_cfg, normalized(), ['profile'])           # no sample at all


def test_recombination_config_points_at_the_day0_sample():
    cfg = normalized(base_name='t3', strain='7172_curated', _day0_base='d0')
    text = rp.snakemake_config_text(cfg, recombination=True)
    assert 'recombination: true' in text and 'strain: "7172_curated"' in text
    assert 'day0_bed: "results/d0/_pipeline/pretelomeric_labels/pretelomeric_regions_7172_curated_simp.bed"' in text
    assert 'x_element_lib: "results/d0/_pipeline/pretelomeric_labels/clustered_x_elements_{strain}.fasta"' in text


def test_create_ref_needs_an_existing_dorado_image_in_docker_mode():
    base = dict(base_name='s', strain='x')
    assert exits(rp.validate_cfg, normalized(**base), ['create_ref'])                   # no image
    assert exits(rp.validate_cfg, normalized(**base, dorado_image='/no/such.sif'), ['create_ref'])
    assert not exits(rp.validate_cfg, normalized(**base, dorado_image=__file__), ['create_ref'])
    assert not exits(rp.validate_cfg, normalized(**base, dorado_mode='local'), ['create_ref'])


def test_label_libraries_are_checked_and_patched_in():
    base = dict(base_name='s', strain='x')
    assert exits(rp.validate_cfg, normalized(**base, label_x_fasta='no/such.fasta'), ['label_regions'])
    y = '_pipeline/references/repeatmasker_6991_all_y_primes.chr16Lfix.fasta'
    assert not exits(rp.validate_cfg, normalized(**base, label_yprime_fasta=y), ['label_regions'])
    script = open(os.path.join(rp.PIPELINE_DIR, '_pipeline', 'label_regions.sh')).read()
    out = rp.substitute_var(script, 'YPRIMES_FASTA', y, 'label_regions.sh')
    assert f'YPRIMES_FASTA="{y}"' in out


def test_recombination_prerequisites_are_listed_before_anything_runs():
    import tempfile
    with tempfile.TemporaryDirectory() as root:
        cfg = normalized(base_name='d0', strain='S', timepoint_samples=['t3'])
        problems = rp.check_recombination_prereqs(cfg, rp.step_samples(cfg), root)
        assert len(problems) == 6, problems            # reference, BED, 3 libraries, t3's input
        assert any('spacer library' in p for p in problems) and any('t3' in p for p in problems)
        labels = os.path.join(root, 'results', 'd0', '_pipeline', 'pretelomeric_labels')
        os.makedirs(os.path.join(labels, 'pairings_for_spacers', 'S_pairings'))
        os.makedirs(os.path.join(root, 'results', 'd0', '_pipeline', 'assembly_S'))
        os.makedirs(os.path.join(root, 'samples_dorado_basecalled'))
        for f in ('assembly_S/assembly_S_dorado_reference.fasta',
                  'pretelomeric_labels/pretelomeric_regions_S_simp.bed',
                  'pretelomeric_labels/extracted_yprimes_S.fasta',
                  'pretelomeric_labels/clustered_x_elements_S.fasta',
                  'pretelomeric_labels/pairings_for_spacers/S_pairings/chr1L.fasta'):
            open(os.path.join(root, 'results', 'd0', '_pipeline', f), 'w').close()
        open(os.path.join(root, 'samples_dorado_basecalled', 't3.fq.gz'), 'w').close()
        assert rp.check_recombination_prereqs(cfg, ['t3'], root) == []


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
