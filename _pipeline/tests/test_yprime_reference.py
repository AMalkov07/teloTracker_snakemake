"""Unit tests for the expected Y' counts (yprime_reference.py) and how run_pipeline.py picks them.
Run: python _pipeline/tests/test_yprime_reference.py
"""
import os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(HERE, '..', 'scripts'), os.path.join(HERE, '..', '..')]
import yprime_reference as yr
import run_pipeline as rp

REFS = os.path.join(HERE, '..', 'references')


def curated_bed(name):
    return os.path.join(REFS, 'day0', name, f'pretelomeric_regions_{name}_simp.bed')


def test_curated_beds_give_the_old_table_counts():
    # the values the hard-coded tables in y_prime_analysis.py held for these strains
    c = yr.load_expected_counts(curated_bed('6991_curated'))
    assert len(c) == 32 and sum(c.values()) == 34
    assert c['chr4R'] == 7 and c['chr12R'] == 7 and c['chr14L'] == 5 and c['chr1L'] == 0
    c = yr.load_expected_counts(curated_bed('7172_curated'))
    assert sum(c.values()) == 35 and c['chr12R'] == 8 and c['chr14L'] == 3 and c['chr16L'] == 3
    c = yr.load_expected_counts(curated_bed('7302_curated'))
    assert sum(c.values()) == 36 and c['chr12R'] == 6 and c['chr13L'] == 4


def test_bed_counts_skip_its_and_unlabelled_ends():
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as fh:
        fh.write('chr4\t1\t10\tchr4R_anchor\t+\t10\n'
                 'chr4\t20\t30\tchr4R_Y_Prime_1\t+\t11\n'
                 'chr4\t40\t50\tITS_chr4R_Y_Prime_1-2\t+\t11\n'
                 'chr4\t60\t70\tchr4R_Y_Prime_2\t+\t11\n'
                 'chr1\t1\t10\tchr1L_anchor\t-\t10\n')
    c = yr.load_expected_counts(fh.name)
    os.unlink(fh.name)
    assert c == {'chr1L': 0, 'chr4R': 2}          # chr1L labelled with no Y'; other ends unknown


def test_shipped_counts_tables_cover_all_ends():
    for name in ('7021', '7372', 'KRLT3_7637', 'tlc1'):
        c = yr.load_expected_counts(os.path.join(REFS, 'yprime_expected_counts', f'{name}.tsv'))
        assert len(c) == 32
    assert sum(yr.load_expected_counts(os.path.join(REFS, 'yprime_expected_counts', 'tlc1.tsv')).values()) == 9


def test_incomplete_counts_table_is_an_error():
    with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as fh:
        fh.write('chr_end\texpected_y_primes\nchr1L\t0\n')
    try:
        yr.load_expected_counts(fh.name)
        raise AssertionError('accepted a table without all 32 ends')
    except ValueError:
        pass
    finally:
        os.unlink(fh.name)


def test_resolve_order():
    path, label = yr.resolve('auto', '7172_curated', None, REFS)
    assert path.endswith('pretelomeric_regions_7172_curated_simp.bed') and 'matched by strain' in label
    assert yr.resolve('auto', '7172', None, REFS)[0].endswith('7172_curated_simp.bed')   # <name>_curated
    assert yr.resolve('auto', '7372', None, REFS)[0].endswith('yprime_expected_counts/7372.tsv')
    # a strain with nothing shipped: the day-0 BED in recombination, else no reference
    bed = curated_bed('6991_curated')
    assert yr.resolve('auto', '6991_day0_with_selection', bed, REFS) == (bed, 'day-0 label BED')
    assert yr.resolve('auto', '6991_day0_with_selection', None, REFS)[0] is None
    assert yr.resolve('none', '7172_curated', bed, REFS)[0] is None
    assert yr.resolve('tlc1', '6991_curated', None, REFS)[0].endswith('tlc1.tsv')
    for bad in ('KRLT36', 'no_such.bed'):
        try:
            yr.resolve(bad, None, None, REFS)
            raise AssertionError(f'accepted {bad}')
        except ValueError:
            pass


def test_driver_writes_the_reference_into_the_config():
    cfg = rp.normalize_cfg({'base_name': 's', 'strain': '7302_curated'})
    text = rp.snakemake_config_text(cfg, recombination=False)
    assert ('yprime_reference_path: "_pipeline/references/day0/7302_curated/'
            'pretelomeric_regions_7302_curated_simp.bed"') in text
    text = rp.snakemake_config_text(rp.normalize_cfg({'base_name': 's'}), recombination=False)
    assert 'yprime_reference_path: ""' in text
    text = rp.snakemake_config_text(rp.normalize_cfg({'base_name': 's', 'yprime_reference': 'tlc1'}),
                                    recombination=False)
    assert 'yprime_expected_counts/tlc1.tsv' in text
    bad = rp.normalize_cfg({'base_name': 's', 'yprime_reference': 'KRLT36'})
    try:
        rp.validate_cfg(bad, ['profile'])
        raise AssertionError('accepted an unknown yprime_reference')
    except SystemExit:
        pass


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
