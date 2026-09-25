"""Unit tests for the Y' boundary check. Run: python _pipeline/tests/test_yprime_boundaries.py

Synthetic Y' elements on an R arm (anchor at the low coordinate), BLASTed for real:
five clean copies of one element, one copy carrying a 77-bp anchor-proximal overhang (the
chr16L_Y_Prime_1 defect), one copy truncated by 900 bp (a mis-assembly), and one with a single
different first base (BLAST alignment-end noise, which must NOT be trimmed).
"""
import os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import yprime_boundaries as yb

random.seed(7)
CORE = ''.join(random.choice('ACGT') for _ in range(3000))
JUNK = ''.join(random.choice('ACGT') for _ in range(77))
SPACER = 'N' * 500


def mutate(seq, n):
    s = list(seq)
    for i in random.sample(range(len(s)), n):
        s[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[s[i]]
    return ''.join(s)


def build(elements):
    """Lay elements out one per chromosome end on an R arm; return (regions, ref)."""
    ref, regions = {}, {}
    for i, (name, seq) in enumerate(elements):
        chrom = f'c{i}'
        ref[chrom] = SPACER + seq + SPACER
        regions[f'chr{i + 1}R'] = {'anchor': [], 'x_prime': [], 'y_prime': [
            {'chr': chrom, 'start': len(SPACER), 'end': len(SPACER) + len(seq), 'label': name}]}
    return regions, ref


def run(elements):
    regions, ref = build(elements)
    els = yb.yprime_elements(regions, ref)
    res = yb.measure(els, yb.blast_all_vs_all(els, threads=2))
    by_label = {}
    for r in res:
        by_label[regions[r['chr_end']]['y_prime'][r['index']]['label']] = r
    return regions, res, by_label


CLEAN = [(f'clean{i}', mutate(CORE, 3)) for i in range(5)]


def test_overhang_is_trimmed_and_nothing_else():
    regions, res, lab = run(CLEAN + [('over', JUNK + mutate(CORE, 3))])
    assert lab['over']['trim'] == 77, lab['over']
    assert all(lab[f'clean{i}']['trim'] == 0 for i in range(5))


def test_trim_moves_the_anchor_proximal_end():
    regions, res, lab = run(CLEAN + [('over', JUNK + mutate(CORE, 3))])
    r = lab['over']
    before = dict(regions[r['chr_end']]['y_prime'][r['index']])
    yb.apply_trims(regions, res)
    after = regions[r['chr_end']]['y_prime'][r['index']]
    # R arm: the anchor side is the LOW coordinate, so only start moves
    assert after['start'] == before['start'] + 77 and after['end'] == before['end']


def test_left_arm_trims_the_high_coordinate():
    over = JUNK + mutate(CORE, 3)
    regions, ref = build(CLEAN + [('over', over)])
    # re-home the overhanging copy on an L arm: stored reverse-complemented, anchor at high coord
    rc = over.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    ref['cL'] = SPACER + rc + SPACER
    regions['chr9L'] = {'anchor': [], 'x_prime': [], 'y_prime': [
        {'chr': 'cL', 'start': len(SPACER), 'end': len(SPACER) + len(rc), 'label': 'overL'}]}
    del regions['chr6R']
    els = yb.yprime_elements(regions, ref)
    res = yb.measure(els, yb.blast_all_vs_all(els, threads=2))
    r = [x for x in res if x['chr_end'] == 'chr9L'][0]
    assert r['trim'] == 77, r
    yb.apply_trims(regions, res)
    y = regions['chr9L']['y_prime'][0]
    assert y['end'] == len(SPACER) + len(rc) - 77 and y['start'] == len(SPACER)


def test_truncated_copy_is_flagged_not_extended():
    regions, res, lab = run(CLEAN + [('trunc', mutate(CORE, 3)[:2100])])
    assert lab['trunc']['trim'] == 0
    assert lab['trunc']['flag'].startswith('shorter_than_partners_by_'), lab['trunc']


def test_one_base_difference_is_not_trimmed():
    first = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[CORE[0]]
    regions, res, lab = run(CLEAN + [('onebase', first + mutate(CORE, 3)[1:])])
    assert all(r['trim'] == 0 for r in res)


def test_lone_element_is_untouched():
    regions, res, lab = run([('solo', CORE)])
    assert res[0]['trim'] == 0 and res[0]['flag'] == '' and res[0]['n_partners'] == 0


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
