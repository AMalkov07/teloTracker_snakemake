"""Unit tests for Y' clustering with very few elements. Run: python _pipeline/tests/test_yprime_clustering.py

A strain can carry 0, 1 or 2 Y' variants (6212's reference has exactly one). Those cases used
to crash the labelling run (no tree to build) or, at 2 variants, split them unconditionally.
"""
import os, random, subprocess, sys, tempfile
import numpy as np
SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts')
sys.path.insert(0, SCRIPTS)
import cluster_yprimes_paper_method as cy


def D(sim):
    """Distance matrix from a symmetric similarity matrix (percent)."""
    return 100.0 - np.asarray(sim, dtype=float)


def test_zero_and_one_element_do_not_crash():
    for fn in (cy.find_clusters_silhouette, cy.find_clusters_threshold):
        k, Z, scores, labels = fn(np.zeros((0, 0)))
        assert (k, Z, len(labels)) == (1, None, 0)
        k, Z, scores, labels = fn(np.zeros((1, 1)))
        assert (k, Z, list(labels)) == (1, None, [1])


def test_two_near_identical_variants_are_one_group():
    # previously silhouette forced k=2 whenever no k could be scored
    k, _, _, labels = cy.find_clusters_silhouette(D([[100, 99.5], [99.5, 100]]), fallback_threshold=97.0)
    assert k == 1 and len(set(labels)) == 1


def test_two_divergent_variants_are_two_groups():
    k, _, _, labels = cy.find_clusters_silhouette(D([[100, 85.0], [85.0, 100]]), fallback_threshold=97.0)
    assert k == 2 and len(set(labels)) == 2


def test_fallback_uses_the_configured_threshold():
    sim = D([[100, 98.0], [98.0, 100]])
    assert cy.find_clusters_silhouette(sim, fallback_threshold=97.0)[0] == 1
    assert cy.find_clusters_silhouette(sim, fallback_threshold=99.0)[0] == 2


def test_three_or_more_still_uses_silhouette():
    sim = D([[100, 99.8, 80, 80], [99.8, 100, 80, 80], [80, 80, 100, 99.7], [80, 80, 99.7, 100]])
    k, Z, scores, labels = cy.find_clusters_silhouette(sim)
    assert Z is not None and scores and k == 2


def _run_main(n_seqs):
    random.seed(n_seqs)
    with tempfile.TemporaryDirectory() as td:
        fa, out = os.path.join(td, 'in.fasta'), os.path.join(td, 'out.fasta')
        with open(fa, 'w') as fh:
            for i in range(n_seqs):
                seq = ''.join(random.choice('ACGT') for _ in range(600))
                fh.write(f'>Y_Prime_chr{i + 1}L1#Long/Solo/ID1_Gray\n{seq}\n')
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, 'cluster_yprimes_paper_method.py'), fa,
                            '--output-fasta', out, '--output-dir', td],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stdout[-800:] + r.stderr[-800:]
        return [l for l in open(out) if l.startswith('>')]


def test_main_with_no_yprimes_writes_empty_fasta():
    assert _run_main(0) == []


def test_main_with_one_yprime_writes_one_group():
    headers = _run_main(1)
    assert len(headers) == 1 and 'ID1' in headers[0]


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
