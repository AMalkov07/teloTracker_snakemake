"""Unit tests for the recombinant Y' variant caller. Run: python _pipeline/tests/test_yprime_variants.py

Synthetic templates stand in for a strain's Y's: A is random, B differs from A at ~3 % of sites, and
C carries B's allele at half of those sites and A's at the other half (plus a few of its own) -- the
trap where a pair-only test sees an A/B SNP mix in what is really a pure copy of a third element.
"""
import os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import yprime_variants as yv

L = 5200
rng = random.Random(7)
A = ''.join(rng.choice('ACGT') for _ in range(L))


def mutate_at(seq, sites, alleles):
    s = list(seq)
    for p, b in zip(sites, alleles):
        s[p] = b
    return ''.join(s)


def other(b):
    return {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[b]


SITES = sorted(rng.sample(range(60, L - 60), 160))
B = mutate_at(A, SITES, [other(A[p]) for p in SITES])
own = sorted(rng.sample([p for p in range(60, L - 60) if p not in SITES], 30))
C = mutate_at(mutate_at(A, SITES[::2], [B[p] for p in SITES[::2]]), own, [other(A[p]) for p in own])
# B2 = B except for private sites in its first 1000 bp: identical to B over the rest
early = sorted(rng.sample([p for p in range(60, 1000) if p not in SITES], 15))
B2 = mutate_at(B, early, [other(B[p]) for p in early])


def noisy(seq, rate=0.02, seed=0):
    """Nanopore-like errors: substitutions, 1-bp insertions and deletions."""
    r = random.Random(seed)
    out = []
    for b in seq:
        x = r.random()
        if x < rate / 3:
            out.append(other(b))
        elif x < 2 * rate / 3:
            out += [b, r.choice('ACGT')]
        elif x < rate:
            continue
        else:
            out.append(b)
    return ''.join(out)


def klass(templates):
    return yv.SizeClass('Short', yv.merge_templates(templates))


CLS = klass({'chrA1': A, 'chrB1': B, 'chrC1': C, 'chrD1': B2})


def call(seq):
    return yv.call_copy(seq, [CLS])


def test_pure_copy_with_read_errors_is_pure():
    for seed in range(5):
        c = call(noisy(A, seed=seed))
        assert c['call'] == 'pure', c
        assert c['signature'] == 'chrA1', c


def test_third_element_is_not_an_AB_mix():
    for seed in range(5):
        c = call(noisy(C, seed=seed))
        assert c['call'] == 'pure' and c['signature'] == 'chrC1', c


def test_chimera_is_called_with_its_switch():
    chim = A[:2000] + B[2000:]
    c = call(noisy(chim, seed=3))
    assert c['call'] == 'recomb', c
    assert c['signature'] == 'chrA1>chrB1|chrD1', c     # B and B2 are identical after 1000 bp


def test_isolated_sites_do_not_make_a_switch():
    one_event = mutate_at(A, [3000, 3003], [other(A[3000]), other(A[3003])])
    assert call(one_event)['call'] == 'pure'
    two_apart = mutate_at(A, [1500, 4000], [B[1500] if 1500 in SITES else other(A[1500]),
                                            other(A[4000])])
    assert call(two_apart)['call'] == 'pure'


def test_consensus_recovers_the_exact_switch_interval():
    chim = A[:2000] + B[2000:]
    copies = [noisy(chim, seed=s) for s in range(30)]
    cons = yv.consensus(copies, A)
    assert cons == chim
    d = yv.describe(cons, CLS)
    assert d['n_segments'] == 2 and d['signature'] == 'chrA1>chrB1|chrD1', d
    lo, hi = map(int, d['switches'].split('-'))
    last_a = max(p for p in SITES if p < 2000) + 1          # 1-based
    first_b = min(p for p in SITES if p >= 2000) + 1
    assert (lo, hi) == (last_a, first_b), (lo, hi, last_a, first_b)
    assert d['unexplained'] == 0


def test_near_identical_templates_are_merged():
    g = yv.merge_templates({'chrA1': A, 'chrA2': mutate_at(A, [100], [other(A[100])]), 'chrB1': B})
    names = sorted(x['name'] for x in g)
    assert names == ['chrA1+1', 'chrB1'], names


def test_partial_copy_is_partial():
    assert call(A[:1500])['call'] == 'partial'


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('ok', name)
