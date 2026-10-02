"""Unit tests for ends-only polishing (polish_ends.py).
Run: python _pipeline/tests/test_polish_ends.py
"""
import os, random, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'scripts'))
import polish_ends as pe

random.seed(1)
rnd = lambda n: ''.join(random.choice('ACGT') for _ in range(n))


def test_extension_lengths():
    base = rnd(10000)
    assert pe.extension_lengths(rnd(300) + base + rnd(700), base) == (300, 700)
    assert pe.extension_lengths(base, base) == (0, 0)


def test_segments_shrink_flank_on_short_contig():
    base = rnd(6000)
    ext = {'chr1_extended': rnd(100) + base + rnd(200)}
    segs = pe.segments(ext, {'chr1': base}, 5000)
    assert [s[:3] for s in segs] == [('1L', 'chr1_extended', 'L'), ('1R', 'chr1_extended', 'R')]
    (_, _, _, lx, fl, s0, e0), (_, _, _, rx, fr, s1, e1) = segs
    assert (lx, rx, fl, fr) == (100, 200, 3000, 3000)
    assert e0 <= s1                                   # the two arms never overlap


def test_junction_survives_edits_away_from_flank_middle():
    seg = rnd(3000)                                   # extension 2000 + flank 1000
    polished = seg[:500] + 'GG' + seg[500:1200] + seg[1203:]   # insertion + deletion in the extension
    i, j = pe.find_junction(seg, polished, 'L', 1000)
    assert seg[i:i + pe.K] == polished[j:j + pe.K]
    assert polished[:j] + seg[i:] == polished         # splice reproduces the polished arm
    rseg = rnd(3000)                                  # flank 1000 + extension 2000
    rpol = rseg[:2200] + 'T' + rseg[2200:]
    i, j = pe.find_junction(rseg, rpol, 'R', 1000)
    assert rseg[:i] + rpol[j:] == rpol


def test_prepare_and_splice_end_to_end():
    base = rnd(20000)
    lx, rx = rnd(4000), rnd(6000)
    with tempfile.TemporaryDirectory() as d:
        ext_fa, base_fa = os.path.join(d, 'ext.fa'), os.path.join(d, 'base.fa')
        pe.write_fasta(ext_fa, [('chr1_extended', lx + base + rx)])
        pe.write_fasta(base_fa, [('chr1', base)])
        tsv = os.path.join(d, 'probe.tsv')
        with open(tsv, 'w') as f:
            f.write('read_id\tchr_end\n' + ''.join(f'L{i}\t1L\n' for i in range(6)) +
                    ''.join(f'R{i}\t1R\n' for i in range(3)) + 'X\t1L\nX\t1R\n')
        fq = os.path.join(d, 'reads.fq')
        with open(fq, 'w') as f:
            for rid in [f'L{i}' for i in range(6)] + [f'R{i}' for i in range(3)] + ['X', 'other']:
                f.write(f'@{rid} x\nACGT\n+\nIIII\n')
        out = os.path.join(d, 'ends')
        run = lambda *a: subprocess.run([sys.executable, os.path.join(HERE, '..', 'scripts', 'polish_ends.py'), *a],
                                        check=True, capture_output=True, text=True)
        run('prepare', '--extended', ext_fa, '--base', base_fa, '--tsv', tsv, '--reads', fq,
            '--flank', '2000', '--min-reads', '5', '--outdir', out)
        man = open(os.path.join(out, 'manifest.tsv')).read().splitlines()
        assert man[1].split('\t')[6:] == ['6', 'yes'] and man[2].split('\t')[6:] == ['3', 'no']
        assert open(os.path.join(out, '1L', 'reads.fastq')).read().count('@') == 6   # X (two arms) excluded
        seg = pe.read_fasta(os.path.join(out, '1L', 'segment.fasta'))['chr1_extended__L']
        assert seg == lx + base[:2000]
        pe.write_fasta(os.path.join(out, '1L', 'pol.fasta'), [('x', seg[:100] + 'AAAA' + seg[100:])])
        run('splice', '--extended', ext_fa, '--outdir', out, '--polished-name', 'pol.fasta',
            '--output', os.path.join(d, 'out.fa'))
        assert pe.read_fasta(os.path.join(d, 'out.fa'))['chr1_extended'] == lx[:100] + 'AAAA' + lx[100:] + base + rx


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn(); print('ok', name)
