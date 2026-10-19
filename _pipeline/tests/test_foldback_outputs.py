"""Unit tests for fold-back detection outputs (detect_foldback_reads.py).
Run: python _pipeline/tests/test_foldback_outputs.py

Snakemake declares both the ID list and the detail table as outputs, so both must exist even when
a sample has no fold-backs or no anchored reads at all (a missing one failed the whole run).
"""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import pandas as pd
import detect_foldback_reads as dfr

BLAST_HEADER = ('read_id\ttotal_read_length\tread_bp_used_for_match\tmatch_start_on_read\tmatch_end_on_read\t'
                'anchor_name\ttotal_anchor_length\tmatch_start_on_anchor\tmatch_end_on_anchor\tpident\tbitscore\tevalue\n')


def test_no_foldbacks_still_writes_both_files():
    with tempfile.TemporaryDirectory() as d:
        ids, detail = os.path.join(d, 'ids.txt'), os.path.join(d, 'detail.tsv')
        dfr.write_outputs({}, [], ids, detail)
        assert open(ids).read() == ''
        df = pd.read_csv(detail, sep='\t')
        assert df.empty and list(df.columns) == dfr.DETAIL_COLUMNS


def test_foldbacks_are_written_sorted():
    with tempfile.TemporaryDirectory() as d:
        ids, detail = os.path.join(d, 'ids.txt'), os.path.join(d, 'detail.tsv')
        row = dict.fromkeys(dfr.DETAIL_COLUMNS, 1)
        dfr.write_outputs({'r2': 'chr1L_anchor', 'r1': 'chr1L_anchor'}, [row], ids, detail)
        assert open(ids).read().split() == ['r1', 'r2']
        assert len(pd.read_csv(detail, sep='\t')) == 1


def test_no_anchored_reads_writes_both_files():
    with tempfile.TemporaryDirectory() as d:
        blast, reads = os.path.join(d, 'blast.tsv'), os.path.join(d, 'reads.fasta')
        ids, detail = os.path.join(d, 'ids.txt'), os.path.join(d, 'detail.tsv')
        open(blast, 'w').write(BLAST_HEADER)          # header only: no anchored reads
        open(reads, 'w').write('>r1\nACGT\n')
        open(reads + '.fai', 'w').write('r1\t4\t4\t4\t5\n')
        argv = sys.argv
        sys.argv = ['detect_foldback_reads.py', '--blast-tsv', blast, '--reads-fasta', reads,
                    '--anchors-fasta', reads, '--output', ids, '--detail-tsv', detail]
        try:
            dfr.main()
        finally:
            sys.argv = argv
        assert os.path.exists(ids) and os.path.exists(detail)


if __name__ == '__main__':
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn(); print('PASS', name)
            except Exception as e:  # noqa
                failed += 1; print('FAIL', name, type(e).__name__, e)
    sys.exit(1 if failed else 0)
