"""Cut each assembled contig into its two halves so profiling can treat them like reads.

A whole chromosome carries the anchors of both of its ends, and a sequence hitting two different
ends' anchors is dropped as a chimera (filter_for_reads_with_anchors.py). Each half holds one end:
<contig>__start is the first half as written, <contig>__end the second half, both in the contig's
own orientation, so the telomere stays at the sequence end as it would in a read.

Usage: split_contig_ends.py <contigs.fasta[.gz]> <out.fasta>
"""
import gzip
import sys


def read_fasta(path):
    opener = gzip.open if path.endswith('.gz') else open
    name, seq = None, []
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, ''.join(seq)


def halves(name, seq):
    mid = len(seq) // 2
    return [(f'{name}__start', seq[:mid]), (f'{name}__end', seq[mid:])]


def main():
    src, out = sys.argv[1:3]
    n = 0
    with open(out, 'w') as fh:
        for name, seq in read_fasta(src):
            for piece, s in halves(name, seq):
                fh.write(f'>{piece}\n{s}\n')
            n += 1
    print(f'split_contig_ends.py: {n} contigs -> {2 * n} halves in {out}')


if __name__ == '__main__':
    main()
