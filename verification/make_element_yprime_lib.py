#!/usr/bin/env python3
"""Build an ELEMENT-LEVEL Y' library: one entry per reference Y' element, each with its
own ID, from a day-0 assembly + its pretelomeric_regions_*_simp.bed.

The pipeline's normal library is condensed (near-identical elements share an entry) and
then clustered, so a read copy can only ever be resolved to a group. Giving every element
its own entry makes the pipeline report, for each Y' copy on a read, the single reference
element it matches best -- the finest assignment RepeatMasker can make. Every grouping
scheme is then a relabelling of that, which is what lets different schemes be scored
against each other on exactly the same read-to-element assignment.

IDs are written as `E-<end>-<n>` (no underscore), so --y-prime-id-level family and variant
give the same answer; `origin` keeps the pipeline's own `chr4R4` location syntax so the
expected day-0 array per chromosome end still resolves.

Usage: make_element_yprime_lib.py <assembly.fasta> <simp.bed> <out.fasta> [<out_table.tsv>]
"""
import re, sys
from collections import defaultdict

asm, bed, out = sys.argv[1:4]
table = sys.argv[4] if len(sys.argv) > 4 else None

seqs, name = {}, None
for line in open(asm):
    if line.startswith('>'):
        name = line[1:].split()[0]; seqs[name] = []
    else:
        seqs[name].append(line.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
# the BED names bare chromosomes; the assembly suffixes them (chr2 -> chr2_extended)
for k in list(seqs):
    base = k.split('_')[0]
    seqs.setdefault(base, seqs[k])

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(COMP)[::-1]

rows = []
for line in open(bed):
    f = line.rstrip('\n').split('\t')
    if len(f) < 5: continue
    chrom, start, end, nm, strand = f[0], int(f[1]), int(f[2]), f[3], f[4]
    m = re.match(r'^(chr\w+?[LR])_Y_Prime_(\d+)$', nm)     # ITS_* and X/anchor lines skipped
    if not m: continue
    ce, idx = m.group(1), int(m.group(2))
    s = seqs[chrom][start:end]
    if strand == '-': s = rc(s)
    rows.append((ce, idx, s))

rows.sort(key=lambda r: (r[0], r[1]))
by_end = defaultdict(list)
for ce, idx, s in rows: by_end[ce].append(idx)

with open(out, 'w') as fh:
    for ce, idx, s in rows:
        cls = 'Long' if len(s) >= 6000 else 'Short'
        tand = 'Tandem' if len(by_end[ce]) > 1 else 'Solo'
        fh.write(f'>Y_Prime_{ce}{idx}#{cls}/{tand}/E-{ce}-{idx}\n{s}\n')

if table:
    with open(table, 'w') as fh:
        fh.write('element\tchr_end\tarray_index\tlength\n')
        for ce, idx, s in rows:
            fh.write(f'E-{ce}-{idx}\t{ce}\t{idx}\t{len(s)}\n')

print(f'{len(rows)} Y\' elements across {len(by_end)} chromosome ends -> {out}')
for ce in sorted(by_end, key=lambda c: (int(re.search(r"\d+", c).group()), c)):
    print(f'  {ce:<8} {len(by_end[ce])}')
