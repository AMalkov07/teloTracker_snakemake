#!/usr/bin/env python3
"""Trim an over-extended Y' element boundary back to where its near-identical partners start.

Background: the curated 6991 library's `Y_Prime_chr16L1` is 6737 bp while `chr7R1` is 6656 and
`chr14L1` is 6657, and the three are 99.94-99.97 % identical over their shared 6655 bp. The excess
sits at chr16L1's anchor-proximal end and is not Y' sequence (0 % telomeric repeat; it aligns into
chr7R's annotated x_variable_element). The labelling step reproduces it in every assembly, and the
coverage-penalised clustering similarity then turns the offset into a spurious group split.

This measures the offset instead of assuming it: the element is BLASTed against the partners it is
>= 99 % identical to, and the trim is the query position where those alignments consistently begin.
Writes a CORRECTED COPY of the BED -- the input is never modified.

Usage: fix_yprime_boundary.py <assembly.fasta> <simp.bed> <out.bed> [--element chr16L_Y_Prime_1]
                              [--min-donors 2] [--max-trim 300] [--report <tsv>]
"""
import os, re, subprocess, sys, tempfile
from collections import defaultdict

asm, bed, out_bed = sys.argv[1:4]
def opt(f, d):
    return sys.argv[sys.argv.index(f) + 1] if f in sys.argv else d
TARGET    = opt('--element', 'chr16L_Y_Prime_1')
MIN_DON   = int(opt('--min-donors', 2))
MAX_TRIM  = int(opt('--max-trim', 300))
REPORT    = opt('--report', None)

seqs, name = {}, None
for line in open(asm):
    if line.startswith('>'): name = line[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(line.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
for k in list(seqs): seqs.setdefault(k.split('_')[0], seqs[k])

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
rc = lambda s: s.translate(COMP)[::-1]

# --- read the BED; Y' rows are 1-based inclusive as written by write_bed_simplified -------
rows = [l.rstrip('\n').split('\t') for l in open(bed) if l.strip()]
yp = {}
for i, f in enumerate(rows):
    if len(f) < 5: continue
    m = re.match(r'^(chr\w+?[LR])_Y_Prime_(\d+)$', f[3])
    if m: yp[f[3]] = (i, f[0], int(f[1]), int(f[2]), f[4], m.group(1), int(m.group(2)))

if TARGET not in yp:
    print(f'{os.path.basename(bed)}: {TARGET} absent -- nothing to do'); sys.exit(0)

def elem_seq(nm):
    _, chrom, s, e, strand, _, _ = yp[nm]
    sub = seqs[chrom][s - 1:e]          # 1-based inclusive, matching the BED
    return rc(sub) if strand == '-' else sub   # element position 1 = anchor-proximal

# --- BLAST the target against every other element ----------------------------------------
with tempfile.TemporaryDirectory(prefix='ypfix_') as td:
    db = os.path.join(td, 'db.fasta')
    with open(db, 'w') as fh:
        for nm in yp:
            if nm != TARGET: fh.write(f'>{nm}\n{elem_seq(nm)}\n')
    q = os.path.join(td, 'q.fasta')
    open(q, 'w').write(f'>{TARGET}\n{elem_seq(TARGET)}\n')
    subprocess.run(['makeblastdb', '-in', db, '-dbtype', 'nucl', '-out', os.path.join(td, 'db')],
                   check=True, capture_output=True)
    r = subprocess.run(['blastn', '-query', q, '-db', os.path.join(td, 'db'),
                        '-outfmt', '6 sseqid pident length qlen slen qstart qend sstart send',
                        '-evalue', '1e-10'], check=True, capture_output=True, text=True)

best = {}
for line in r.stdout.splitlines():
    f = line.split('\t')
    s, pid, ln = f[0], float(f[1]), int(f[2])
    if s not in best or ln > best[s][1]: best[s] = (pid, ln, int(f[3]), int(f[4]), int(f[5]), int(f[6]))

qlen = len(elem_seq(TARGET))
donors = []
for s, (pid, ln, ql, sl, qs, qe) in best.items():
    if pid < 99.0: continue
    if ln < 0.90 * max(ql, sl): continue          # keeps Short Y' out of a Long partner set
    donors.append((s, qs - 1, pid, sl))           # qs-1 = bp of query before the shared start

trim = 0
if len(donors) >= MIN_DON:
    offs = sorted(d[1] for d in donors)
    cand = offs[len(offs) // 2]                   # median implied offset
    agree = [o for o in offs if abs(o - cand) <= 5]
    if len(agree) >= MIN_DON and 0 < cand <= MAX_TRIM:
        trim = cand

print(f'{os.path.basename(bed)}')
print(f'  {TARGET}: {qlen} bp, {len(donors)} partner(s) >=99 % identical')
for s, off, pid, sl in sorted(donors, key=lambda d: -d[2])[:6]:
    print(f'    {s:<22} {pid:6.3f}%  {sl:>5} bp  implies trimming {off} bp')
if not trim:
    print('  -> no confident offset; BED copied unchanged')
else:
    print(f'  -> trimming {trim} bp from the anchor-proximal end ({qlen} -> {qlen - trim} bp)')

# --- rewrite: shrink the Y' at its anchor-proximal end, give the bases to ITS_0-1 ---------
if trim:
    idx, chrom, s, e, strand, ce, num = yp[TARGET]
    if strand == '-':   new_s, new_e = s, e - trim        # L arm: anchor side is the high coord
    else:               new_s, new_e = s + trim, e
    rows[idx][1], rows[idx][2] = str(new_s), str(new_e)
    if len(rows[idx]) > 5: rows[idx][5] = str(new_e - new_s + 1)
    its = f'ITS_{ce}_Y_Prime_{num-1}-{num}' if num > 1 else f'ITS_{ce}_Y_Prime_0-1'
    for j, f in enumerate(rows):
        if len(f) > 3 and f[3] == its:
            if strand == '-': rows[j][1] = str(new_e + 1)
            else:             rows[j][2] = str(new_s - 1)
            if len(rows[j]) > 5: rows[j][5] = str(int(rows[j][2]) - int(rows[j][1]) + 1)
            print(f'  -> {its} extended to {rows[j][1]}-{rows[j][2]} ({rows[j][5]} bp)')
            break
    else:
        print(f'  -> note: no {its} row to extend (the bases fall to the X element instead)')

with open(out_bed, 'w') as fh:
    for f in rows: fh.write('\t'.join(f) + '\n')
if REPORT:
    with open(REPORT, 'a') as fh:
        fh.write(f'{os.path.basename(bed)}\t{TARGET}\t{qlen}\t{trim}\t{qlen-trim}\t{len(donors)}\n')
print(f'  written {out_bed}')
