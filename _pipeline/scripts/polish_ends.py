#!/usr/bin/env python3
"""Ends-only polishing for create_ref.sh (POLISH_SCOPE="ends").

The extended reference is <L extension> + <base contig> + <R extension>: the base reference
(e.g. 6991_only_to_anchors.fasta) with each arm extended to the telomere from one scaffold read.
Genome-wide polishing with every read is fine when the library covers the whole genome, but an
end-enriched library leaves chromosome interiors at a few x, and polishing at that depth corrupts
them (KRLT43: Flye removed 1.33 of 12.25 Mb). In ends mode each arm is polished on its own:

  prepare  cut one segment per arm (extension + FLANK bp of the base contig) and write the reads
           anchored to that arm (chr_end column of the probe TSV) beside it
  (create_ref.sh then polishes every segment with only its own reads: Flye, then dorado)
  splice   put the polished segments back onto the unpolished base contigs

The junction is placed at a k-mer taken from the middle of the flank that occurs exactly once in
the polished segment, so the join is always between sequence the two pieces agree on.

Usage:
  polish_ends.py prepare --extended EXT.fa --base BASE.fa --tsv PROBE.tsv --reads TRIMMED.fastq
                         --flank 5000 --min-reads 5 --outdir DIR
  polish_ends.py splice  --extended EXT.fa --outdir DIR --polished-name NAME --output OUT.fa
"""
import argparse
import csv
import os
import re
import sys

K = 40              # junction k-mer length
PROBE_STEP = 50     # spacing between junction probes tried along the flank


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs


def write_fasta(path, records):
    with open(path, 'w') as f:
        for name, seq in records:
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + '\n')


def extension_lengths(ext_seq, base_seq):
    """(left, right) extension lengths of an extended contig around its base contig."""
    p = ext_seq.find(base_seq[:2000])
    if p < 0 or ext_seq[p:p + len(base_seq)] != base_seq:
        raise ValueError('base contig not found intact inside the extended contig')
    return p, len(ext_seq) - p - len(base_seq)


def segments(ext, base, flank):
    """Per arm: (end, contig, side, ext_len, flank, start, stop) in extended coordinates.
    The flank is shrunk if the two arms' segments would meet (short base contig)."""
    out = []
    for contig in sorted(ext):
        b = contig[:-len('_extended')] if contig.endswith('_extended') else contig
        if b not in base:
            raise ValueError(f'{contig}: no base contig {b} in the base reference')
        L = len(ext[contig])
        lx, rx = extension_lengths(ext[contig], base[b])
        f = min(flank, (L - lx - rx) // 2)
        chrom = b[3:] if b.startswith('chr') else b
        out.append((f'{chrom}L', contig, 'L', lx, f, 0, lx + f))
        out.append((f'{chrom}R', contig, 'R', rx, f, L - rx - f, L))
    return out


def find_junction(orig_seg, polished, side, flank):
    """Splice point: (index in orig_seg, index in polished) of a K-mer from the flank that
    occurs once in each. Probes start mid-flank and walk inward (towards the interior, where
    the flank sequence is copied unchanged) then outward. None if no probe is unique."""
    n = len(orig_seg)
    if side == 'L':
        lo, hi = n - flank, n - K          # flank = last `flank` bp of the left segment
    else:
        lo, hi = 0, flank - K              # flank = first `flank` bp of the right segment
    if hi < lo:
        return None
    mid = (lo + hi) // 2
    inward = range(mid, hi + 1, PROBE_STEP) if side == 'L' else range(mid, lo - 1, -PROBE_STEP)
    outward = range(mid, lo - 1, -PROBE_STEP) if side == 'L' else range(mid, hi + 1, PROBE_STEP)
    for i in list(inward) + list(outward):
        kmer = orig_seg[i:i + K]
        if 'N' in kmer or orig_seg.count(kmer) != 1:
            continue
        j = polished.find(kmer)
        if j >= 0 and polished.find(kmer, j + 1) < 0:
            return i, j
    return None


def cmd_prepare(a):
    ext, base = read_fasta(a.extended), read_fasta(a.base)
    segs = segments(ext, base, a.flank)
    arm_reads = {}
    with open(a.tsv) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r.get('chr_end'):
                arm_reads.setdefault(r['read_id'], set()).add(r['chr_end'])
    # a read anchored to two arms polishes neither (multi-end chimeras are already removed upstream)
    read_arm = {rid: next(iter(e)) for rid, e in arm_reads.items() if len(e) == 1}
    os.makedirs(a.outdir, exist_ok=True)
    handles, counts = {}, {}
    for end, contig, side, lx, f, s, e in segs:
        d = os.path.join(a.outdir, end)
        os.makedirs(d, exist_ok=True)
        write_fasta(os.path.join(d, 'segment.fasta'), [(f'{contig}__{side}', ext[contig][s:e])])
        handles[end] = open(os.path.join(d, 'reads.fastq'), 'w')
        counts[end] = 0
    with open(a.reads) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s, p, q = f.readline(), f.readline(), f.readline()
            end = read_arm.get(h[1:].split()[0])
            if end in handles:
                handles[end].write(h + s + '+\n' + q)
                counts[end] += 1
    for h in handles.values():
        h.close()
    with open(os.path.join(a.outdir, 'manifest.tsv'), 'w') as m:
        m.write('end\tcontig\tside\textension_bp\tflank_bp\tsegment_bp\tanchored_reads\tpolish\n')
        for end, contig, side, lx, f, s, e in segs:
            ok = counts[end] >= a.min_reads
            m.write(f'{end}\t{contig}\t{side}\t{lx}\t{f}\t{e - s}\t{counts[end]}\t{"yes" if ok else "no"}\n')
            print(f'  {end:<4} {contig} {side}: extension {lx:,} bp + flank {f:,} bp, '
                  f'{counts[end]} anchored reads' + ('' if ok else f'  -> NOT polished (< {a.min_reads} reads)'))


def cmd_splice(a):
    ext = read_fasta(a.extended)
    with open(os.path.join(a.outdir, 'manifest.tsv')) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    by_contig = {}
    for r in rows:
        by_contig.setdefault(r['contig'], {})[r['side']] = r
    out, report = [], []
    for contig in sorted(ext, key=lambda c: [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', c)]):
        seq = ext[contig]
        L = len(seq)
        left, right = seq[:0], seq[L:]
        cut_l, cut_r = 0, L
        for side in ('L', 'R'):
            r = by_contig[contig][side]
            seg_len, flank = int(r['segment_bp']), int(r['flank_bp'])
            orig = seq[:seg_len] if side == 'L' else seq[L - seg_len:]
            status, pol = 'unpolished', None
            ppath = os.path.join(a.outdir, r['end'], a.polished_name)
            if r['polish'] == 'yes' and os.path.exists(ppath):
                p = read_fasta(ppath)
                pol = next(iter(p.values()), None) if p else None
                status = 'polished' if pol else 'polisher returned nothing'
            elif r['polish'] == 'yes':
                status = f'missing {a.polished_name}'
            junc = find_junction(orig, pol, side, flank) if pol else None
            if pol and junc is None:
                status = 'no unique junction k-mer; left unpolished'
            if junc:
                i, j = junc
                if side == 'L':
                    left, cut_l = pol[:j], i
                    new_len = j + (seg_len - i)
                else:
                    right, cut_r = pol[j:], L - seg_len + i
                    new_len = i + (len(pol) - j)
            else:
                new_len = seg_len
                if side == 'L':
                    cut_l = 0
                else:
                    cut_r = L
            report.append((r['end'], contig, side, r['anchored_reads'], seg_len, new_len, status))
        out.append((contig, left + seq[cut_l:cut_r] + right))
    write_fasta(a.output, out)
    with open(os.path.join(a.outdir, f'splice_report_{os.path.splitext(a.polished_name)[0]}.tsv'), 'w') as f:
        f.write('end\tcontig\tside\tanchored_reads\tsegment_bp\tpolished_segment_bp\tstatus\n')
        for row in report:
            f.write('\t'.join(map(str, row)) + '\n')
    n_pol = sum(r[6] == 'polished' for r in report)
    print(f'  spliced {n_pol}/{len(report)} polished arms into {a.output}')
    for row in report:
        if row[6] != 'polished':
            print(f'  {row[0]:<4} {row[6]}')
        elif abs(row[5] - row[4]) > 0.05 * row[4]:
            print(f'  {row[0]:<4} WARNING polished arm is {row[5]:,} bp vs {row[4]:,} bp before polishing')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--extended', required=True)
    p.add_argument('--base', required=True)
    p.add_argument('--tsv', required=True, help='post_y_prime_probe.tsv (read_id, chr_end)')
    p.add_argument('--reads', required=True, help='adapter-trimmed FASTQ')
    p.add_argument('--flank', type=int, default=5000)
    p.add_argument('--min-reads', type=int, default=5)
    p.add_argument('--outdir', required=True)
    s = sub.add_parser('splice')
    s.add_argument('--extended', required=True)
    s.add_argument('--outdir', required=True)
    s.add_argument('--polished-name', required=True, help='file name of each arm\'s polished segment')
    s.add_argument('--output', required=True)
    a = ap.parse_args()
    {'prepare': cmd_prepare, 'splice': cmd_splice}[a.cmd](a)


if __name__ == '__main__':
    sys.exit(main())
