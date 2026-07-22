#!/usr/bin/env python3
"""
detect_foldback_reads.py -- flag fold-back (hairpin) reads.

A fold-back read is one molecule whose complementary strand was sequenced
immediately after the template, so the anchor + subtelomere appear TWICE in
inverted orientation. A real molecule cannot contain its own anchor twice.

Such a read is dangerous: create_ref picks one read per chr_arm as a scaffold
and grafts its soft-clipped overhang onto the reference. A fold-back carries a
second, spurious overhang that can win extend_reference_multi's max(len)
comparison and be grafted onto the WRONG chromosome arm. On 7172 this replaced
chr16R's subtelomere with an inverted copy of chr16L's, making every chr16R
read disagree with its own reference (100% false recombination).

WHY A SEPARATE, PERMISSIVE BLAST IS NEEDED
------------------------------------------
The main anchor BLAST runs with -min_raw_gapped_score, which filters on the
PRELIMINARY gapped score (before traceback), not the score it reports. The
second copy of a fold-back is always later in the read, and ONT accuracy decays
~2 percentage points per 100 kb, so the gappier second copy is preferentially
discarded before traceback ever runs. The 7172 chr16L scaffold is the worked
example: its reverse copy has a final raw score of 9010 -- comfortably above the
5000 cutoff -- yet never appears in the strict table. Measured on 7172/7858/7871
the cliff sits between 5000 and 4500; at 3000 the copy is recovered with no
false positives (1,199/1,200 control reads unaffected).

SCOPE: DETECTION ONLY
---------------------
Each read is compared ONLY against the anchor it was already assigned. This
pass therefore cannot discover a new anchor, cannot reassign a read, and cannot
make an unanchored read anchored -- the anchored set can only shrink by the
exclusions this produces. It is also ~7.6x less work than re-BLASTing every
read (12,576 anchored vs 95,816 total on 7172) and ~32x fewer anchor
comparisons, which matters because these subtelomeric reads are highly
repetitive (mid_occ ~4,400).

CRITERIA (validated on 32 scaffolds x 3 strains)
-----------------------------------------------
  same anchor, >=2 hits, each >= --min-len, in OPPOSITE orientations.

No positional constraint: the hairpin point can fall anywhere. On 7172 the
culprit's copies sat at 21% and 78% of the read, while other fold-backs had
copies only ~5 kb apart, so "near opposite ends" would miss most of them.

Length, not identity, is the discriminator: real second copies ran ~5,070 bp
while the largest repeat-noise fragment was 565 bp (~9x separation). Identity
is a poor proxy -- detected second copies range from 91.2% to 99.4%.
"""
import argparse
import collections
import os
import subprocess
import sys
import tempfile

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description='Detect fold-back (hairpin) reads')
    p.add_argument('--blast-tsv', required=True,
                   help='raw anchor BLAST table (used only for read->anchor assignment)')
    p.add_argument('--reads-fasta', required=True,
                   help='the same FASTA the anchor BLAST used (needs a .fai)')
    p.add_argument('--anchors-fasta', required=True)
    p.add_argument('--output', required=True, help='text file of fold-back read IDs')
    p.add_argument('--detail-tsv', default=None, help='optional per-hit detail for review')
    p.add_argument('--min-len', type=int, default=2000,
                   help='each copy must align at least this many bp (default 2000)')
    p.add_argument('--perc-identity', type=float, default=85.0)
    p.add_argument('--threads', type=int, default=4)
    return p.parse_args()


def load_fai(fai_path):
    idx = {}
    with open(fai_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 5:
                idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return idx


def fetch(fh, rec):
    """Pull one sequence using its .fai record."""
    length, offset, line_bases, line_width = rec
    fh.seek(offset)
    n_lines = (length + line_bases - 1) // line_bases
    raw = fh.read(n_lines * line_width).decode()
    return ''.join(raw.split())[:length]


def write_fasta(path, records):
    with open(path, 'w') as out:
        for name, seq in records:
            out.write(f'>{name}\n{seq}\n')


def main():
    a = parse_args()

    fai_path = a.reads_fasta + '.fai'
    if not os.path.exists(fai_path):
        sys.exit(f'ERROR: missing FASTA index {fai_path} (run: samtools faidx {a.reads_fasta})')

    # ---- read -> assigned anchor (same significance filter the main script uses)
    df = pd.read_csv(a.blast_tsv, sep='\t')
    df = df[df['read_bp_used_for_match'] > df['total_anchor_length'] / 2]
    if df.empty:
        print('No anchored reads; nothing to check.')
        open(a.output, 'w').close()
        return
    best = df.sort_values('bitscore', ascending=False).drop_duplicates('read_id')
    assigned = dict(zip(best['read_id'], best['anchor_name']))
    by_anchor = collections.defaultdict(list)
    for rid, anc in assigned.items():
        by_anchor[anc].append(rid)
    print(f'anchored reads to screen: {len(assigned):,} across {len(by_anchor)} anchors')

    # ---- anchor sequences
    anchors, name, buf = {}, None, []
    with open(a.anchors_fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    anchors[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        anchors[name] = ''.join(buf)

    idx = load_fai(fai_path)
    foldbacks, detail = {}, []

    with open(a.reads_fasta, 'rb') as reads_fh, tempfile.TemporaryDirectory() as tmp:
        subj = os.path.join(tmp, 'anchor.fa')
        qry = os.path.join(tmp, 'reads.fa')
        for anchor_name in sorted(by_anchor):
            if anchor_name not in anchors:
                print(f'  WARNING: {anchor_name} not in anchors FASTA, skipping')
                continue
            rids = [r for r in by_anchor[anchor_name] if r in idx]
            if not rids:
                continue
            write_fasta(subj, [(anchor_name, anchors[anchor_name])])
            write_fasta(qry, [(r, fetch(reads_fh, idx[r])) for r in rids])

            # NOTE: -min_raw_gapped_score deliberately omitted -- that flag is
            # exactly what hides the degraded second copy.
            cmd = ['blastn', '-query', qry, '-subject', subj,
                   '-task', 'dc-megablast',
                   '-perc_identity', str(a.perc_identity),
                   '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send',
                   '-num_threads', str(a.threads)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                sys.exit(f'ERROR: blastn failed for {anchor_name}:\n{res.stderr[:500]}')

            hits = collections.defaultdict(list)
            for line in res.stdout.splitlines():
                f = line.split('\t')
                if len(f) < 8:
                    continue
                read, pid, ln = f[0], float(f[2]), int(f[3])
                if ln < a.min_len:
                    continue
                s_start, s_end = int(f[6]), int(f[7])
                hits[read].append({
                    'pident': pid, 'length': ln,
                    'read_start': int(f[4]), 'read_end': int(f[5]),
                    'forward': s_start < s_end,
                })

            n_here = 0
            for read, hs in hits.items():
                if len({h['forward'] for h in hs}) > 1:   # both orientations present
                    foldbacks[read] = anchor_name
                    n_here += 1

                    # Hairpin symmetry: a folded molecule is a palindrome about
                    # the fold point, so a feature d bp from the read start
                    # should reappear d bp from the read end. Measured on 115
                    # 7172 fold-backs the median offset is 162bp (0.23% of read
                    # length) and 80% fall within 1%. REPORTED ONLY, never used
                    # to gate: ~7% are genuine but truncated (the complement
                    # strand stopped early), and gating would discard those.
                    # A large asymmetry also hints the read is a ligation
                    # chimera of two molecules sharing an anchor rather than a
                    # true hairpin.
                    read_len = idx[read][0]
                    ordered = sorted(hs, key=lambda h: min(h['read_start'], h['read_end']))
                    d_start = min(ordered[0]['read_start'], ordered[0]['read_end'])
                    d_end = read_len - max(ordered[-1]['read_start'], ordered[-1]['read_end'])
                    asym = abs(d_start - d_end)
                    asym_pct = 100.0 * asym / max(1, read_len)

                    for h in sorted(hs, key=lambda x: -x['length']):
                        detail.append({
                            'read_id': read, 'anchor_name': anchor_name,
                            'read_length': read_len,
                            'pident': h['pident'], 'align_len': h['length'],
                            'read_start': h['read_start'], 'read_end': h['read_end'],
                            'orientation': 'forward' if h['forward'] else 'reverse',
                            'hairpin_asymmetry_bp': asym,
                            'hairpin_asymmetry_pct': round(asym_pct, 3),
                            'hairpin_confidence': ('high' if asym_pct <= 2
                                                   else 'medium' if asym_pct <= 10
                                                   else 'low'),
                        })
            print(f'  {anchor_name:<16} {len(rids):>6} reads  fold-backs: {n_here}')

    with open(a.output, 'w') as out:
        for r in sorted(foldbacks):
            out.write(r + '\n')
    if a.detail_tsv and detail:
        pd.DataFrame(detail).to_csv(a.detail_tsv, sep='\t', index=False)

    n = len(foldbacks)
    print(f'\nfold-back reads detected: {n} of {len(assigned):,} '
          f'({100 * n / max(1, len(assigned)):.2f}%)')
    print(f'wrote {a.output}')


if __name__ == '__main__':
    main()
