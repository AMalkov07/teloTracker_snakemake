"""
Correct over-extended Y' element boundaries by measuring them against near-identical partners.

Why: the labelling step takes each Y' element's extent from a BLAST hit against a library, so
an element inherits whatever boundary the library entry happens to have. The known case is
chr16L_Y_Prime_1, which comes out 77 bp too long at its anchor-proximal end in every 6991,
7172 and 7302 day-0 reference. The extra 77 bp is X-element sequence (0 % telomeric repeat),
chr7R_1 / chr14L_1 / chr16L_1 are otherwise identical to within 1-3 bp, and the clustering
similarity penalises the unaligned overhang, so at a 99 % cut chr16L splits into its own group
and 806/806 chr7R reads get matched to the wrong element.

How: every Y' is extracted anchor -> telomere and BLASTed against every other Y'. Partners are
elements >= 99 % identical over >= 90 % of the longer length (which keeps Short Y' out of a Long
partner set). Where a partner alignment starts some bases into the query, those bases are an
anchor-proximal overhang the partner does not have. The element is trimmed only when at least
two partners agree on the overhang to within 5 bp and it is 20-300 bp long -- so a real
element whose partners disagree, or one with no partners, is never touched.

Elements clearly SHORTER than near-identical partners are flagged but never extended: that is
what a mis-assembled or truncated copy looks like (6991_day0_with_selection chr14L_Y_Prime_1 is
5,720 bp against 6,654 for its partners), and stretching it to consensus would hide the defect.

Generalises verification/fix_yprime_boundary.py, which measured the same thing for one element
and patched a BED after the fact. Here the correction is applied to the in-memory regions of
label_pretelomeric_regions.py before any file is written, so the BED, simplified BED (whose ITS
rows are derived from the gaps), GFF3 and TSV all agree by construction.

Coordinates follow label_pretelomeric_regions.py: 0-based half-open [start, end).
"""

import os
import subprocess
import tempfile
from statistics import median

MIN_IDENTITY = 99.0      # partner: percent identity of the best HSP
MIN_COVERAGE = 0.90      # partner: HSP length as a fraction of the LONGER element
MIN_PARTNERS = 2         # partners that must agree on the overhang
AGREE_BP = 5             # ... to within this many bp
MAX_TRIM = 300           # never trim more than this
MIN_TRIM = 20            # ... nor less. A 1-bp "overhang" is BLAST declining to extend through
                         # a mismatched terminal base, not an over-extension (chr7R_1 and
                         # chr15R_1 show exactly 1 bp on every day-0 reference). The known
                         # defect is 77 bp.
SHORT_MIN_DEFICIT = 500  # flag an element this much shorter than its full-coverage partners.
                         # Short-class Y' lengths vary naturally: the curated ID2 copies alone
                         # span 5,060-5,487 bp, and chr12L_1 is ~294 bp shorter than its
                         # partners in every independent 6991/7172/7302 assembly -- biology.
                         # The known mis-assembly (with_selection chr14L_1) is 932 bp short.
SHORT_MIN_QCOV = 0.95    # ... where the partner alignment covers this much of the element

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def _rc(seq):
    return seq.translate(_COMP)[::-1]


def yprime_elements(chr_end_regions, ref_sequences):
    """One record per Y' element, sequence oriented anchor -> telomere.

    Numbering matches write_bed_simplified: Y_Prime_1 is the copy nearest the anchor
    (L arm: highest coordinate; R arm: lowest).
    """
    out = []
    for chr_end, regions in chr_end_regions.items():
        yps = regions.get('y_prime', [])
        arm = chr_end[-1]
        order = sorted(range(len(yps)), key=lambda i: yps[i]['start'])
        for rank, i in enumerate(order):
            r = yps[i]
            seq = ref_sequences.get(r['chr'])
            if seq is None:
                continue
            num = len(order) - rank if arm == 'L' else rank + 1
            sub = seq[r['start']:r['end']]
            out.append({
                'name': f'{chr_end}_Y_Prime_{num}', 'chr_end': chr_end, 'index': i, 'arm': arm,
                'chr': r['chr'], 'start': r['start'], 'end': r['end'],
                'seq': _rc(sub) if arm == 'L' else sub,
            })
    return out


def blast_all_vs_all(elements, threads=4):
    """Best plus/plus HSP per (query, subject) pair, self hits excluded.

    Returns {(q, s): (pident, length, qlen, slen, qstart, qend, sstart, send)}.
    Only plus/plus: every element is already oriented anchor -> telomere, so a minus-strand
    HSP would be an inverted match, not a partner.
    """
    if len(elements) < 2:
        return {}
    with tempfile.TemporaryDirectory(prefix='ypbound_') as td:
        fa = os.path.join(td, 'yp.fasta')
        with open(fa, 'w') as fh:
            for e in elements:
                fh.write(f">{e['name']}\n{e['seq']}\n")
        db = os.path.join(td, 'db')
        subprocess.run(['makeblastdb', '-in', fa, '-dbtype', 'nucl', '-out', db],
                       check=True, capture_output=True)
        r = subprocess.run(['blastn', '-query', fa, '-db', db, '-evalue', '1e-10',
                            '-max_target_seqs', '1000', '-num_threads', str(threads),
                            '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend sstart send'],
                           check=True, capture_output=True, text=True)
    best = {}
    for line in r.stdout.splitlines():
        f = line.split('\t')
        q, s = f[0], f[1]
        if q == s:
            continue
        pid, ln, ql, sl, qs, qe, ss, se = float(f[2]), *map(int, f[3:10])
        if ss > se:
            continue
        if (q, s) not in best or ln > best[(q, s)][1]:
            best[(q, s)] = (pid, ln, ql, sl, qs, qe, ss, se)
    return best


def measure(elements, hits):
    """Decide a trim (and a short flag) for every element. Returns a list of dicts."""
    results = []
    for e in elements:
        q, qlen = e['name'], len(e['seq'])
        overhangs, short_partner_lens = [], []
        n_partners = 0
        for (qq, s), (pid, ln, ql, sl, qs, qe, ss, se) in hits.items():
            if qq != q or pid < MIN_IDENTITY:
                continue
            if ln >= MIN_COVERAGE * max(ql, sl):
                n_partners += 1
                # overhang = bases of the query before the shared start, but only when the
                # partner's own alignment starts at its position 1 -- otherwise the partner,
                # not the query, is the one carrying extra sequence
                if ss <= 1 + AGREE_BP:
                    overhangs.append(qs - 1)
                else:
                    overhangs.append(0)
            if ln >= SHORT_MIN_QCOV * ql and sl - ql >= SHORT_MIN_DEFICIT:
                short_partner_lens.append(sl)

        trim, agreeing = 0, 0
        if len(overhangs) >= MIN_PARTNERS:
            cand = sorted(overhangs)[len(overhangs) // 2]
            agreeing = sum(1 for o in overhangs if abs(o - cand) <= AGREE_BP)
            if MIN_TRIM <= cand <= MAX_TRIM and agreeing >= MIN_PARTNERS:
                trim = cand

        flag = ''
        if trim:
            flag = 'trimmed'
        elif short_partner_lens and len(short_partner_lens) >= MIN_PARTNERS:
            flag = f'shorter_than_partners_by_{int(median(short_partner_lens) - qlen)}bp'
        results.append({**{k: e[k] for k in ('name', 'chr_end', 'index', 'arm', 'chr', 'start', 'end')},
                        'length': qlen, 'n_partners': n_partners, 'n_agreeing': agreeing,
                        'trim': trim, 'new_length': qlen - trim, 'flag': flag})
    return results


def apply_trims(chr_end_regions, results):
    """Shrink each trimmed element at its anchor-proximal end, in place.

    L arm: the anchor side is the HIGH coordinate, so `end` moves down.
    R arm: the anchor side is the LOW coordinate, so `start` moves up.
    The freed bases fall into the gap before the element, which write_bed_simplified reports
    as ITS_<end>_Y_Prime_(n-1)-n (ITS_0-1 for the first copy).
    """
    for r in results:
        if not r['trim']:
            continue
        region = chr_end_regions[r['chr_end']]['y_prime'][r['index']]
        if r['arm'] == 'L':
            region['end'] -= r['trim']
        else:
            region['start'] += r['trim']
        region['boundary_trim'] = r['trim']
    return chr_end_regions


def write_provenance(results, path):
    cols = ['name', 'chr', 'start', 'end', 'length', 'n_partners', 'n_agreeing', 'trim',
            'new_length', 'flag']
    with open(path, 'w') as fh:
        fh.write('# Y\' boundary check: partners are >= %.0f %% identical over >= %.0f %% of the '
                 'longer element; a trim needs >= %d partners agreeing within %d bp, %d-%d bp.\n'
                 % (MIN_IDENTITY, MIN_COVERAGE * 100, MIN_PARTNERS, AGREE_BP, MIN_TRIM, MAX_TRIM))
        fh.write('# start/end are the ORIGINAL 0-based half-open coordinates; lengths in bp.\n')
        fh.write('\t'.join(cols) + '\n')
        for r in sorted(results, key=lambda x: x['name']):
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')


def correct_yprime_boundaries(chr_end_regions, ref_sequences, provenance_path=None, threads=4):
    """Measure, apply and report. Returns (chr_end_regions, results)."""
    elements = yprime_elements(chr_end_regions, ref_sequences)
    results = measure(elements, blast_all_vs_all(elements, threads=threads))
    apply_trims(chr_end_regions, results)
    if provenance_path:
        write_provenance(results, provenance_path)
    return chr_end_regions, results
