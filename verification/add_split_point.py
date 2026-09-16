#!/usr/bin/env python3
"""Add the Y' split point (anchor half / telomere half boundary) to the cut99 combined tables.

For a recombinant read, the junction divides the Y' into an anchor-side half and a
telomere-side half (already computed by flank_identities.py). This adds, for every row that
has a located junction:

  split_anchor_bp   length of the anchor-side half, in bp
  split_telo_bp     length of the telomere-side half, in bp
  split_total_bp     the two summed (approx. the read's Y' span)
  split_pct_from_anchor   position of the cut as % of the way through the Y', measured from
                          the anchor side (0% = right at the anchor, 100% = right at the
                          telomere). This is anchor-side-relative regardless of telo_side,
                          since anchor_half_bp/telo_half_bp are already biological-side-named.

Rows with no located junction (no junction / reference defect) get blank split columns.

Usage: add_split_point.py --bundle-root <dir with <strain>/ and <strain>/<sample>/ subdirs>
"""
import argparse, csv, glob, os, re

p = argparse.ArgumentParser()
p.add_argument('--bundle-root', required=True)
a = p.parse_args()
B = a.bundle_root

def load_flank(path):
    d = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        d[r['read_id']] = r
    return d

# gather every per-sample flank_identities.tsv, keyed by sample name (from the filename)
flanks = {}
for f in glob.glob(f'{B}/*/**/*_flank_identities.tsv', recursive=True) + glob.glob(f'{B}/*/*_flank_identities.tsv'):
    sample = os.path.basename(f).replace('_flank_identities.tsv', '')
    flanks[sample] = load_flank(f)

def split_cols(sample, read_id):
    fl = flanks.get(sample, {}).get(read_id)
    if not fl or not fl.get('anchor_half_bp') or not fl.get('telo_half_bp'):
        return {'split_anchor_bp': '', 'split_telo_bp': '', 'split_total_bp': '',
                'split_pct_from_anchor': ''}
    a_bp, t_bp = int(fl['anchor_half_bp']), int(fl['telo_half_bp'])
    tot = a_bp + t_bp
    pct = round(100 * a_bp / tot, 1) if tot else ''
    return {'split_anchor_bp': a_bp, 'split_telo_bp': t_bp, 'split_total_bp': tot,
            'split_pct_from_anchor': pct}

NEW = ['split_anchor_bp', 'split_telo_bp', 'split_total_bp', 'split_pct_from_anchor']

for tsv in glob.glob(f'{B}/*/ALL_*_combined.tsv'):
    rows = list(csv.DictReader(open(tsv), delimiter='\t'))
    if not rows: continue
    for r in rows:
        r.update(split_cols(r['sample'], r['read_id']))
    cols = [c for c in rows[0].keys() if c not in NEW]
    anchor_at = cols.index('donor_group') + 1
    cols = cols[:anchor_at] + NEW + cols[anchor_at:]
    with open(tsv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
        for r in rows: w.writerow(r)

    # rebuild the matching .md
    md = tsv.replace('.tsv', '.md')
    strain = os.path.basename(tsv).replace('ALL_', '').replace('_combined.tsv', '')
    lines = [f'# {strain}: every cut99 Y\' group mismatch', '',
             '| sample | read | end | expected | donor | split (anchor/telo bp) | cut at % from anchor | anchor margin | telo margin | evidence |',
             '|---|---|---|---|---|---|---|---|---|---|']
    order = {'strong': 0, 'weak': 1, 'FAILS': 2, 'no junction': 3, 'reference defect': 4}
    for r in sorted(rows, key=lambda x: (order.get(x['evidence'], 9), x['sample'], x['chr_end'])):
        split = f"{r['split_anchor_bp']}/{r['split_telo_bp']}" if r['split_anchor_bp'] else '—'
        pct = f"{r['split_pct_from_anchor']}%" if r['split_pct_from_anchor'] != '' else '—'
        samp = r['sample'].replace('6991_day0', '').replace('7172_day0_with_selection', '7172') \
                          .replace('7302_day0_with_selection', '7302') or 'day0'
        lines.append(f"| {samp} | {r['read_id']} | {r['chr_end']} | {r['expected_element']} | "
                     f"{r['donor_element']} | {split} | {pct} | {r['anchor_margin'] or '—'} | "
                     f"{r['telo_margin'] or '—'} | {r['evidence']} |")
    lines += ['', '`cut at % from anchor`: where the anchor-side half ends and the telomere-side',
              'half begins, as a percentage of the read\'s total Y\' span, measured from the',
              'anchor. 0% = the cut sits right at the anchor (an almost entirely foreign Y\');',
              '100% = right at the telomere (an almost entirely native Y\'). Located from the',
              'sliding-window scan (`scan_recombinant_junctions.py`), quantised to its 150 bp',
              'step, so treat it as approximate rather than base-pair precise. Blank for rows',
              'with no located junction (`no junction`, `reference defect`).', '']
    open(md, 'w').write('\n'.join(lines))
    print(f'  {tsv}  (+{len(NEW)} cols, {len(rows)} rows)')
    print(f'  {md}  rebuilt')
