#!/usr/bin/env python3
"""
verify_day0_reference.py -- Part A verification of a pipeline-built day-0
reference against a curated (hand-checked) reference for the same strain.

Four checks, each written to its own TSV plus a combined Markdown report:

  Depth 1  per-end Y' counts            <prefix>_counts.tsv
  Depth 2  per-feature coordinates      <prefix>_features.tsv
           (anchor-relative offsets + lengths, so two assemblies with
            different absolute coordinates can be compared)
  Depth 3  Y' grouping concordance      <prefix>_yprime_assignment.tsv
           (label-independent: Adjusted Rand Index + per-element mapping +
            BLAST of our unique variants against the curated variants)
  Prov.    library provenance           <prefix>_provenance.tsv
           (which Y' library the run used; BED <-> library bijection)

Usage:
  python verify_day0_reference.py \
      --ours-bed     .../pretelomeric_regions_<ref>_simp.bed \
      --ours-lib     .../extracted_yprimes_<ref>.fasta \
      --curated-bed  .../7302_features/7302_final_features.bed \
      --curated-lib  .../7302_features/repeatmasker_7302_all_y_primes.fasta \
      --strain 7302 --ref-name 7302_day0_with_selection \
      --out-prefix verification/reports/partA/7302_day0_with_selection \
      [--probe-blast ..._probe_blast.txt] [--working-lib ..._working_yprimes.fasta] \
      [--run-config <run-dir _pipeline/config.yaml>] [--summary-tsv partA_summary.tsv]
"""

import argparse
import hashlib
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from recombination_utils import read_fasta, run_blast
from analyze_features import parse_y_prime_header

# ---------------------------------------------------------------------------
# Tolerances (bp). Offsets accumulate small ITS-length differences along a
# tandem array, so the offset tolerance is looser than the length tolerance.
# ---------------------------------------------------------------------------
LENGTH_TOL = {
    'y_prime': ('pct', 2.0, 50),   # max(2 %, 50 bp)
    'x_core': ('pct', 2.0, 50),
    'x_variable': ('pct', 2.0, 50),
    'its': ('abs', 25, 25),
    'spacer': ('abs', 100, 100),
    'anchor': ('abs', 5, 5),
    'telomere': ('abs', 10**9, 10**9),   # informational only
}
OFFSET_TOL = 100
MAJOR_BP = 500        # a length/offset difference above this is a real labeling/assembly defect
ARI_PASS = 0.85
CHR_END_RE = re.compile(r'(chr\d+[LR])')


# ---------------------------------------------------------------------------
# BED loading (both files are 6-column, 1-based inclusive, same vocabulary)
# ---------------------------------------------------------------------------

def classify(name):
    if 'Telomere_Repeat' in name:
        return 'telomere'
    if name.startswith('ITS_'):
        return 'its'
    if '_x_variable_element' in name:
        return 'x_variable'
    if '_x_core_element' in name:
        return 'x_core'
    if '_space_between_anchor' in name:
        return 'spacer'
    if '_Y_Prime_' in name:
        return 'y_prime'
    if name.endswith('_anchor'):
        return 'anchor'
    return 'unknown'


def feature_key(name, chr_end):
    """Strip the chr_end prefix so the same feature has the same key in both
    files: 'chr4R_Y_Prime_3' -> 'Y_Prime_3', 'ITS_chr12_Y_Prime_0-1' ->
    'ITS_Y_Prime_0-1'."""
    if name.startswith('ITS_'):
        return re.sub(r'^ITS_chr\d+[LR]?_', 'ITS_', name)
    if chr_end and name.startswith(chr_end + '_'):
        return name[len(chr_end) + 1:]
    return name


def load_bed(path):
    feats = []
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        if len(p) < 4 or line.startswith('#'):
            continue
        name = p[3]
        m = CHR_END_RE.search(name)
        feats.append({
            'contig': p[0], 'start': int(p[1]), 'end': int(p[2]), 'name': name,
            'strand': p[4] if len(p) > 4 else '',
            'length': int(p[5]) if len(p) > 5 and p[5].isdigit() else int(p[2]) - int(p[1]) + 1,
            'chr_end': m.group(1) if m else None, 'ftype': classify(name),
        })
    # Resolve arm-less names (curated 'ITS_chr12_Y_Prime_0-1') by the nearest
    # Y' feature on the same contig.
    for f in feats:
        if f['chr_end'] is None:
            same = [g for g in feats if g['contig'] == f['contig'] and g['chr_end'] and g['ftype'] == 'y_prime']
            if same:
                f['chr_end'] = min(same, key=lambda g: abs(g['start'] - f['start']))['chr_end']
    for f in feats:
        f['key'] = feature_key(f['name'], f['chr_end'])
    return feats


def by_end(feats):
    d = defaultdict(list)
    for f in feats:
        if f['chr_end']:
            d[f['chr_end']].append(f)
    return d


def end_sort_key(ce):
    m = re.match(r'chr(\d+)([LR])', ce)
    return (int(m.group(1)), m.group(2)) if m else (99, ce)


# ---------------------------------------------------------------------------
# Y' library parsing  (header: >Y_Prime_chr13L1,3;chr14L3,4,5#Short/Tandem/ID2_Red)
# ---------------------------------------------------------------------------

def parse_origin_locations(origin):
    out = []
    body = origin.replace('Y_Prime_', '', 1)
    for grp in body.split(';'):
        m = re.match(r'(chr\d+[LR])([\d,]+)', grp.strip())
        if not m:
            continue
        for p in m.group(2).split(','):
            if p:
                out.append((m.group(1), int(p)))
    return out


def parse_lib(path, level='family'):
    """Return ({(chr_end,pos): rec}, {header: rec}, {header: seq}).

    level='family'  -> rec['id'] is the ID number   (ID2)
    level='variant' -> rec['id'] is the full group   (ID2_Red-Light)
    The curated libraries carry two levels: the ID number is a sequence family
    and the colour shade (-Light/-Dark/-Neutral) a variant within it. The
    pipeline's own libraries have a single level (colour is a function of ID)."""
    seqs = read_fasta(path)
    elem, recs = {}, {}
    for header in seqs:
        info = parse_y_prime_header(header)
        cls = header.split('#', 1)[1] if '#' in header else ''
        size = cls.split('/')[0] if cls else ''
        gid = info['color_group'] if level == 'variant' else info['id']
        rec = {'header': header, 'id': gid, 'family': info['id'], 'color': info['color_group'], 'size': size,
               'locations': parse_origin_locations(info['origin'])}
        recs[header] = rec
        for loc in rec['locations']:
            elem[loc] = rec
    return elem, recs, seqs


# ---------------------------------------------------------------------------
# Depth 1 -- counts
# ---------------------------------------------------------------------------

def yprime_counts(feats):
    c = Counter()
    for f in feats:
        if f['ftype'] == 'y_prime':
            c[f['chr_end']] += 1
    return c


def probe_counts(probe_blast):
    if not probe_blast or not os.path.exists(probe_blast):
        return {}
    try:
        from label_pretelomeric_regions import count_yprimes_from_probe
        return count_yprimes_from_probe(probe_blast, '')
    except Exception as e:  # pragma: no cover
        print(f'  (probe count unavailable: {e})', file=sys.stderr)
        return {}


def depth1(cur, ours, probe):
    cc, oc = yprime_counts(cur), yprime_counts(ours)
    ends = sorted(set(cc) | set(oc) | set(probe), key=end_sort_key)
    rows = []
    for ce in ends:
        n_c, n_o, n_p = cc.get(ce, 0), oc.get(ce, 0), probe.get(ce, '')
        rows.append({'chr_end': ce, 'n_curated': n_c, 'n_ours': n_o, 'n_probe_ours': n_p,
                     'diff_ours_minus_curated': n_o - n_c,
                     'verdict': 'match' if n_c == n_o else 'MISMATCH'})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Depth 2 -- anchor-relative feature coordinates
# ---------------------------------------------------------------------------

def anchor_relative(feats_one_end, chr_end):
    anchors = [f for f in feats_one_end if f['ftype'] == 'anchor']
    if not anchors:
        return None
    a = anchors[0]
    arm = chr_end[-1]
    out = {}
    for f in feats_one_end:
        off = (f['start'] - a['end']) if arm == 'R' else (a['start'] - f['end'])
        out[f['key']] = {'offset': off, 'length': f['length'], 'ftype': f['ftype'], 'name': f['name']}
    return out


def within(ftype, dlen):
    mode, val, floor = LENGTH_TOL.get(ftype, ('abs', 50, 50))
    return abs(dlen) <= (max(floor, val) if mode == 'abs' else floor)


def depth2(cur_by_end, ours_by_end):
    rows = []
    for ce in sorted(set(cur_by_end) | set(ours_by_end), key=end_sort_key):
        c = anchor_relative(cur_by_end.get(ce, []), ce)
        o = anchor_relative(ours_by_end.get(ce, []), ce)
        if c is None or o is None:
            rows.append({'chr_end': ce, 'feature': '(anchor)', 'ftype': 'anchor', 'verdict':
                         'missing_anchor_in_' + ('ours' if o is None else 'curated'), 'severity': 'major'})
            continue
        for key in sorted(set(c) | set(o), key=lambda k: (c.get(k) or o.get(k))['offset']):
            fc, fo = c.get(key), o.get(key)
            ftype = (fc or fo)['ftype']
            if fc and fo:
                dlen = fo['length'] - fc['length']
                doff = fo['offset'] - fc['offset']
                if ftype in ('y_prime', 'x_core', 'x_variable'):
                    len_ok = abs(dlen) <= max(50, 0.02 * fc['length'])
                else:
                    len_ok = within(ftype, dlen)
                off_ok = abs(doff) <= OFFSET_TOL or ftype == 'anchor'
                verdict = 'ok' if (len_ok and off_ok) else ('length_diff' if not len_ok and off_ok
                          else 'offset_diff' if len_ok else 'length+offset_diff')
                severity = '' if verdict == 'ok' else ('major' if (abs(dlen) > MAJOR_BP or abs(doff) > MAJOR_BP) else 'minor')
                rows.append({'chr_end': ce, 'feature': key, 'ftype': ftype,
                             'length_curated': fc['length'], 'length_ours': fo['length'], 'dlen': dlen,
                             'offset_curated': fc['offset'], 'offset_ours': fo['offset'], 'doffset': doff,
                             'verdict': verdict, 'severity': severity})
            elif fc:
                v = 'expected_absent_in_ours' if ftype == 'telomere' else 'missing_in_ours'
                rows.append({'chr_end': ce, 'feature': key, 'ftype': ftype,
                             'length_curated': fc['length'], 'offset_curated': fc['offset'],
                             'verdict': v, 'severity': '' if ftype == 'telomere' else
                             ('minor' if (ftype == 'its' and fc['length'] < 100) else 'major')})
            else:
                rows.append({'chr_end': ce, 'feature': key, 'ftype': ftype,
                             'length_ours': fo['length'], 'offset_ours': fo['offset'],
                             'verdict': 'missing_in_curated',
                             'severity': 'minor' if (ftype == 'its' and fo['length'] < 100) else 'major'})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Depth 3 -- grouping concordance
# ---------------------------------------------------------------------------

def depth3(our_elem, cur_elem, cur_elem_fam, our_lib, cur_lib, cur_lib_path, our_seqs, working_lib):
    """cur_elem is at the finer 'variant' level (ID2_Red-Light); cur_elem_fam at
    the 'family' level (ID2). Verdicts and the mapping are computed at the
    variant level; both ARIs are reported."""
    common = sorted(set(our_elem) & set(cur_elem), key=lambda e: (end_sort_key(e[0]), e[1]))
    only_ours = sorted(set(our_elem) - set(cur_elem))
    only_cur = sorted(set(cur_elem) - set(our_elem))

    our_lab = [our_elem[e]['id'] for e in common]
    cur_lab = [cur_elem[e]['id'] for e in common]
    fam_lab = [cur_elem_fam[e]['id'] for e in common]
    ari = ari_fam = float('nan')
    if common:
        from sklearn.metrics import adjusted_rand_score
        ari = adjusted_rand_score(cur_lab, our_lab)
        ari_fam = adjusted_rand_score(fam_lab, our_lab)

    our_groups = defaultdict(set)
    cur_groups = defaultdict(set)
    for e, o, c in zip(common, our_lab, cur_lab):
        our_groups[o].add(e)
        cur_groups[c].add(e)

    mapping = {}
    for o, members in our_groups.items():
        ov = Counter(cur_elem[e]['id'] for e in members)
        best, n = ov.most_common(1)[0]
        mapping[o] = {'curated_id': best, 'overlap': n, 'size': len(members),
                      'purity': round(n / len(members), 3),
                      'curated_ids_seen': ','.join(f'{k}:{v}' for k, v in ov.most_common())}

    # inherited (pre-clustering) IDs -- labeling BLASTs every Y' against the
    # 6991 curated library, so these are the curated-scheme IDs by best hit
    inherited = {}
    if working_lib and os.path.exists(working_lib):
        w_elem, _, _ = parse_lib(working_lib, level='variant')
        inherited = {e: r['id'] for e, r in w_elem.items()}

    rows = []
    for e, o, c in zip(common, our_lab, cur_lab):
        O, C = our_groups[o], cur_groups[c]
        if O == C:
            verdict = 'concordant'
        elif C < O:
            verdict = 'split_in_curated'
        elif O < C:
            verdict = 'merged_in_curated'
        else:
            verdict = 'CONFLICT'
        inh = inherited.get(e, '')
        rows.append({'chr_end': e[0], 'pos': e[1], 'element': f'{e[0]}_Y_Prime_{e[1]}',
                     'our_id': o, 'our_size': our_elem[e]['size'], 'our_group_n': len(O),
                     'curated_variant': c, 'curated_family': cur_elem_fam[e]['id'],
                     'curated_size': cur_elem[e]['size'], 'curated_group_n': len(C),
                     'size_class_agrees': our_elem[e]['size'] == cur_elem[e]['size'],
                     'inherited_variant_preclustering': inh,
                     'inherited_matches_curated': (inh == c) if inh else '',
                     'verdict': verdict})
    for e in only_ours:
        rows.append({'chr_end': e[0], 'pos': e[1], 'element': f'{e[0]}_Y_Prime_{e[1]}',
                     'our_id': our_elem[e]['id'], 'our_size': our_elem[e]['size'], 'verdict': 'only_in_ours'})
    for e in only_cur:
        rows.append({'chr_end': e[0], 'pos': e[1], 'element': f'{e[0]}_Y_Prime_{e[1]}',
                     'curated_variant': cur_elem[e]['id'], 'curated_size': cur_elem[e]['size'], 'verdict': 'only_in_curated'})
    per_elem = pd.DataFrame(rows)

    contingency = pd.crosstab(pd.Series(our_lab, name='our_id'), pd.Series(cur_lab, name='curated_variant')) if common else pd.DataFrame()

    # How different are the curated variants that our clustering merged?
    merged_rows = []
    with tempfile.TemporaryDirectory() as tmp:
        cur_seqs = read_fasta(cur_lib_path)
        rep = {}
        for h in cur_seqs:
            rep.setdefault(cur_lib[h]['id'], h)
        q = os.path.join(tmp, 'cur_reps.fa')
        with open(q, 'w') as fh:
            for v, h in rep.items():
                fh.write(f'>{h.split()[0]}\n{cur_seqs[h]}\n')
        df = run_blast(q, cur_lib_path, tmp, label='curated_all_vs_all', min_pident=75.0)
        hdr_to_var = {h.split()[0]: cur_lib[h]['id'] for h in cur_seqs}
        best = {}
        for _, r in df.iterrows():
            a, b = hdr_to_var.get(r['qseqid']), hdr_to_var.get(r['sseqid'])
            if a is None or b is None or a == b:
                continue
            key = tuple(sorted((a, b)))
            if key not in best or r['bitscore'] > best[key]['bitscore']:
                best[key] = {'bitscore': r['bitscore'], 'pident': r['pident'],
                             'coverage_pct': round(100.0 * r['length'] / r['qlen'], 1)}
        for o, members in our_groups.items():
            variants = sorted({cur_elem[e]['id'] for e in members})
            for i in range(len(variants)):
                for j in range(i + 1, len(variants)):
                    key = (variants[i], variants[j])
                    b = best.get(key, {})
                    merged_rows.append({'our_id': o, 'curated_variant_a': variants[i], 'curated_variant_b': variants[j],
                                        'pident': b.get('pident', ''), 'coverage_pct': b.get('coverage_pct', ''),
                                        'size_a': next(cur_elem[e]['size'] for e in members if cur_elem[e]['id'] == variants[i]),
                                        'size_b': next(cur_elem[e]['size'] for e in members if cur_elem[e]['id'] == variants[j])})
    merged_df = pd.DataFrame(merged_rows)

    # sequence confirmation: our unique variants vs curated variants
    seq_rows = []
    with tempfile.TemporaryDirectory() as tmp:
        q = os.path.join(tmp, 'ours.fa')
        with open(q, 'w') as fh:
            for i, (h, s) in enumerate(our_seqs.items()):
                fh.write(f'>q{i}\n{s}\n')
        idx = {f'q{i}': h for i, h in enumerate(our_seqs)}
        df = run_blast(q, cur_lib_path, tmp, label='ours_vs_curated', min_pident=75.0)
        for qid, h in idx.items():
            sub = df[df['qseqid'] == qid]
            rec = our_lib[h]
            if sub.empty:
                seq_rows.append({'our_header': h, 'our_id': rec['id'], 'our_size': rec['size'], 'best_hit': '', 'verdict': 'NO_HIT'})
                continue
            b = sub.sort_values('bitscore', ascending=False).iloc[0]
            cov = round(100.0 * b['length'] / b['qlen'], 1)
            hit_hdr = b['sseqid']
            cur_rec = next((r for hh, r in cur_lib.items() if hh.split()[0] == hit_hdr), None)
            cur_id = cur_rec['id'] if cur_rec else ''
            cur_size = cur_rec['size'] if cur_rec else ''
            expected = mapping.get(rec['id'], {}).get('curated_id', '')
            seq_rows.append({'our_header': h, 'our_id': rec['id'], 'our_size': rec['size'],
                             'best_hit': hit_hdr, 'hit_curated_id': cur_id, 'hit_size': cur_size,
                             'pident': b['pident'], 'coverage_pct': cov, 'bitscore': b['bitscore'],
                             'mapped_curated_id': expected,
                             'verdict': 'ok' if (cur_id == expected or not expected) else 'hit_id_differs_from_mapping'})
    seq_df = pd.DataFrame(seq_rows)
    return {'ari': ari, 'ari_family': ari_fam, 'n_common': len(common), 'per_elem': per_elem,
            'contingency': contingency, 'mapping': mapping, 'seq': seq_df, 'merged': merged_df,
            'n_our_clusters': len(our_groups), 'n_curated_variants': len(cur_groups),
            'n_curated_families': len(set(fam_lab)),
            'n_conflict': int((per_elem['verdict'] == 'CONFLICT').sum()) if not per_elem.empty else 0,
            'n_split': int((per_elem['verdict'] == 'split_in_curated').sum()) if not per_elem.empty else 0,
            'n_merged': int((per_elem['verdict'] == 'merged_in_curated').sum()) if not per_elem.empty else 0,
            'n_concordant': int((per_elem['verdict'] == 'concordant').sum()) if not per_elem.empty else 0,
            'inherited_frac': (per_elem['inherited_matches_curated'].replace('', pd.NA).dropna().astype(bool).mean()
                               if not per_elem.empty and 'inherited_matches_curated' in per_elem else float('nan'))}


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def parse_run_config(path):
    cfg = {'strain': '', 'y_prime_lib': '', 'y_prime_lib_override': ''}
    if not path or not os.path.exists(path):
        return cfg
    for line in open(path):
        s = line.split('#', 1)[0].strip()
        for k in cfg:
            m = re.match(rf'^{k}\s*:\s*"?([^"]*)"?\s*$', s)
            if m:
                cfg[k] = m.group(1).strip()
    return cfg


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def provenance(run_config, ours_lib, ours_bed_feats, our_elem):
    cfg = parse_run_config(run_config)
    resolved = cfg['y_prime_lib_override'] or cfg['y_prime_lib'].replace('{strain}', cfg['strain'])
    bed_elems = set()
    for f in ours_bed_feats:
        if f['ftype'] == 'y_prime':
            m = re.match(r'(chr\d+[LR])_Y_Prime_(\d+)$', f['name'])
            if m:
                bed_elems.add((m.group(1), int(m.group(2))))
    lib_elems = set(our_elem)
    rows = [
        {'check': 'run_config', 'value': run_config or '(not provided)', 'ok': bool(cfg['strain'])},
        {'check': 'config_strain', 'value': cfg['strain'], 'ok': True},
        {'check': 'y_prime_lib_override_present', 'value': bool(cfg['y_prime_lib_override']),
         'ok': not cfg['y_prime_lib_override']},
        {'check': 'y_prime_lib_resolved', 'value': resolved,
         'ok': (os.path.basename(resolved) == os.path.basename(ours_lib)) if resolved else False},
        {'check': 'ours_lib_md5', 'value': md5(ours_lib), 'ok': True},
        {'check': 'bed_yprimes_not_in_lib', 'value': ','.join(f'{c}_{p}' for c, p in sorted(bed_elems - lib_elems)) or '-',
         'ok': not (bed_elems - lib_elems)},
        {'check': 'lib_yprimes_not_in_bed', 'value': ','.join(f'{c}_{p}' for c, p in sorted(lib_elems - bed_elems)) or '-',
         'ok': not (lib_elems - bed_elems)},
        {'check': 'n_bed_yprimes', 'value': len(bed_elems), 'ok': True},
        {'check': 'n_lib_yprime_locations', 'value': len(lib_elems), 'ok': True},
        {'check': 'n_lib_unique_variants', 'value': len({r['header'] for r in our_elem.values()}), 'ok': True},
    ]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def md_table(df, max_rows=60):
    if df is None or df.empty:
        return '_(none)_\n'
    d = df.head(max_rows)
    cols = list(d.columns)
    out = ['| ' + ' | '.join(str(c) for c in cols) + ' |', '|' + '---|' * len(cols)]
    for _, r in d.iterrows():
        out.append('| ' + ' | '.join('' if pd.isna(v) else str(v) for v in r.values) + ' |')
    if len(df) > max_rows:
        out.append(f'\n_({len(df) - max_rows} more rows in the TSV)_')
    return '\n'.join(out) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ours-bed', required=True)
    ap.add_argument('--ours-lib', required=True)
    ap.add_argument('--curated-bed', required=True)
    ap.add_argument('--curated-lib', required=True)
    ap.add_argument('--strain', required=True)
    ap.add_argument('--ref-name', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--probe-blast', default='')
    ap.add_argument('--working-lib', default='', help='pre-clustering *_working_yprimes.fasta (inherited IDs)')
    ap.add_argument('--run-config', default='')
    ap.add_argument('--summary-tsv', default='', help='append a one-line verdict row here')
    a = ap.parse_args()

    os.makedirs(os.path.dirname(a.out_prefix) or '.', exist_ok=True)
    cur = load_bed(a.curated_bed)
    ours = load_bed(a.ours_bed)
    our_elem, our_lib, our_seqs = parse_lib(a.ours_lib)
    cur_elem, cur_lib, _ = parse_lib(a.curated_lib, level='variant')
    cur_elem_fam, _, _ = parse_lib(a.curated_lib, level='family')

    d1 = depth1(cur, ours, probe_counts(a.probe_blast))
    d2 = depth2(by_end(cur), by_end(ours))
    d3 = depth3(our_elem, cur_elem, cur_elem_fam, our_lib, cur_lib, a.curated_lib, our_seqs, a.working_lib)
    pv = provenance(a.run_config, a.ours_lib, ours, our_elem)

    d1.to_csv(f'{a.out_prefix}_counts.tsv', sep='\t', index=False)
    d2.to_csv(f'{a.out_prefix}_features.tsv', sep='\t', index=False)
    d3['per_elem'].to_csv(f'{a.out_prefix}_yprime_assignment.tsv', sep='\t', index=False)
    d3['seq'].to_csv(f'{a.out_prefix}_yprime_sequence_check.tsv', sep='\t', index=False)
    d3['merged'].to_csv(f'{a.out_prefix}_yprime_merged_variants.tsv', sep='\t', index=False)
    pv.to_csv(f'{a.out_prefix}_provenance.tsv', sep='\t', index=False)

    n_count_mm = int((d1['verdict'] == 'MISMATCH').sum())
    bad_feats = d2[~d2['verdict'].isin(['ok', 'expected_absent_in_ours'])]
    n_feat_bad = len(bad_feats)
    n_feat_major = int((bad_feats['severity'] == 'major').sum()) if not bad_feats.empty else 0
    seq_bad = d3['seq'][d3['seq']['verdict'] != 'ok'] if not d3['seq'].empty else d3['seq']
    prov_ok = bool(pv['ok'].all())
    a1 = 'PASS' if n_count_mm == 0 else 'FAIL'
    a2 = 'PASS' if n_feat_major == 0 else 'FAIL'
    a3 = 'PASS' if (d3['ari'] >= ARI_PASS and d3['n_conflict'] == 0) else 'FAIL'
    a4 = 'PASS' if prov_ok else 'FAIL'

    L = [f'# Day-0 reference verification: {a.ref_name} (strain {a.strain})\n',
         f'- curated bed: `{a.curated_bed}`', f'- curated lib: `{a.curated_lib}`',
         f'- our bed: `{a.ours_bed}`', f'- our lib: `{a.ours_lib}`\n',
         '## Verdicts\n',
         '| check | result | detail |', '|---|---|---|',
         f"| A1 per-end Y' counts | **{a1}** | {n_count_mm} end(s) differ |",
         f'| A2 feature coordinates | **{a2}** | {n_feat_major} major, {n_feat_bad - n_feat_major} minor difference(s) |',
         f"| A3 Y' grouping | **{a3}** | ARI vs curated variants={d3['ari']:.3f}, vs curated families={d3['ari_family']:.3f} "
         f"over {d3['n_common']} elements ({d3['n_our_clusters']} our clusters vs {d3['n_curated_variants']} curated variants / "
         f"{d3['n_curated_families']} families); {d3['n_concordant']} concordant, {d3['n_split']} split_in_curated, "
         f"{d3['n_merged']} merged_in_curated, {d3['n_conflict']} CONFLICT |",
         f'| A4 library provenance | **{a4}** | {int((~pv["ok"]).sum())} failing check(s) |\n',
         "## A1 -- per-end Y' counts (mismatches only)\n", md_table(d1[d1['verdict'] == 'MISMATCH']),
         '\nFull table: `' + os.path.basename(a.out_prefix) + '_counts.tsv`\n',
         '## A2 -- features outside tolerance or missing\n', md_table(bad_feats),
         "\n## A3 -- Y' grouping concordance\n",
         f"Adjusted Rand Index vs curated **variants** (ID + colour shade) = **{d3['ari']:.3f}**; "
         f"vs curated **families** (ID number only) = **{d3['ari_family']:.3f}** "
         f"(1.0 = identical partition; pass >= {ARI_PASS} at the variant level).\n",
         "Verdicts: `split_in_curated` = curated is finer (our cluster = union of several curated variants); "
         "`merged_in_curated` = ours is finer; `CONFLICT` = the element is grouped with different partners in the two schemes.\n",
         '### Contingency (rows = our ID, cols = curated variant)\n',
         md_table(d3['contingency'].reset_index()) if not d3['contingency'].empty else '_(none)_\n',
         '### Cluster mapping ours -> curated variant\n',
         md_table(pd.DataFrame([{'our_id': k, **v} for k, v in sorted(d3['mapping'].items())])),
         '### Curated variants that our clustering merged (pairwise identity of representatives)\n',
         md_table(d3['merged']),
         '### Elements not concordant\n',
         md_table(d3['per_elem'][d3['per_elem']['verdict'] != 'concordant']),
         f"\nInherited (pre-clustering, 6991-curated-scheme) ID equals this strain's curated ID for "
         f"{d3['inherited_frac']*100:.0f}% of elements.\n" if d3['inherited_frac'] == d3['inherited_frac'] else '',
         '### Sequence confirmation (our unique variants BLASTed vs curated variants)\n',
         (f"min pident {d3['seq']['pident'].min():.2f}, min coverage {d3['seq']['coverage_pct'].min():.1f}%; "
          f"{len(seq_bad)} variant(s) whose best hit disagrees with the cluster mapping\n") if not d3['seq'].empty and 'pident' in d3['seq'] else '',
         md_table(seq_bad),
         '\n## A4 -- library provenance\n', md_table(pv)]
    with open(f'{a.out_prefix}_day0_verification.md', 'w') as fh:
        fh.write('\n'.join(x for x in L if x is not None))

    if a.summary_tsv:
        row = {'ref': a.ref_name, 'strain': a.strain, 'A1_counts': a1, 'n_count_mismatch': n_count_mm,
               'A2_features': a2, 'n_feature_major': n_feat_major, 'n_feature_minor': n_feat_bad - n_feat_major,
               'A3_grouping': a3, 'ari_variant': round(d3['ari'], 3), 'ari_family': round(d3['ari_family'], 3),
               'n_our_clusters': d3['n_our_clusters'], 'n_curated_variants': d3['n_curated_variants'],
               'n_common_elements': d3['n_common'], 'n_concordant': d3['n_concordant'], 'n_split_in_curated': d3['n_split'],
               'n_merged_in_curated': d3['n_merged'], 'n_conflict': d3['n_conflict'],
               'seq_min_pident': round(d3['seq']['pident'].min(), 2) if 'pident' in d3['seq'] else '',
               'inherited_id_match_frac': round(d3['inherited_frac'], 3) if d3['inherited_frac'] == d3['inherited_frac'] else '',
               'A4_provenance': a4}
        hdr = not os.path.exists(a.summary_tsv)
        pd.DataFrame([row]).to_csv(a.summary_tsv, sep='\t', index=False, mode='a', header=hdr)

    print(f'{a.ref_name}: A1={a1} ({n_count_mm} ends)  A2={a2} ({n_feat_major} major/{n_feat_bad - n_feat_major} minor)  '
          f"A3={a3} (ARI variant={d3['ari']:.3f} family={d3['ari_family']:.3f}, conflicts={d3['n_conflict']})  A4={a4}")


if __name__ == '__main__':
    main()
