"""
Build a single self-contained HTML report for one sample (or a comparison across samples).

The report is one file with no external assets: charts are inline SVG, sorting is inline JS.
That matters because these runs live on HPC scratch -- the file can be opened directly over
VSCode Remote / SSHFS, or scp'd anywhere, with no network access and no plotting libraries.

Usage:
  # one sample
  python make_sample_report.py --pipeline-dir results/{base}/_pipeline --base-name {base} \
      --output results/{base}/_pipeline/{base}_report.html

  # several samples side by side (also writes each sample's own report next to it)
  python make_sample_report.py --results-dir results \
      --samples SAMPLE_A SAMPLE_B SAMPLE_C --compare-output results/panel_report.html
"""

import argparse
import glob
import html
import math
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd


# --------------------------------------------------------------------------------------
# tolerant parsing helpers
# --------------------------------------------------------------------------------------

def to_num(value):
    """Return float(value) or None. Tolerates '', 'NA', '0.0' in int-looking columns."""
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(out) else out


def to_int(value):
    """int() that survives '0.0' -- the trap that silently empties a join."""
    num = to_num(value)
    return None if num is None else int(num)


def read_tsv(path):
    """Read a pipeline TSV. Skips leading '#' comment lines (read_summary.tsv has five)
    without using pandas' comment= , which would also truncate any field containing '#'
    -- Y' library headers legitimately do (e.g. 'chr14L3,4#Short/Tandem/ID2_Red')."""
    if not path or not os.path.isfile(path):
        return None
    try:
        skip = 0
        with open(path) as fh:
            for line in fh:
                if line.startswith('#'):
                    skip += 1
                else:
                    break
        return pd.read_csv(path, sep='\t', skiprows=skip, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return None                          # expected for ends the pipeline skipped
    except Exception as exc:                                    # noqa: BLE001
        print(f'  WARNING: could not read {path}: {exc}')
        return None


def pct(part, whole):
    return 0.0 if not whole else 100.0 * part / whole


def quantile(sorted_values, q):
    if not sorted_values:
        return None
    idx = min(int(len(sorted_values) * q), len(sorted_values) - 1)
    return sorted_values[idx]


def describe(values):
    """Summary stats for a numeric list."""
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return {}
    total = sum(vals)
    return {
        'n': len(vals),
        'min': vals[0],
        'p10': quantile(vals, 0.10),
        'median': quantile(vals, 0.50),
        'mean': total / len(vals),
        'p90': quantile(vals, 0.90),
        'p99': quantile(vals, 0.99),
        'max': vals[-1],
    }


# --------------------------------------------------------------------------------------
# inline SVG charts -- no external plotting dependency
# --------------------------------------------------------------------------------------

def _svg_open(width, height, title=''):
    lab = f'<title>{html.escape(title)}</title>' if title else ''
    return (f'<svg class="chart" viewBox="0 0 {width} {height}" width="100%" '
            f'preserveAspectRatio="xMidYMid meet" role="img">{lab}')


def svg_barh(rows, value_max=None, unit='', width=680, row_h=17, label_w=58, color='var(--accent)'):
    """Horizontal bars: rows = [(label, value, optional_note), ...]."""
    rows = list(rows)
    if not rows:
        return '<p class="muted">no data</p>'
    if value_max is None:
        value_max = max((r[1] for r in rows if r[1] is not None), default=0)
    value_max = value_max or 1
    bar_w = width - label_w - 86
    height = len(rows) * row_h + 10
    out = [_svg_open(width, height)]
    for i, row in enumerate(rows):
        label, value = row[0], row[1]
        note = row[2] if len(row) > 2 else ''
        y = i * row_h + 5
        out.append(f'<text x="0" y="{y + 11}" class="lbl">{html.escape(str(label))}</text>')
        if value is None:
            out.append(f'<text x="{label_w}" y="{y + 11}" class="lbl muted-t">{html.escape(note or "skipped")}</text>')
            continue
        w = max(1.0, bar_w * (value / value_max))
        out.append(f'<rect x="{label_w}" y="{y + 2}" width="{w:.1f}" height="{row_h - 6}" '
                   f'rx="2" fill="{color}"/>')
        txt = note if note else f'{value:,.1f}{unit}' if isinstance(value, float) else f'{value:,}{unit}'
        out.append(f'<text x="{label_w + w + 6}" y="{y + 11}" class="val">{html.escape(txt)}</text>')
    out.append('</svg>')
    return ''.join(out)


def svg_hist(values, bins=30, width=680, height=170, color='var(--accent)', log_y=False, xlabel=''):
    """Histogram of a numeric list."""
    vals = [v for v in values if v is not None]
    if not vals:
        return '<p class="muted">no data</p>'
    lo, hi = min(vals), max(vals)
    if hi <= lo:
        hi = lo + 1
    step = (hi - lo) / bins
    counts = [0] * bins
    for v in vals:
        counts[min(int((v - lo) / step), bins - 1)] += 1
    top = max(counts) or 1
    pad_l, pad_b = 42, 24
    plot_w, plot_h = width - pad_l - 8, height - pad_b - 10
    out = [_svg_open(width, height)]
    out.append(f'<line x1="{pad_l}" y1="{height - pad_b}" x2="{width - 8}" y2="{height - pad_b}" class="axis"/>')
    out.append(f'<line x1="{pad_l}" y1="10" x2="{pad_l}" y2="{height - pad_b}" class="axis"/>')
    bw = plot_w / bins
    for i, c in enumerate(counts):
        if not c:
            continue
        frac = (math.log10(c + 1) / math.log10(top + 1)) if log_y else (c / top)
        h = max(1.0, plot_h * frac)
        x = pad_l + i * bw
        out.append(f'<rect x="{x:.1f}" y="{height - pad_b - h:.1f}" width="{max(1.0, bw - 1):.1f}" '
                   f'height="{h:.1f}" fill="{color}"><title>{lo + i * step:,.0f}-{lo + (i + 1) * step:,.0f}: {c:,}</title></rect>')
    out.append(f'<text x="{pad_l}" y="{height - 6}" class="tick">{lo:,.0f}</text>')
    out.append(f'<text x="{width - 8}" y="{height - 6}" class="tick" text-anchor="end">{hi:,.0f}</text>')
    out.append(f'<text x="0" y="16" class="tick">{top:,}{" (log)" if log_y else ""}</text>')
    if xlabel:
        out.append(f'<text x="{pad_l + plot_w / 2}" y="{height - 6}" class="tick" text-anchor="middle">{html.escape(xlabel)}</text>')
    out.append('</svg>')
    return ''.join(out)


def svg_stacked(segments, width=680, height=26):
    """One stacked proportion bar: segments = [(label, value, color), ...]."""
    total = sum(s[1] for s in segments) or 1
    out = [_svg_open(width, height)]
    x = 0.0
    for label, value, color in segments:
        w = width * value / total
        if w > 0:
            out.append(f'<rect x="{x:.1f}" y="0" width="{w:.1f}" height="16" fill="{color}">'
                       f'<title>{html.escape(label)}: {value:,} ({pct(value, total):.1f}%)</title></rect>')
            if w > 46:
                out.append(f'<text x="{x + w / 2:.1f}" y="12" class="seg">{pct(value, total):.0f}%</text>')
        x += w
    out.append('</svg>')
    return ''.join(out)


# --------------------------------------------------------------------------------------
# data collection
# --------------------------------------------------------------------------------------

YP_RE = re.compile(r'(ID\d+):(\d+)-(\d+)')


# Every input the report draws on, with the pipeline stage that produces it. The report is
# written even when a run failed part-way, so it must say plainly what is missing rather than
# rendering an empty section that looks like a real zero.
INPUT_SPEC = [
    ('read_summary', "Read summary",          '{base}_read_summary.tsv',                            'read_summary'),
    ('probe',        "Y' probe table",        '{base}_post_y_prime_probe.tsv',                      'y_prime_analysis'),
    ('recomb',       "Recombination summary", 'recombination/{base}_recombination_summary.tsv',     'recombination_summary'),
    ('features',     "Per-end features",      'recombination/{base}_*_features.tsv',                'recombination_alignment'),
    ('events',       "Recombination events",  'recombination_events/{base}_all_events_summary.tsv', 'extract_recombination_events'),
    ('ypstats',      "Y' statistics",         'graphs/stats_for_y_primes/*_stats_y_prime.txt',      'single_sample_plots'),
    ('tracks',       "Track plots",           'graphs/recombination_tracks/*_tracks.png',           'recombination_track_plots'),
    ('onion',        "Onion-skin summary",    'recombination_events/{base}_onion_skin_summary.tsv', 'onion_skin'),
]


def audit_inputs(pipeline_dir, base_name):
    """Which expected inputs exist. Never raises -- a missing directory is a finding."""
    found = []
    for key, label, pattern, stage in INPUT_SPEC:
        rel = pattern.format(base=base_name)
        matches = [m for m in glob.glob(os.path.join(pipeline_dir, rel))
                   if os.path.isfile(m) and os.path.getsize(m) > 0]
        found.append({'key': key, 'label': label, 'stage': stage, 'rel': rel,
                      'n': len(matches), 'ok': bool(matches)})
    return found


def collect(pipeline_dir, base_name):
    """Gather everything the report needs from one sample's _pipeline directory."""
    d = pipeline_dir
    data = {'base_name': base_name, 'pipeline_dir': d, 'warnings': []}
    data['audit'] = audit_inputs(d, base_name)

    # -- read_summary.tsv: per-end anchored / qualifying counts, plus header comments ----
    rs_path = os.path.join(d, f'{base_name}_read_summary.tsv')
    data['total_reads'] = None
    data['filtered_reads'] = None
    if os.path.isfile(rs_path):
        comments = []
        with open(rs_path) as fh:
            for line in fh:
                if line.startswith('#'):
                    comments.append(line.strip('#').strip())
                else:
                    break
        for c in comments:
            m = re.search(r'reads after filter_reads\.py:\s*([\d,]+)', c)
            if m:
                data['filtered_reads'] = int(m.group(1).replace(',', ''))
        data['read_summary_comments'] = comments
        df = read_tsv(rs_path)
        if df is not None:
            df = df[df['chr_end'] != 'ALL'] if 'chr_end' in df.columns else df
            data['per_end_reads'] = [
                {'chr_end': r.get('chr_end', ''),
                 'anchored': to_int(r.get('anchored_reads')) or 0,
                 'qualifying': to_int(r.get('qualifying_reads')) or 0,
                 'pct_qualifying': to_num(r.get('pct_qualifying')) or 0.0}
                for _, r in df.iterrows()]
    else:
        data['warnings'].append(f'missing read_summary: {rs_path}')
        data['per_end_reads'] = []
        data['read_summary_comments'] = []

    # total reads from the FASTA index if present (cheap line count)
    fai = os.path.join(d, f'{base_name}.fasta.fai')
    if os.path.isfile(fai):
        try:
            with open(fai, 'rb') as fh:
                data['total_reads'] = sum(1 for _ in fh)
        except OSError:
            pass

    # -- probe table: per-read telomere length and Y' delta ------------------------------
    probe = read_tsv(os.path.join(d, f'{base_name}_post_y_prime_probe.tsv'))
    data['probe_rows'] = 0
    data['telo_lengths'] = []
    data['yp_delta'] = []            # qualifying reads only (matches stats_y_prime.txt)
    data['yp_delta_all'] = []        # every anchored read, for context
    data['probe_by_read'] = {}
    if probe is not None:
        data['probe_rows'] = len(probe)
        for _, r in probe.iterrows():
            rid = r.get('read_id', '')
            telo = to_num(r.get('repeat_length'))
            delta = to_num(r.get('y_primes_relative_to_ref'))
            probe_n = to_int(r.get('y_prime_probe_count'))
            ref_n = to_int(r.get('reference_y_primes'))
            # 'qualifying' == the pipeline's own definition in read_summary.py:
            # adapter found after the telomere AND repeat_length >= 30. stats_y_prime.txt
            # uses this subset, so the report must too or the two disagree.
            adapter_ok = str(r.get('Adapter_After_Telomere', '')).strip().lower() in ('true', '1')
            qualifying = adapter_ok and telo is not None and telo >= 30
            if telo is not None and telo > 0:
                data['telo_lengths'].append(telo)
            if delta is not None:
                data['yp_delta_all'].append(delta)
                if qualifying:
                    data['yp_delta'].append(delta)
            data['probe_by_read'][rid] = {
                'telo': telo, 'delta': delta, 'probe_n': probe_n, 'ref_n': ref_n,
                'chr_end': r.get('chr_end', ''),
                'adapter_after_telo': r.get('Adapter_After_Telomere', ''),
            }
    else:
        data['warnings'].append('missing *_post_y_prime_probe.tsv (Y\' and telomere panels will be empty)')

    # -- recombination summary: per-end rates -------------------------------------------
    rec = read_tsv(os.path.join(d, 'recombination', f'{base_name}_recombination_summary.tsv'))
    data['per_end_recomb'] = []
    data['recomb_reads'] = 0
    data['recomb_events'] = 0
    data['conf_weighted'] = None
    data['has_v3'], data['rc_weighted'], data['dc_weighted'] = False, None, None
    if rec is not None:
        conf_acc = rc_acc = dc_acc = 0.0
        v3_n = 0
        for _, r in rec.iterrows():
            status = r.get('status', '')
            total = to_int(r.get('total_reads')) or 0
            nrec = to_num(r.get('n_recombination'))
            row = {
                'chr_end': r.get('chr_end', ''),
                'status': status,
                'skip_reason': r.get('skip_reason', ''),
                'total_reads': total,
                'n_recombination': nrec,
                'pct': to_num(r.get('pct_recombination')),
                'n_spacer_switch': to_num(r.get('n_spacer_switch')),
                'n_x_element_switch': to_num(r.get('n_x_element_switch')),
                'n_y_prime_change': to_num(r.get('n_y_prime_change')),
                'mean_confidence': to_num(r.get('mean_confidence')),
                'mean_rc': to_num(r.get('mean_recombination_confidence')),
                'mean_dc': to_num(r.get('mean_donor_confidence')),
                'n_confident_donor': to_num(r.get('n_confident_donor')),
                'n_complex_events': to_num(r.get('n_complex_events')),
                'source': r.get('most_common_source', ''),
            }
            data['per_end_recomb'].append(row)
            if status == 'analyzed' and total:
                data['recomb_reads'] += total
                data['recomb_events'] += int(nrec or 0)
                conf_acc += (row['mean_confidence'] or 0.0) * total
                if row['mean_dc'] is not None and nrec:
                    rc_acc += (row['mean_rc'] or 0.0) * nrec
                    dc_acc += row['mean_dc'] * nrec
                    v3_n += nrec
        if data['recomb_reads']:
            data['conf_weighted'] = conf_acc / data['recomb_reads']
        # v3 scores are averaged over RECOMBINANT reads, so weight by n_recombination
        if v3_n:
            data['has_v3'] = True
            data['rc_weighted'] = rc_acc / v3_n
            data['dc_weighted'] = dc_acc / v3_n
    else:
        data['warnings'].append('missing recombination_summary.tsv (recombination panel will be empty)')

    # -- per-read features: Y' arrays, positions, confidence -----------------------------
    data['reads'] = []
    data['yp_counts'] = []
    data['gaps'] = []
    data['elem_len'] = defaultdict(list)
    data['status_counts'] = Counter()
    data['overlaps'] = 0
    for fpath in sorted(glob.glob(os.path.join(d, 'recombination', f'{base_name}_*_features.tsv'))):
        df = read_tsv(fpath)
        if df is None or 'read_id' not in df.columns:
            continue
        for _, r in df.iterrows():
            arr = [x for x in (r.get('y_prime_observed_array') or '').split(',') if x]
            ivs = [(m.group(1), int(m.group(2)), int(m.group(3)))
                   for m in YP_RE.finditer(r.get('y_prime_positions') or '')]
            ivs.sort(key=lambda t: t[1])
            for idn, a, b in ivs:
                data['elem_len'][idn].append(b - a)
            for i in range(1, len(ivs)):
                gap = ivs[i][1] - ivs[i - 1][2]
                data['gaps'].append(gap)
                if gap < 0:
                    data['overlaps'] += 1
            status = r.get('y_prime_recombination_status', '')
            if status:
                data['status_counts'][status] += 1
            rid = r.get('read_id', '')
            pb = data['probe_by_read'].get(rid, {})
            data['yp_counts'].append(len(arr))
            data['reads'].append({
                'read_id': rid,
                'chr_end': r.get('chr_end', ''),
                'read_length': to_int(r.get('read_length')),
                'telo': pb.get('telo'),
                'yp_n': len(arr),
                'yp_array': r.get('y_prime_observed_array', ''),
                'yp_delta': pb.get('delta'),
                'status': status,
                'recomb': r.get('recombination_detected', ''),
                'source': r.get('recombination_source', ''),
                'conf': to_num(r.get('overall_confidence')),
                'rc': to_num(r.get('recombination_confidence')),
                'dc': to_num(r.get('donor_confidence')),
                'complex': r.get('is_complex_event', ''),
                'compatible': r.get('y_prime_compatible_ends', ''),
            })
    if not data['reads']:
        data['warnings'].append('no *_features.tsv found (per-read table will be empty)')

    # -- onion skin: per-end summary of how gained Y' arrays were built -------------------
    onion = read_tsv(os.path.join(d, 'recombination_events', f'{base_name}_onion_skin_summary.tsv'))
    data['onion'] = [] if onion is None else onion.to_dict('records')
    return data


def classify_survivor(data):
    """Heuristic Type I / Type II signature. Deliberately reports the evidence, not a verdict."""
    deltas = data.get('yp_delta') or []
    telo = data.get('telo_lengths') or []
    if not deltas or not telo:
        return None
    positive = sum(1 for d in deltas if d > 0)
    pct_pos = pct(positive, len(deltas))
    telo_med = quantile(sorted(telo), 0.50)
    if pct_pos >= 40 and telo_med is not None and telo_med < 300:
        label, cls = 'Type I signature', 'type1'
        why = f'{pct_pos:.0f}% of read-ends gained Y&prime; while telomeres stayed short (median {telo_med:,.0f} bp)'
    elif pct_pos < 25 and telo_med is not None and telo_med >= 300:
        label, cls = 'Type II signature', 'type2'
        why = f'telomeres elongated (median {telo_med:,.0f} bp) with little Y&prime; gain ({pct_pos:.0f}% of read-ends)'
    else:
        label, cls = 'Mixed / unclear', 'mixed'
        why = f'{pct_pos:.0f}% Y&prime;-positive read-ends, telomere median {telo_med:,.0f} bp'
    return {'label': label, 'cls': cls, 'why': why, 'pct_pos': pct_pos, 'telo_med': telo_med}


# --------------------------------------------------------------------------------------
# HTML
# --------------------------------------------------------------------------------------

CSS = """
:root{--bg:#fbfbfa;--fg:#1c1b18;--muted:#6b6862;--line:#e4e2dd;--card:#fff;
--accent:#3d6b9e;--accent2:#b06a3b;--good:#2f7d52;--warn:#b0873b;--bad:#a4453a;--chip:#f0eee9;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
--bg:#1a1917;--fg:#e9e7e2;--muted:#9a968e;--line:#33312d;--card:#222120;
--accent:#7aa5d2;--accent2:#d69a6b;--good:#6bbb8c;--warn:#d9b775;--bad:#d98b80;--chip:#2b2a27;}}
:root[data-theme="dark"]{--bg:#1a1917;--fg:#e9e7e2;--muted:#9a968e;--line:#33312d;--card:#222120;
--accent:#7aa5d2;--accent2:#d69a6b;--good:#6bbb8c;--warn:#d9b775;--bad:#d98b80;--chip:#2b2a27;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;}
.wrap{max-width:1080px;margin:0 auto;padding:24px 16px 80px;}
h1{font-size:23px;margin:0 0 4px;font-weight:650;letter-spacing:-.01em}
h2{font-size:16px;margin:34px 0 10px;font-weight:600;padding-bottom:6px;border-bottom:1px solid var(--line)}
h3{font-size:13px;margin:20px 0 6px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
.sub{color:var(--muted);font-size:13px;margin:0 0 18px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:12px 0}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(142px,1fr));gap:10px;margin:14px 0}
.tile{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:11px 13px}
.tile .k{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
.tile .v{font-size:21px;font-weight:650;margin-top:3px;font-variant-numeric:tabular-nums}
.tile .s{font-size:11.5px;color:var(--muted);margin-top:2px}
table{border-collapse:collapse;width:100%;font-size:12.5px;font-variant-numeric:tabular-nums}
th,td{text-align:right;padding:5px 7px;border-bottom:1px solid var(--line);white-space:nowrap}
th:first-child,td:first-child{text-align:left}
th{font-weight:600;color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.03em;
position:sticky;top:0;background:var(--bg);cursor:pointer;user-select:none}
th:hover{color:var(--fg)}
tbody tr:hover{background:var(--chip)}
.scroll{max-height:540px;overflow:auto;border:1px solid var(--line);border-radius:8px}
.muted{color:var(--muted)} .mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11.5px}
.chip{display:inline-block;padding:2px 8px;border-radius:99px;background:var(--chip);font-size:11.5px;margin-right:5px}
.type1{background:color-mix(in srgb,var(--accent) 20%,transparent);color:var(--accent)}
.type2{background:color-mix(in srgb,var(--accent2) 22%,transparent);color:var(--accent2)}
.mixed{background:var(--chip)}
.chart .lbl{font-size:10.5px;fill:var(--muted)} .chart .val{font-size:10.5px;fill:var(--fg)}
.chart .tick{font-size:10px;fill:var(--muted)} .chart .axis{stroke:var(--line);stroke-width:1}
.chart .seg{font-size:10px;fill:#fff;text-anchor:middle;font-weight:600}
.chart .muted-t{fill:var(--muted)}
.banner{background:color-mix(in srgb,var(--warn) 18%,transparent);border:1px solid var(--warn);border-radius:8px;padding:11px 14px;margin:14px 0;font-size:13px}
.warn{border-left:3px solid var(--warn);padding-left:11px;color:var(--muted);font-size:12.5px;margin:8px 0}
.legend{font-size:11.5px;color:var(--muted);margin:6px 0 0}
.legend i{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:4px}
.note{font-size:12.5px;color:var(--muted);margin:6px 0 0;line-height:1.55}
code{background:var(--chip);padding:1px 5px;border-radius:3px;font-size:11.5px}
"""

JS = """
document.querySelectorAll('table.sortable').forEach(function(t){
  t.querySelectorAll('th').forEach(function(th,i){
    th.addEventListener('click',function(){
      var tb=t.tBodies[0], rows=Array.from(tb.rows);
      var asc=!(th.dataset.asc==='1');
      t.querySelectorAll('th').forEach(function(o){o.dataset.asc='';});
      th.dataset.asc=asc?'1':'0';
      rows.sort(function(a,b){
        var x=a.cells[i].dataset.v!==undefined?a.cells[i].dataset.v:a.cells[i].innerText;
        var y=b.cells[i].dataset.v!==undefined?b.cells[i].dataset.v:b.cells[i].innerText;
        var nx=parseFloat(x), ny=parseFloat(y);
        if(!isNaN(nx)&&!isNaN(ny)){return asc?nx-ny:ny-nx;}
        return asc?String(x).localeCompare(String(y)):String(y).localeCompare(String(x));
      });
      rows.forEach(function(r){tb.appendChild(r);});
    });
  });
});
"""


def tile(key, value, sub=''):
    s = f'<div class="s">{sub}</div>' if sub else ''
    return f'<div class="tile"><div class="k">{key}</div><div class="v">{value}</div>{s}</div>'


def fmt(value, dec=0, dash='&mdash;'):
    if value is None:
        return dash
    return f'{value:,.{dec}f}'


def sort_key_end(name):
    """chr10L -> (10, 'L') so ends order naturally rather than lexically."""
    m = re.match(r'(?:chr)?(\d+)([LR])?', str(name))
    return (int(m.group(1)), m.group(2) or '') if m else (999, str(name))


def render_sample(data):
    """Return the HTML body for one sample."""
    b = html.escape(data['base_name'])
    out = []

    # ---- headline tiles ---------------------------------------------------------------
    anchored = sum(r['anchored'] for r in data.get('per_end_reads', []))
    qualifying = sum(r['qualifying'] for r in data.get('per_end_reads', []))
    sig = classify_survivor(data)
    overall = pct(data['recomb_events'], data['recomb_reads']) if data['recomb_reads'] else None
    telo_stats = describe(data['telo_lengths'])
    deltas = data['yp_delta']
    n_pos = sum(1 for d in deltas if d > 0)
    n_neg = sum(1 for d in deltas if d < 0)
    n_neu = sum(1 for d in deltas if d == 0)

    # --- completeness banner: the first thing on the page when anything is missing ---
    audit = data.get('audit', [])
    missing = [a for a in audit if not a['ok']]
    if missing:
        out.append('<div class="banner"><b>PARTIAL REPORT</b> &mdash; '
                   f'{len(missing)} of {len(audit)} expected inputs are missing, so sections '
                   'below are incomplete. Missing values render as &mdash;, never as zero. '
                   'See <i>Pipeline completeness</i> at the foot of the page.</div>')

    out.append('<div class="tiles">')
    # NOTE: {base}.fasta is the post-filter read set, so its .fai counts filtered reads,
    # not raw basecalled reads. The raw count is not recorded in any pipeline output.
    n_filtered = data.get('filtered_reads') or data.get('total_reads')
    # A tile whose source input is missing must read as a dash, never as 0 -- a zero here
    # would be indistinguishable from a real measurement of zero.
    have = {a['key']: a['ok'] for a in audit}
    out.append(tile('Reads after filter', fmt(n_filtered), 'filter_reads.py output'))
    out.append(tile('Anchored reads', f'{anchored:,}' if have.get('read_summary') else '&mdash;',
                    (f'{pct(anchored, n_filtered or 0):.2f}% of filtered'
                     if have.get('read_summary') and n_filtered else 'read summary missing')))
    out.append(tile('Telomere-anchored', f'{qualifying:,}' if have.get('read_summary') else '&mdash;',
                    (f'{pct(qualifying, anchored):.1f}% of anchored'
                     if have.get('read_summary') and anchored else '')))
    out.append(tile('Analyzed (recomb.)',
                    f'{data["recomb_reads"]:,}' if have.get('recomb') else '&mdash;',
                    (f'{data["recomb_events"]:,} events' if have.get('recomb')
                     else 'recombination summary missing')))
    if data['has_v3']:
        conf_sub = f'donor conf {data["dc_weighted"]:.2f} &middot; call conf {data["rc_weighted"]:.2f}'
    elif data['conf_weighted'] is not None:
        conf_sub = f'conf {data["conf_weighted"]:.3f} (old score)'
    else:
        conf_sub = ''
    out.append(tile('Recombination', f'{overall:.1f}%' if overall is not None else '&mdash;', conf_sub))
    out.append(tile('Telomere median', f'{telo_stats["median"]:,.0f} bp' if telo_stats else '&mdash;',
                    f'p90 {telo_stats["p90"]:,.0f}' if telo_stats else ''))
    out.append('</div>')

    if sig:
        out.append(f'<div class="card"><span class="chip {sig["cls"]}">{sig["label"]}</span>'
                   f'<span class="muted">{sig["why"]}</span>'
                   f'<p class="note">Heuristic, from two independent measurements: Y&prime; gain per read-end '
                   f'and telomere repeat length. Type I amplifies Y&prime; and keeps short telomeres; '
                   f'Type II elongates the terminal tract and leaves Y&prime; copy number alone. '
                   f'Treat as a prompt to look, not a verdict.</p></div>')

    # ---- read accounting --------------------------------------------------------------
    out.append('<h2>Read accounting</h2>')
    if data.get('read_summary_comments'):
        out.append('<p class="note">' + '<br>'.join(html.escape(c) for c in data['read_summary_comments']) + '</p>')
    ends = sorted(data.get('per_end_reads', []), key=lambda r: sort_key_end(r['chr_end']))
    if ends:
        out.append('<h3>Anchored reads per chromosome end</h3>')
        out.append(svg_barh([(e['chr_end'], e['anchored'], f"{e['anchored']:,}") for e in ends]))
        out.append('<h3>Telomere-anchored fraction per end</h3>')
        out.append('<p class="note">Qualifying = adapter found after the telomere AND repeat length &ge; 30 bp. '
                   'A low fraction means short or absent telomere tracts, which is itself informative &mdash; '
                   'Type I survivors run low here.</p>')
        out.append(svg_barh([(e['chr_end'], e['pct_qualifying'],
                              f"{e['pct_qualifying']:.1f}%  ({e['qualifying']:,}/{e['anchored']:,})") for e in ends],
                            value_max=100, color='var(--accent2)'))

    # ---- recombination ----------------------------------------------------------------
    out.append('<h2>Recombination</h2>')
    rec_ends = sorted(data.get('per_end_recomb', []), key=lambda r: sort_key_end(r['chr_end']))
    analyzed = [r for r in rec_ends if r['status'] == 'analyzed']
    skipped = [r for r in rec_ends if r['status'] != 'analyzed']
    if analyzed:
        out.append(f'<p class="note">Overall <b>{overall:.1f}%</b> across {data["recomb_reads"]:,} analyzed reads '
                   f'({data["recomb_events"]:,} recombination events) over {len(analyzed)} of {len(rec_ends)} ends. '
                   f'This weighted figure is not in the summary TSV &mdash; that file is per-end only.</p>')
        out.append(svg_barh([(r['chr_end'], r['pct'], f"{r['pct']:.1f}%  (n={r['total_reads']:,})")
                             for r in analyzed], value_max=100))
        out.append('<h3>Per-end detail</h3>')
        v3 = data['has_v3']
        out.append('<div class="scroll"><table class="sortable"><thead><tr>'
                   '<th>end</th><th>reads</th><th>recomb</th><th>%</th><th>spacer sw</th>'
                   '<th>X sw</th><th>Y&prime; chg</th><th>complex</th>'
                   + ('<th>call conf</th><th>donor conf</th><th>confident donors</th>' if v3 else '<th>conf</th>')
                   + '<th>top source</th></tr></thead><tbody>')
        for r in analyzed:
            out.append(
                f'<tr><td>{html.escape(r["chr_end"])}</td>'
                f'<td data-v="{r["total_reads"]}">{r["total_reads"]:,}</td>'
                f'<td data-v="{r["n_recombination"] or 0}">{fmt(r["n_recombination"])}</td>'
                f'<td data-v="{r["pct"] or 0}">{fmt(r["pct"], 1)}%</td>'
                f'<td data-v="{r["n_spacer_switch"] or 0}">{fmt(r["n_spacer_switch"])}</td>'
                f'<td data-v="{r["n_x_element_switch"] or 0}">{fmt(r["n_x_element_switch"])}</td>'
                f'<td data-v="{r["n_y_prime_change"] or 0}">{fmt(r["n_y_prime_change"])}</td>'
                f'<td data-v="{r["n_complex_events"] or 0}">{fmt(r["n_complex_events"])}</td>'
                + (f'<td data-v="{r["mean_rc"] or 0}">{fmt(r["mean_rc"], 2)}</td>'
                   f'<td data-v="{r["mean_dc"] or 0}">{fmt(r["mean_dc"], 2)}</td>'
                   f'<td data-v="{r["n_confident_donor"] or 0}">{fmt(r["n_confident_donor"])}</td>' if v3 else
                   f'<td data-v="{r["mean_confidence"] or 0}">{fmt(r["mean_confidence"], 3)}</td>')
                + f'<td>{html.escape(r["source"] or "")}</td></tr>')
        out.append('</tbody></table></div>')
        if v3:
            out.append('<p class="note"><b>Two confidence scores</b>, both averaged over the recombinant '
                       'reads only. <b>call conf</b> (<code>recombination_confidence</code>): is the read '
                       'really changed? It rises with independent evidence (more gained Y&prime; copies, a '
                       'strong spacer or X-element switch) and drops for a Loss on a read that never reaches '
                       'the telomere. <b>donor conf</b> (<code>donor_confidence</code>): is the named donor '
                       'right? It rises with a long, unique Y&prime; fingerprint, agreeing evidence and a clear '
                       'vote margin, and is 0 when the donor is ambiguous. <b>confident donors</b>: reads with '
                       'donor conf &ge; 0.5. Every read lists its components in <code>confidence_basis</code>.</p>')
        else:
            out.append('<p class="note"><b>On <code>conf</code>:</b> the per-end mean averages every read. '
                   'Non-recombinant reads contribute a fixed 0.95, so the mean mostly tracks the recombination '
                   'rate: an end with few recombinants sits near 0.95 whatever the evidence. For recombinant '
                   'reads the score reflects how confidently the <i>donor end</i> was identified, not whether '
                   'recombination occurred &mdash; long tandem arrays of a common Y&prime; ID leave many '
                   'compatible donors, so it falls as copy number rises. Read it alongside the recombination '
                   'rate, not as an independent quality measure.</p>')
    if skipped:
        out.append('<h3>Skipped ends</h3><table><thead><tr><th>end</th><th>reason</th></tr></thead><tbody>')
        for r in skipped:
            out.append(f'<tr><td>{html.escape(r["chr_end"])}</td>'
                       f'<td style="text-align:left">{html.escape(r["skip_reason"] or r["status"])}</td></tr>')
        out.append('</tbody></table>')

    # ---- Y' analysis ------------------------------------------------------------------
    out.append("<h2>Y&prime; elements</h2>")
    if deltas:
        out.append('<h3>Y&prime; copy number relative to the day-0 reference</h3>')
        out.append(f'<p class="note">Over the <b>{len(deltas):,} qualifying</b> read-ends '
                   f'(adapter after telomere and repeat &ge; 30 bp), of {len(data["yp_delta_all"]):,} '
                   'anchored read-ends total. This is the same subset the pipeline\'s '
                   '<code>stats_y_prime.txt</code> reports, so the two agree.</p>')
        out.append(svg_stacked([
            ('gained Y&prime;', n_pos, 'var(--accent)'),
            ('unchanged', n_neu, 'var(--line)'),
            ('lost Y&prime;', n_neg, 'var(--bad)')]))
        out.append(f'<p class="legend"><i style="background:var(--accent)"></i>gained {n_pos:,} '
                   f'({pct(n_pos, len(deltas)):.1f}%) &nbsp; '
                   f'<i style="background:var(--line)"></i>unchanged {n_neu:,} ({pct(n_neu, len(deltas)):.1f}%) &nbsp; '
                   f'<i style="background:var(--bad)"></i>lost {n_neg:,} ({pct(n_neg, len(deltas)):.1f}%)</p>')
        pos_vals = [d for d in deltas if d > 0]
        if pos_vals:
            st = describe(pos_vals)
            out.append('<div class="tiles">')
            out.append(tile('Mean gain', f'+{st["mean"]:.2f}', 'Y&prime; copies per read-end'))
            out.append(tile('Median gain', f'+{st["median"]:.0f}'))
            out.append(tile('Max gain', f'+{st["max"]:.0f}'))
            out.append('</div>')
        out.append('<h3>Distribution of Y&prime; change per read-end</h3>')
        out.append(svg_hist(deltas, bins=min(40, int(max(deltas) - min(deltas)) + 1 or 1),
                            xlabel='Y&prime; relative to reference', log_y=True))

    if data['yp_counts']:
        out.append('<h3>Y&prime; copies observed on a single read</h3>')
        cc = Counter(data['yp_counts'])
        rows = [(str(k), cc[k], f'{cc[k]:,} reads') for k in sorted(cc)]
        out.append(svg_barh(rows, row_h=15))

    if data['gaps']:
        gs = describe(data['gaps'])
        out.append('<h3>Spacing between tandem Y&prime; copies</h3>')
        out.append('<div class="tiles">')
        out.append(tile('Median gap', f'{gs["median"]:,.0f} bp', 'internal telomeric sequence'))
        out.append(tile('10th&ndash;90th pct', f'{gs["p10"]:,.0f}&ndash;{gs["p90"]:,.0f}'))
        out.append(tile('Overlapping', f'{data["overlaps"]:,}',
                        f'{pct(data["overlaps"], len(data["gaps"])):.1f}% &mdash; should be ~0'))
        out.append('</div>')
        out.append('<p class="note">Genuine tandem Y&prime; arrays are separated by a short internal telomeric '
                   'tract, so a tight gap distribution around ~170 bp with essentially no overlaps is evidence '
                   'the array is real rather than an alignment artifact.</p>')
        out.append(svg_hist([g for g in data['gaps'] if -500 < g < 1500], bins=40, xlabel='gap (bp)'))

    if data['elem_len']:
        out.append('<h3>Observed element length by Y&prime; ID</h3>')
        out.append('<table class="sortable"><thead><tr><th>ID</th><th>n</th><th>median (bp)</th>'
                   '<th>min</th><th>max</th></tr></thead><tbody>')
        for idn in sorted(data['elem_len']):
            st = describe(data['elem_len'][idn])
            out.append(f'<tr><td>{html.escape(idn)}</td><td data-v="{st["n"]}">{st["n"]:,}</td>'
                       f'<td data-v="{st["median"]}">{st["median"]:,.0f}</td>'
                       f'<td data-v="{st["min"]}">{st["min"]:,.0f}</td>'
                       f'<td data-v="{st["max"]}">{st["max"]:,.0f}</td></tr>')
        out.append('</tbody></table>')
        out.append('<p class="note">Y&prime; length is <b>per ID</b>: Short classes run ~5.0&ndash;5.5 kb and Long '
                   '~6.5&ndash;6.9 kb. Medians far from the library value for that ID suggest mis-assignment.</p>')

    if data['status_counts']:
        out.append('<h3>Y&prime; recombination status</h3><table><thead><tr><th>status</th><th>reads</th>'
                   '<th>%</th></tr></thead><tbody>')
        tot = sum(data['status_counts'].values())
        for k, v in data['status_counts'].most_common():
            out.append(f'<tr><td>{html.escape(k)}</td><td>{v:,}</td><td>{pct(v, tot):.1f}%</td></tr>')
        out.append('</tbody></table>')

    # ---- onion skin -------------------------------------------------------------------
    onion = [r for r in data.get('onion', []) if r.get('chr_end') == 'ALL' or (to_int(r.get('n_gain_like')) or 0) > 0]
    if onion:
        out.append('<h3>Onion skin: how gained Y&prime; arrays were built</h3>')
        out.append('<p class="note">Each gained array is split into donor pieces by the path parser. A '
                   '<b>same-donor repeat</b> is a circle (a donor piece copied in tandem) or one donor '
                   'giving two or more pieces. <b>Circles</b> are graded strong / moderate / weak by '
                   'support; <b>unassigned</b> circles are tandem copies whose donor cannot be named '
                   '(a Y&prime; ID found at several ends). Click an end for its per-read schematic.</p>')
        out.append('<div class="scroll"><table class="sortable"><thead><tr><th>end</th><th>gain-like</th>'
                   '<th>&ge;2 Y&prime;</th><th>1 donor</th><th>multi-donor</th><th>same-donor</th>'
                   '<th>circles</th><th>S/M/W</th><th>unassigned</th><th>top donors</th></tr></thead><tbody>')
        for r in onion:
            end = str(r.get('chr_end', ''))
            png = f'graphs/onion_skin/{data["base_name"]}_{end}_ycopies_schematic.png'
            cell = (html.escape(end) if end == 'ALL' else
                    f'<a href="{html.escape(png)}">{html.escape(end)}</a>')
            n = lambda k: to_int(r.get(k)) or 0
            out.append(
                f'<tr><td>{cell}</td><td data-v="{n("n_gain_like")}">{n("n_gain_like"):,}</td>'
                f'<td data-v="{n("n_gain_2plus")}">{n("n_gain_2plus"):,}</td>'
                f'<td data-v="{n("n_single_donor")}">{n("n_single_donor"):,}</td>'
                f'<td data-v="{n("n_multi_donor")}">{n("n_multi_donor"):,}</td>'
                f'<td data-v="{n("n_same_donor_repeat")}">{n("n_same_donor_repeat"):,}</td>'
                f'<td data-v="{n("n_circle")}">{n("n_circle"):,}</td>'
                f'<td>{n("n_circle_strong")}/{n("n_circle_moderate")}/{n("n_circle_weak")}</td>'
                f'<td data-v="{n("n_circle_unassigned")}">{n("n_circle_unassigned"):,}</td>'
                f'<td style="text-align:left">{html.escape(str(r.get("top_donors") or ""))}</td></tr>')
        out.append('</tbody></table></div>')

    # ---- telomere ---------------------------------------------------------------------
    out.append('<h2>Telomere repeat length</h2>')
    if telo_stats:
        out.append('<div class="tiles">')
        for k, lab in (('n', 'Reads w/ tract'), ('median', 'Median'), ('mean', 'Mean'),
                       ('p90', '90th pct'), ('p99', '99th pct'), ('max', 'Max')):
            out.append(tile(lab, f'{telo_stats[k]:,.0f}' + ('' if k == 'n' else ' bp')))
        out.append('</div>')
        out.append(svg_hist(data['telo_lengths'], bins=45, log_y=True,
                            xlabel='telomere repeat length (bp)', color='var(--accent2)'))
        out.append(f'<p class="note">{telo_stats["n"]:,} of {data["probe_rows"]:,} probe-table reads carry a '
                   'detectable tract. A low fraction is expected for Type I survivors, whose telomeres are '
                   'often too short to measure.</p>')
    else:
        out.append('<p class="muted">No telomere length data available.</p>')

    # ---- per-read table ---------------------------------------------------------------
    out.append('<h2>Per-read detail</h2>')
    if data['reads']:
        rows = sorted(data['reads'], key=lambda r: (-(r['yp_n'] or 0), r['chr_end']))
        v3r = any(r.get('rc') is not None for r in rows)
        out.append(f'<p class="note">{len(rows):,} reads, sorted by Y&prime; copy number. Click any header to '
                   're-sort. <code>compatible ends</code> shows the donor ambiguity behind a low confidence.</p>')
        out.append('<div class="scroll"><table class="sortable"><thead><tr>'
                   '<th>read</th><th>end</th><th>len</th><th>telo</th><th>Y&prime;</th><th>&Delta;Y&prime;</th>'
                   + ('<th>call</th><th>donor</th>' if v3r else '<th>conf</th>')
                   + '<th>status</th><th>source</th><th>array</th><th>compatible ends</th>'
                   '</tr></thead><tbody>')
        for r in rows:
            out.append(
                f'<tr><td class="mono">{html.escape((r["read_id"] or "")[:8])}</td>'
                f'<td>{html.escape(r["chr_end"] or "")}</td>'
                f'<td data-v="{r["read_length"] or 0}">{fmt(r["read_length"])}</td>'
                f'<td data-v="{r["telo"] or 0}">{fmt(r["telo"])}</td>'
                f'<td data-v="{r["yp_n"]}">{r["yp_n"]}</td>'
                f'<td data-v="{r["yp_delta"] if r["yp_delta"] is not None else 0}">'
                f'{("+" if (r["yp_delta"] or 0) > 0 else "") + fmt(r["yp_delta"])}</td>'
                + (f'<td data-v="{r["rc"] or 0}">{fmt(r["rc"], 2)}</td>'
                   f'<td data-v="{r["dc"] if r["dc"] is not None else -1}">{fmt(r["dc"], 2)}</td>' if v3r else
                   f'<td data-v="{r["conf"] or 0}">{fmt(r["conf"], 3)}</td>')
                + f'<td style="text-align:left">{html.escape(r["status"] or "")}</td>'
                f'<td>{html.escape(r["source"] or "")}</td>'
                f'<td class="mono" style="text-align:left">{html.escape((r["yp_array"] or "")[:60])}</td>'
                f'<td class="mono" style="text-align:left">{html.escape((r["compatible"] or "")[:50])}</td></tr>')
        out.append('</tbody></table></div>')
    else:
        out.append('<p class="muted">No per-read feature data available.</p>')

    # ---- pipeline completeness --------------------------------------------------------
    out.append('<h2>Pipeline completeness</h2>')
    out.append('<p class="note">What this report could and could not read. A report is written '
               'even when a run fails part-way, so this table is how you tell an empty section '
               'from a genuine zero.</p>')
    out.append('<table><thead><tr><th>input</th><th>produced by</th><th>files</th>'
               '<th>status</th></tr></thead><tbody>')
    for a in data.get('audit', []):
        mark = ('<span style="color:var(--good)">present</span>' if a['ok']
                else '<span style="color:var(--bad)">MISSING</span>')
        out.append(f'<tr><td>{html.escape(a["label"])}</td>'
                   f'<td style="text-align:left"><code>{html.escape(a["stage"])}</code></td>'
                   f'<td>{a["n"] or "&mdash;"}</td><td>{mark}</td></tr>')
    out.append('</tbody></table>')
    for w in data['warnings']:
        out.append(f'<div class="warn">{html.escape(w)}</div>')

    return ''.join(out)


def page(title, body, subtitle=''):
    return (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{html.escape(title)}</title><style>{CSS}</style></head><body><div class="wrap">'
            f'<h1>{html.escape(title)}</h1><p class="sub">{subtitle}</p>{body}'
            f'</div><script>{JS}</script></body></html>')


def render_comparison(all_data):
    """Compact side-by-side panel across samples."""
    out = ['<h2>Panel overview</h2>']
    out.append('<div class="scroll"><table class="sortable"><thead><tr>'
               '<th>sample</th><th>anchored</th><th>analyzed</th><th>recomb %</th><th>donor conf</th>'
               '<th>Y&prime;+ %</th><th>mean gain</th><th>max gain</th><th>telo median</th>'
               '<th>telo p90</th><th>signature</th></tr></thead><tbody>')
    for d in all_data:
        anchored = sum(r['anchored'] for r in d.get('per_end_reads', []))
        overall = pct(d['recomb_events'], d['recomb_reads']) if d['recomb_reads'] else None
        deltas = d['yp_delta']
        pos = [x for x in deltas if x > 0]
        pct_pos = pct(len(pos), len(deltas)) if deltas else None
        gain = describe(pos)
        telo = describe(d['telo_lengths'])
        sig = classify_survivor(d)
        conf_cell = (f'<td data-v="{d["dc_weighted"] or 0}">{fmt(d["dc_weighted"], 2)}</td>' if d.get('has_v3') else
                     f'<td data-v="{d["conf_weighted"] or 0}">{fmt(d["conf_weighted"], 3)} (old)</td>')
        out.append(
            f'<tr><td>{html.escape(d["base_name"])}</td>'
            f'<td data-v="{anchored}">{anchored:,}</td>'
            f'<td data-v="{d["recomb_reads"]}">{d["recomb_reads"]:,}</td>'
            f'<td data-v="{overall or 0}">{fmt(overall, 1)}%</td>'
            + conf_cell +
            f'<td data-v="{pct_pos or 0}">{fmt(pct_pos, 1)}%</td>'
            f'<td data-v="{gain.get("mean") or 0}">{"+" + fmt(gain.get("mean"), 2) if gain else "&mdash;"}</td>'
            f'<td data-v="{gain.get("max") or 0}">{"+" + fmt(gain.get("max")) if gain else "&mdash;"}</td>'
            f'<td data-v="{telo.get("median") or 0}">{fmt(telo.get("median"))}</td>'
            f'<td data-v="{telo.get("p90") or 0}">{fmt(telo.get("p90"))}</td>'
            f'<td><span class="chip {sig["cls"] if sig else "mixed"}">'
            f'{sig["label"] if sig else "&mdash;"}</span></td></tr>')
    out.append('</tbody></table></div>')

    labels = [d['base_name'] for d in all_data]
    out.append('<h3>Recombination rate</h3>')
    out.append(svg_barh([(l[:26], pct(d['recomb_events'], d['recomb_reads']) if d['recomb_reads'] else None,
                          f"{pct(d['recomb_events'], d['recomb_reads']):.1f}%" if d['recomb_reads'] else '')
                         for l, d in zip(labels, all_data)], value_max=100, label_w=190))
    out.append('<h3>Y&prime;-positive read-ends</h3>')
    out.append(svg_barh([(l[:26], pct(sum(1 for x in d['yp_delta'] if x > 0), len(d['yp_delta']))
                          if d['yp_delta'] else None,
                          f"{pct(sum(1 for x in d['yp_delta'] if x > 0), len(d['yp_delta'])):.1f}%"
                          if d['yp_delta'] else '')
                         for l, d in zip(labels, all_data)], value_max=100, label_w=190))
    out.append('<h3>Median telomere repeat length</h3>')
    out.append(svg_barh([(l[:26], describe(d['telo_lengths']).get('median'),
                          f"{describe(d['telo_lengths']).get('median', 0):,.0f} bp" if d['telo_lengths'] else '')
                         for l, d in zip(labels, all_data)], label_w=190, color='var(--accent2)'))
    out.append('<p class="note">Type I and Type II separate on these last two charts: Type I amplifies Y&prime; '
               'while telomeres stay short; Type II is the mirror image.</p>')
    return ''.join(out)


# --------------------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description='Build a self-contained HTML report for a sample.')
    p.add_argument('--pipeline-dir', help='results/{base}/_pipeline for a single sample')
    p.add_argument('--base-name', help='sample base name (single-sample mode)')
    p.add_argument('--output', help='output HTML path (single-sample mode)')
    p.add_argument('--results-dir', help='results/ directory (multi-sample mode)')
    p.add_argument('--samples', nargs='+', help='sample base names (multi-sample mode)')
    p.add_argument('--compare-output', help='write a panel comparison HTML here')
    return p.parse_args()


def build_one(pipeline_dir, base_name, output):
    print(f'  collecting: {base_name}')
    data = collect(pipeline_dir, base_name)
    body = render_sample(data)
    sub = (f'day-0 comparison &middot; generated {datetime.now():%Y-%m-%d %H:%M} &middot; '
           f'<span class="mono">{html.escape(pipeline_dir)}</span>')
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, 'w') as fh:
        fh.write(page(f'{base_name} — TeloTracker report', body, sub))
    print(f'  wrote: {output}  ({os.path.getsize(output) / 1024:.0f} KB)')
    return data


def main():
    args = parse_args()

    if args.samples:
        if not args.results_dir:
            print('ERROR: --samples requires --results-dir')
            sys.exit(1)
        collected = []
        for s in args.samples:
            pdir = os.path.join(args.results_dir, s, '_pipeline')
            if not os.path.isdir(pdir):
                print(f'  SKIP {s}: no {pdir}')
                continue
            collected.append(build_one(pdir, s, os.path.join(pdir, f'{s}_report.html')))
        if args.compare_output and collected:
            body = render_comparison(collected)
            os.makedirs(os.path.dirname(os.path.abspath(args.compare_output)), exist_ok=True)
            with open(args.compare_output, 'w') as fh:
                fh.write(page('TeloTracker panel comparison', body,
                              f'{len(collected)} samples &middot; generated {datetime.now():%Y-%m-%d %H:%M}'))
            print(f'  wrote: {args.compare_output}')
        return

    if not (args.pipeline_dir and args.base_name):
        print('ERROR: need --pipeline-dir and --base-name (or --results-dir with --samples)')
        sys.exit(1)
    out = args.output or os.path.join(args.pipeline_dir, f'{args.base_name}_report.html')
    build_one(args.pipeline_dir, args.base_name, out)


if __name__ == '__main__':
    main()
