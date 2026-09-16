#!/usr/bin/env python3
"""Directional tiled-window recombination scan, per the user's described method:

Non-overlapping 300bp windows starting at the Y' anchor-side boundary. Each window is
compared, at the SAME absolute offset, to (a) the read's own/recipient reference and
(b) every other library element (excluding same-cut99-group members). While windows keep
scoring >=99% identity against own, keep going. Once 3 CONSECUTIVE windows score <99%
against own, check whether those same 3 windows all score higher against one single other
element at that same offset -- if so, that is the breakpoint (end of the 3rd window) and
that element is the candidate donor. If not, keep extending the streak and re-check with the
latest 3 windows; a window scoring >=99% against own resets the streak.

Once a breakpoint + donor are found: global-align the read's first part against a
same-length slice from the START of the recipient reference, and the read's second part
against a same-length slice from the END of the donor reference. Report both identities;
"high" (pass) is judged at a >=90% threshold, but exact numbers are always reported.
"""
import sys, csv, json, edlib

WINDOW = 300
OWN_THRESHOLD = 99.0
PASS_THRESHOLD = 90.0

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(COMP)[::-1]

def load_fasta(f):
    d = {}; n = None; b = []
    for l in open(f):
        if l.startswith('>'):
            if n: d[n] = ''.join(b)
            n = l[1:].split()[0]; b = []
        else: b.append(l.strip())
    if n: d[n] = ''.join(b)
    return d

def nw_identity(a, b):
    if not a or not b: return 0.0
    res = edlib.align(a, b, mode="NW", task="distance")
    return 100 * (1 - res['editDistance'] / max(len(a), len(b)))

def best_strand_sub(sub, own_ref):
    best = None
    for r in (sub, rc(sub)):
        res = edlib.align(own_ref, r, mode="HW", task="distance")
        ident = 100 * (1 - res['editDistance'] / len(own_ref))
        if best is None or ident > best[0]:
            best = (ident, r)
    return best[1]

def window_id(read_win, ref, offset):
    """identity of read_win vs ref[offset:offset+len(read_win)], or None if ref too short"""
    ref_win = ref[offset:offset + len(read_win)]
    if len(ref_win) < len(read_win) * 0.8:
        return None
    return nw_identity(read_win, ref_win)

def scan_read(oriented, own_ref, candidates):
    """candidates: {name: seq}, already excludes own element + own group.
    Returns dict with breakpoint info, or None if no breakpoint found."""
    n_windows = len(oriented) // WINDOW
    own_streak = []  # list of window indices currently failing vs own
    for i in range(n_windows):
        off = i * WINDOW
        read_win = oriented[off:off + WINDOW]
        own_id = window_id(read_win, own_ref, off)
        if own_id is None:
            break  # ran out of own reference length
        if own_id >= OWN_THRESHOLD:
            own_streak = []
            continue
        own_streak.append(i)
        if len(own_streak) < 3:
            continue
        # check the latest 3 failing windows against every candidate
        last3 = own_streak[-3:]
        for name, dref in candidates.items():
            ok = True
            per_window = []
            for w in last3:
                woff = w * WINDOW
                rw = oriented[woff:woff + WINDOW]
                d_id = window_id(rw, dref, woff)
                o_id = window_id(rw, own_ref, woff)
                if d_id is None or o_id is None or not (d_id > o_id):
                    ok = False
                    break
                per_window.append((w, o_id, d_id))
            if ok:
                breakpoint = (last3[-1] + 1) * WINDOW
                return {'donor': name, 'breakpoint': breakpoint,
                        'trigger_windows': per_window}
    return None

def run(tag, reads_path, elems, groups, mismatches_rows):
    reads = load_fasta(reads_path)
    results = []
    for row in mismatches_rows:
        rid = row['read_id']
        own_elem = row['own_elem'] if 'own_elem' in row else row['true_element'].replace('E-', '')
        seq = reads.get(rid)
        if seq is None:
            continue
        ys, ye = int(row['yp_start']), int(row['yp_end'])
        core = seq[ys:ye]
        oriented = best_strand_sub(core, elems[own_elem])
        own_group = groups.get('E-' + own_elem)
        candidates = {name: s for name, s in elems.items()
                      if name != own_elem and groups.get('E-' + name) != own_group}
        found = scan_read(oriented, elems[own_elem], candidates)
        if found is None:
            results.append({'strain': tag, 'read_id': rid, 'own_elem': own_elem,
                             'breakpoint_found': False})
            continue
        bp = found['breakpoint']
        donor = found['donor']
        first_read = oriented[:bp]
        second_read = oriented[bp:]
        first_own_ref = elems[own_elem][:len(first_read)]
        second_donor_ref = elems[donor][-len(second_read):] if len(second_read) <= len(elems[donor]) else elems[donor]
        id_first = nw_identity(first_read, first_own_ref)
        id_second = nw_identity(second_read, second_donor_ref)
        passed = id_first >= PASS_THRESHOLD and id_second >= PASS_THRESHOLD
        results.append({'strain': tag, 'read_id': rid, 'own_elem': own_elem,
                         'breakpoint_found': True, 'breakpoint_bp': bp, 'donor': donor,
                         'trigger_windows': found['trigger_windows'],
                         'first_half_len': len(first_read), 'second_half_len': len(second_read),
                         'id_vs_recipient_first_half': round(id_first, 2),
                         'id_vs_donor_second_half': round(id_second, 2),
                         'pass': passed})
    return results

all_results = []

# 6991
elems58 = load_fasta('all_elements_58.fasta')
groups6991 = json.load(open('/home/andrey/teloTracker_snakemake/verification/reports/cut99_6991_day0/groups.json'))['groups']
master6991 = list(csv.DictReader(open('all58/master_results.tsv'), delimiter='\t'))
coords6991 = {r['read_id']: r for r in csv.DictReader(open('all58/all_mismatches.tsv'), delimiter='\t')}
mid6991 = [r for r in master6991 if "partial junction" in r['verdict']]
for r in mid6991:
    r['yp_start'] = coords6991[r['read_id']]['yp_start']
    r['yp_end'] = coords6991[r['read_id']]['yp_end']
all_results += run('6991', 'all_mismatch_reads.fasta', elems58, groups6991, mid6991)

# 7172 / 7302
for tag in ('7172', '7302'):
    elems = json.load(open(f'73xx/elems_{tag}.json'))
    groups = json.load(open(f'73xx/groups_{tag}.json'))['groups']
    master = list(csv.DictReader(open(f'73xx/master_{tag}.tsv'), delimiter='\t'))
    mid = [r for r in master if "partial junction" in r['verdict']]
    all_results += run(tag, f'73xx/reads_{tag}.fasta', elems, groups, mid)

print(f'{"strain":<7}{"read_id":<22}{"bp found":<10}{"donor":<12}{"bp":>6}  {"1st vs own":>11}{"2nd vs donor":>13}  pass?')
for r in all_results:
    if not r['breakpoint_found']:
        print(f"{r['strain']:<7}{r['read_id']:<22}{'NO':<10}")
        continue
    print(f"{r['strain']:<7}{r['read_id']:<22}{'yes':<10}{r['donor']:<12}{r['breakpoint_bp']:>6}  {r['id_vs_recipient_first_half']:>10.2f}%{r['id_vs_donor_second_half']:>12.2f}%  {'PASS' if r['pass'] else 'fail'}")

n_bp = sum(1 for r in all_results if r['breakpoint_found'])
n_pass = sum(1 for r in all_results if r.get('pass'))
print(f'\n{len(all_results)} reads tested; {n_bp} found a qualifying breakpoint; {n_pass} passed the final high-identity check')

with open('tiled_scan_results.tsv', 'w', newline='') as out:
    fields = ['strain','read_id','own_elem','breakpoint_found','breakpoint_bp','donor',
              'first_half_len','second_half_len','id_vs_recipient_first_half','id_vs_donor_second_half','pass']
    w = csv.DictWriter(out, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    for r in all_results: w.writerow(r)
