#!/usr/bin/env python3
"""Y'-weighted, state-based recombination scoring (prototype).

Replaces the current reconcile_features() model, which treats the three axes as
independent equal voters, takes base_confidence = max(spacer_conf, x_conf), and
falls back to a hardcoded 0.3 whenever only Y' fires -- i.e. it treats Y'-only
as the WEAKEST case. Measurements on 7172 day-6 argue the opposite:

  * 542 reads (14%) are detected by Y' ALONE -- the single largest event class
  * among reads where both Y' and an upstream axis can name a donor, they
    DISAGREE 57.4% of the time (63 agree vs 85 conflict)
  * 207 reads have an upstream switch while Y' actively reports "No Change"

Y' carries more information (variable copy number, genuinely divergent
sequences) than X-element clusters, which are only ~76-88% distinct between
ends. So Y' is weighted highest, and evidence that cannot be checked against Y'
is penalised.

Crucially, Y' has three distinct states, not two:
    fires    -- Y' detects a change (strongest evidence)
    silent   -- Y' array is gone (Y' Loss); it CANNOT speak. Penalise, but this
                is not evidence against recombination.
    native   -- Y' asserts "No Change". If an upstream axis claims a switch,
                this is evidence AGAINST it, and is treated more harshly.
The current code cannot distinguish silent from native.
"""
import argparse, collections, csv, glob, re

# --- weights: Y' highest (inverts the current spacer 0.9 / y 0.6 / x 0.4) ---
Y_BASE   = 0.85   # Y' is the most trusted axis, so it can score high alone
W_YPRIME = 0.60
W_SPACER = 0.25
W_XELEM  = 0.15

PENALTY_SILENT     = 0.50   # bucket 0: Y' gone, cannot corroborate
PENALTY_CONTRADICT = 0.35   # Y' says native while upstream claims a switch
CONFLICT_CAP       = 0.70   # both fire but name different ends
AGREE_BOOST        = 1.15


def noisy_or(vals):
    p = 1.0
    for v in vals: p *= (1.0 - max(0.0, min(1.0, v)))
    return 1.0 - p


def y_detect_conf(status, n_obs, downstream_consistent):
    if status in ('No Change', '', 'no_data'): return 0.0
    clean = 1.0 if downstream_consistent else 0.7
    support = min(1.0, n_obs / 2.0) if n_obs else 0.5
    return Y_BASE * clean * support


def classify(r, yprime_null_ends):
    """-> (state, upstream_donors{axis:end}, y_state, y_candidates, confs)"""
    ce = r['chr_end']
    sw = lambda v: v.strip().lower() not in ('no_change', 'no_data', '')
    def f(k):
        try: return float(r.get(k) or 0)
        except ValueError: return 0.0

    up, up_confs = {}, []
    if sw(r['spacer_recombination']) and r['spacer_source'] and r['spacer_source'] != ce:
        up['spacer'] = r['spacer_source']; up_confs.append(f('spacer_confidence'))
    if sw(r['x_element_recombination']) and r['x_element_source'] and r['x_element_source'] != ce:
        up['x_element'] = r['x_element_source']; up_confs.append(f('x_element_confidence'))

    status = r['y_prime_recombination_status']
    obs = [x for x in (r['y_prime_observed_array'] or '').split(',') if x]
    try: d = int(r['y_prime_divergence_idx'])
    except ValueError: d = -1
    gained = obs[d:] if d >= 0 else []
    ycand = [e for e in (r['y_prime_compatible_ends'] or '').split(',') if e and e != ce]
    dsc = str(r.get('y_prime_downstream_consistent', 'True')).lower() == 'true'
    yc = y_detect_conf(status, len(obs), dsc)

    ref_has_yp = r['chr_end'] not in yprime_null_ends
    if status == 'No Change' and not ref_has_yp: y_state = 'uninformative'
    elif status == 'No Change':          y_state = 'native'
    elif not gained:                     y_state = 'silent'
    else:                                y_state = 'fires'

    if not up and y_state == 'native':                     state = 'none'
    elif not up and y_state in ('uninformative',):         state = 'none'
    elif not up and y_state != 'native':                   state = 'yprime_only'
    elif up and y_state == 'uninformative':                state = 'yprime_uninformative'
    elif up and y_state == 'native':                       state = 'contradicts_yprime'
    elif up and y_state == 'silent':                       state = 'unverified_by_yprime'
    elif up and ycand and any(v in ycand for v in up.values()): state = 'agree'
    elif up and y_state == 'fires':                        state = 'donor_conflict'
    else:                                                  state = 'other'
    return state, up, y_state, ycand, yc, (max(up_confs) if up_confs else 0.0)


def score(state, up, ycand, yc, upc):
    """-> (detection_confidence, donor, donor_confidence, flag)"""
    if state == 'none':
        return 0.0, '', 0.0, ''
    if state == 'agree':
        donor = next(v for v in up.values() if v in ycand)
        return (min(0.98, noisy_or([yc, upc]) * AGREE_BOOST), donor,
                min(0.95, noisy_or([yc, upc]) * AGREE_BOOST), '')
    if state == 'yprime_only':
        # Y' trusted on its own; donor confidence falls with candidate count
        dc = (yc / len(ycand)) if ycand else 0.0
        return yc, (ycand[0] if len(ycand) == 1 else 'ambiguous'), round(dc, 3), ''
    if state == 'donor_conflict':
        return (min(CONFLICT_CAP, noisy_or([yc, upc])), 'conflict', 0.10,
                'donor_conflict(possible_multi_jump)')
    if state == 'yprime_uninformative':          # Y'-null end: Y' simply cannot help
        return (upc * 0.85, list(up.values())[0], upc * 0.85, 'yprime_uninformative')
    if state == 'unverified_by_yprime':          # bucket 0
        return (upc * PENALTY_SILENT, list(up.values())[0], upc * PENALTY_SILENT,
                'unverified_by_yprime')
    if state == 'contradicts_yprime':
        return (upc * PENALTY_CONTRADICT, list(up.values())[0], upc * PENALTY_CONTRADICT,
                'contradicts_yprime')
    return 0.0, '', 0.0, 'unclassified'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features-dir', required=True)
    ap.add_argument('--lib', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    have = set()
    for line in open(a.lib):
        if line.startswith('>'):
            for t in re.findall(r'chr\d+[LR]', line.split('#')[0]): have.add(t)
    YNULL = {f'chr{n}{s2}' for n in range(1,17) for s2 in 'LR'} - have
    rows = []
    for fp in sorted(glob.glob(f'{a.features_dir}/*_features.tsv')):
        rows.extend(csv.DictReader(open(fp), delimiter='\t'))

    out = []
    agg = collections.defaultdict(lambda: {'n':0,'old':0.0,'det':0.0,'don':0.0})
    for r in rows:
        state, up, y_state, ycand, yc, upc = classify(r, YNULL)
        det, donor, dconf, flag = score(state, up, ycand, yc, upc)
        try: old = float(r.get('overall_confidence') or 0)
        except ValueError: old = 0.0
        s = agg[state]; s['n']+=1; s['old']+=old; s['det']+=det; s['don']+=dconf
        out.append((r['read_id'], r['chr_end'], state, y_state,
                    r['y_prime_recombination_status'],
                    ';'.join(f'{k}={v}' for k,v in up.items()) or '-',
                    r['recombination_source'] or '-', round(old,3),
                    round(det,3), donor or '-', round(dconf,3), flag))

    print(f"{'state':<22} {'n':>5} {'%':>6} | {'OLD conf':>9} | {'NEW detect':>11} {'NEW donor':>10}")
    print('-'*76)
    tot = sum(v['n'] for v in agg.values())
    order = ['agree','yprime_only','yprime_uninformative','donor_conflict',
             'unverified_by_yprime','contradicts_yprime','other','none']
    for k in order:
        if k not in agg: continue
        v = agg[k]; n = v['n']
        print(f"{k:<22} {n:5d} {100*n/tot:5.1f}% | {v['old']/n:9.3f} | "
              f"{v['det']/n:11.3f} {v['don']/n:10.3f}")
    rec = [k for k in order if k not in ('none','other')]
    nrec = sum(agg[k]['n'] for k in rec if k in agg)
    print(f"\nrecombinant reads: {nrec}")
    print('flagged:')
    for k in ('donor_conflict','unverified_by_yprime','contradicts_yprime','yprime_uninformative'):
        if k in agg: print(f"   {agg[k]['n']:5d}  {k}")

    with open(a.out,'w',newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['read_id','chr_end','state','y_state','y_status','upstream',
                    'old_source','old_confidence','detection_confidence','donor',
                    'donor_confidence','flag'])
        w.writerows(out)
    print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()
