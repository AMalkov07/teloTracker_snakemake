#!/usr/bin/env python3
"""Classify how many recombination jumps are needed to explain a read's gained Y' array.

Adds a `jump_model` column:

    no_gained_yprimes  the read lost Y' or has none gained -- nothing to explain
    1_jump             the gained array appears intact in ONE end's native array;
                       that end is reported in `jump_donor`
    2_jumps            the gained array cannot come from one end, but splits into
                       two contiguous blocks that each exist at some end
    >2_jumps           not explainable by two blocks
    no_match           at least one gained element matches no end at all

DELIBERATELY NO DONOR PREDICTION FOR 2_jumps. Enumerating all valid 2-jump
solutions on 7172 day-6 gave a mean of 19.4 distinct donor pairs per read
(max 40) and ZERO reads with a unique solution -- seven different pairs each
explained the same 26 reads. With only 7 Y' cluster IDs across 18 Y'-carrying
ends, splitting a short array into two blocks that each exist somewhere is
nearly always possible, so any single pair reported would be an artifact of
iteration order rather than evidence.

Direct global alignment of individual reads confirms why: a gained element's
cluster is identified confidently (94.7-99.7% identity, always the assigned
cluster), but the margin between candidate ENDS within that cluster is only
0.02-0.99% -- matching the reference's own nearest-competitor divergences
(chr12R/chr4R 0.03%, chr9L/chr10L 0.23%, chr2L/chr6L 0.99%). The ceiling is
structural: donor ends carry near-identical Y' sequences, so no read quality or
matching method can resolve them.

So the 1-jump donor is trustworthy; beyond that, only the jump COUNT is.
"""
import argparse, collections, csv, glob, re


def read_fasta(p):
    n, buf = None, []
    for line in open(p):
        if line.startswith('>'):
            if n: yield n, ''.join(buf)
            n, buf = line[1:].split()[0], []
        else: buf.append(line.strip())
    if n: yield n, ''.join(buf)


def build_native(lib):
    """chr_end -> ordered list of Y' cluster IDs (anchor->telomere)."""
    pos = collections.defaultdict(dict)
    for line in open(lib):
        if not line.startswith('>'): continue
        h = line[1:].strip()
        nm, c = (h.split('#', 1) + [''])[:2]
        q = c.split('/')
        yid = q[2].split('_', 1)[0] if len(q) >= 3 else ''
        if not yid: continue
        for tok in re.sub(r'^Y_Prime_', '', nm.split()[0]).split(';'):
            m = re.match(r'(chr\d+[LR])([\d,]+)$', tok)
            if m:
                for x in m.group(2).split(','):
                    if x: pos[m.group(1)][int(x)] = yid
    return {e: [d[k] for k in sorted(d)] for e, d in pos.items()}


def ends_containing(block, native, exclude):
    """Ends whose native array contains `block` as a contiguous run."""
    out, B = [], len(block)
    for end, nat in native.items():
        if end == exclude or len(nat) < B: continue
        if any(all(nat[s + k] == block[k] for k in range(B))
               for s in range(len(nat) - B + 1)):
            out.append(end)
    return sorted(out)


def classify_jumps(gained, native, exclude):
    """-> (jump_model, donor_for_1_jump_or_empty)"""
    if not gained:
        return 'no_gained_yprimes', ''

    # every element must exist somewhere, else the array is unexplainable
    for el in set(gained):
        if not ends_containing([el], native, exclude):
            return 'no_match', ''

    one = ends_containing(gained, native, exclude)
    if one:
        return '1_jump', ','.join(one)

    # only look for a 2-block split once 1 jump has failed
    for k in range(1, len(gained)):
        if ends_containing(gained[:k], native, exclude) and \
           ends_containing(gained[k:], native, exclude):
            return '2_jumps', ''          # donor deliberately not reported

    return '>2_jumps', ''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features-dir', required=True)
    ap.add_argument('--lib', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    native = build_native(a.lib)
    rows = []
    for fp in sorted(glob.glob(f'{a.features_dir}/*_features.tsv')):
        rows.extend(csv.DictReader(open(fp), delimiter='\t'))

    out, tally = [], collections.Counter()
    for r in rows:
        obs = [x for x in (r['y_prime_observed_array'] or '').split(',') if x]
        try: d = int(r['y_prime_divergence_idx'])
        except ValueError: d = -1
        gained = obs[d:] if (d >= 0 and obs) else []
        model, donor = classify_jumps(gained, native, r['chr_end'])
        tally[model] += 1
        out.append((r['read_id'], r['chr_end'], r['y_prime_recombination_status'],
                    ','.join(gained), model, donor))

    print(f'{"jump_model":<20} {"reads":>7} {"%":>7}')
    print('-' * 36)
    tot = sum(tally.values())
    for k in ('1_jump', '2_jumps', '>2_jumps', 'no_match', 'no_gained_yprimes'):
        if k in tally:
            print(f'{k:<20} {tally[k]:7d} {100*tally[k]/tot:6.1f}%')

    rec = [o for o in out if o[4] not in ('no_gained_yprimes',)]
    if rec:
        print(f'\nof {len(rec)} reads with a gained array:')
        for k in ('1_jump', '2_jumps', '>2_jumps', 'no_match'):
            if k in tally:
                print(f'   {tally[k]:5d} ({100*tally[k]/len(rec):5.1f}%)  {k}')

    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['read_id', 'chr_end', 'y_prime_recombination_status',
                    'gained_array', 'jump_model', 'jump_donor'])
        w.writerows(out)
    print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()
