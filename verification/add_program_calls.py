#!/usr/bin/env python3
"""Add the pipeline's OWN per-read call to the cut99 annotated summaries.

The summaries so far describe what the evidence shows. This adds what the program actually
reported for the same read, so the two can be compared directly:

  program_assigned_yprime  y_prime_observed_array   -- the Y' element(s) the program assigned
  program_status           y_prime_recombination_status
  program_source           recombination_source     -- the end it attributed the Y' to
  program_mechanism        recombination_mechanism
  program_path             y_prime_path             -- its reconstructed origin path

Note the program reports a whole-Y' attribution; it has no notion of a junction part-way
through an element, so a partial recombinant is reported as if the entire Y' came from the
donor.

Usage: add_program_calls.py <annotated_summary.tsv> [...]
"""
import csv, glob, os, sys

RESULTS = '/home/andrey/argon_scratch/telo_sra_runs/results'
FIELDS = [('program_assigned_yprime', 'y_prime_observed_array'),
          ('program_status',          'y_prime_recombination_status'),
          ('program_source',          'recombination_source'),
          ('program_mechanism',       'recombination_mechanism'),
          ('program_path',            'y_prime_path')]

cache = {}
def lookup(sample, chr_end, read_id):
    key = (sample, chr_end)
    if key not in cache:
        f = f'{RESULTS}/{sample}__elemYPfix/_pipeline/recombination/{sample}_{chr_end}_features.tsv'
        d = {}
        if os.path.exists(f):
            for r in csv.DictReader(open(f), delimiter='\t'): d[r['read_id']] = r
        cache[key] = d
    return cache[key].get(read_id, {})

for path in sys.argv[1:]:
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    if not rows: continue
    for r in rows:
        src = lookup(r['sample'], r['chr_end'], r['read_id'])
        for new, old in FIELDS:
            r[new] = src.get(old, '')
    cols = list(rows[0].keys())
    # put the program columns straight after the expected/donor element columns
    for new, _ in FIELDS:
        if new in cols: cols.remove(new)
    anchor = cols.index('donor_group') + 1 if 'donor_group' in cols else len(cols)
    cols = cols[:anchor] + [n for n, _ in FIELDS] + cols[anchor:]
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
        for r in rows: w.writerow(r)
    print(f'  {path}  (+{len(FIELDS)} columns, {len(rows)} rows)')
