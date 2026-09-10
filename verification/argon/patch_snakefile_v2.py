#!/usr/bin/env python3
"""Patch a run-dir _pipeline/Snakefile so rule recombination_analyze passes the
v2 attribution inputs/flags to analyze_features.py:
  + inputs telo_tsv / probe_tsv   (telomere-end confirmation, Y' probe count)
  + params attribution_mode, y_prime_id_level
  + flags  --telo-tsv --probe-tsv --attribution-mode --y-prime-id-level
Idempotent; also upgrades a Snakefile patched by the earlier version (no id-level).
Usage: patch_snakefile_v2.py <path/to/Snakefile>"""
import sys

p = sys.argv[1]
s = open(p).read()
if '--y-prime-id-level' in s:
    print('already patched:', p)
    sys.exit(0)

IN_OLD = '        day0_ref    = DAY0_REF,\n        y_prime_lib = Y_PRIME_LIB,\n    output:\n'
IN_NEW = ('        day0_ref    = DAY0_REF,\n        y_prime_lib = Y_PRIME_LIB,\n'
          '        telo_tsv    = f"{RESULTS}/{BASE}_post_telo_trimming.tsv",\n'
          '        probe_tsv   = f"{RESULTS}/{BASE}_post_y_prime_probe.tsv",\n    output:\n')
PAR_OLD = '        x_element_lib     = X_ELEM_LIB,\n    threads: workflow.cores\n'
PAR_NEW = ('        x_element_lib     = X_ELEM_LIB,\n'
           '        attribution_mode  = config.get("attribution_mode", "v2"),\n'
           '        y_prime_id_level  = config.get("y_prime_id_level", "family"),\n    threads: workflow.cores\n')
FL_OLD = '            --x-element-lib     {params.x_element_lib} \\\n            --chr-end           {params.chr_end} \\\n'
FL_NEW = ('            --x-element-lib     {params.x_element_lib} \\\n'
          '            --telo-tsv          {input.telo_tsv} \\\n'
          '            --probe-tsv         {input.probe_tsv} \\\n'
          '            --attribution-mode  {params.attribution_mode} \\\n'
          '            --y-prime-id-level  {params.y_prime_id_level} \\\n'
          '            --chr-end           {params.chr_end} \\\n')

if '--attribution-mode' in s:
    # patched by the earlier version: add only the id-level param + flag
    a = '        attribution_mode  = config.get("attribution_mode", "v2"),\n'
    b = '            --attribution-mode  {params.attribution_mode} \\\n'
    if s.count(a) != 1 or s.count(b) != 1:
        print('PATTERN NOT FOUND (upgrade):', p); sys.exit(1)
    s = s.replace(a, a + '        y_prime_id_level  = config.get("y_prime_id_level", "family"),\n', 1)
    s = s.replace(b, b + '            --y-prime-id-level  {params.y_prime_id_level} \\\n', 1)
    open(p, 'w').write(s); print('upgraded (id-level):', p); sys.exit(0)

for old, new in ((IN_OLD, IN_NEW), (PAR_OLD, PAR_NEW), (FL_OLD, FL_NEW)):
    if s.count(old) != 1:
        print('PATTERN NOT FOUND (or not unique):', p, '\n', old); sys.exit(1)
    s = s.replace(old, new, 1)
open(p, 'w').write(s)
print('patched:', p)
