#!/usr/bin/env python3
"""Write a read-name translation table from a FASTQ whose deflines carry a second
name token -- SRA dumps look like '@SRR33298449.1 00967a24-b45e-4e3a-ba71-852b2c6762c8 length=922',
so the pipeline's read_id (first token) can be translated back to the original ONT read UUID.

Usage: make_read_id_map.py <reads.fastq[.gz]> <out.tsv>
Output columns: read_id  original_name  (only deflines with a second token are written)."""
import gzip, sys
src, out = sys.argv[1], sys.argv[2]
op = gzip.open if src.endswith('.gz') else open
n = 0
with op(src, 'rt') as fh, open(out, 'w') as o:
    o.write('read_id\toriginal_name\n')
    for i, line in enumerate(fh):
        if i % 4:
            continue
        parts = line[1:].split()
        if len(parts) >= 2 and not parts[1].startswith('length='):
            o.write(f'{parts[0]}\t{parts[1]}\n'); n += 1
print(f'{out}: {n} reads mapped')
