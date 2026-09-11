#!/usr/bin/env python3
"""
read_summary.py — Per-chromosome-end summary of anchored and qualifying reads.

Two counts per chromosome end:

  anchored    Reads carrying that end's anchor. Every row of
              <base>_post_y_prime_probe.tsv is an anchored read, so this is
              simply the row count per chr_end.

  qualifying  Reads with a telomere that was actually found and measured —
              Adapter_After_Telomere is True AND repeat_length >= 30. This is
              the SAME filter single_sample_plots.py applies to build df_graph,
              so the qualifying total equals the N printed on the
              <base>_500bp.png and <base>_facet_by_chr.png titles, and the
              per-end numbers equal the N behind each facet panel.

The gap between the two matters: an anchored read that dies in the subtelomere
never reaches the telomere, so anchored counts overstate how well an end is
characterised — often by a lot at ends carrying long Y' arrays, where reads run
out of length before the chromosome runs out of sequence. Ends with a low
qualifying count are the ones whose telomere and Y' numbers should not be
trusted, which is exactly what this file is for.

Usage:
    python read_summary.py <base_name> <output_dir>

Writes:
    <output_dir>/<base_name>_read_summary.tsv
"""

import os
import sys

import pandas as pd

base_name  = sys.argv[1]
output_dir = sys.argv[2]

input_tsv  = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')
output_tsv = os.path.join(output_dir, f'{base_name}_read_summary.tsv')

print("Starting read_summary.py")
print(f'Opening {input_tsv}...')

df = pd.read_csv(input_tsv, sep='\t')

# Same filter as single_sample_plots.py df_graph
is_qualifying = (df['repeat_length'] >= 30) & (df['Adapter_After_Telomere'] == True)

chr_order = [f'{n}{arm}' for n in range(1, 17) for arm in ('L', 'R')]
present   = [c for c in chr_order if c in df['chr_end'].unique()]
# keep any end not covered by the 1-16 L/R naming rather than silently dropping it
present  += [c for c in df['chr_end'].dropna().unique() if c not in chr_order]

# Total reads that entered the pipeline, for context on how few end up anchored
total_reads = None
fai = os.path.join(output_dir, f'{base_name}.fasta.fai')
if os.path.exists(fai):
    with open(fai) as f:
        total_reads = sum(1 for _ in f)

rows = []
for chr_end in present:
    sub = df[df['chr_end'] == chr_end]
    anchored   = len(sub)
    qualifying = int(is_qualifying[sub.index].sum())
    rows.append({
        'chr_end':          chr_end,
        'anchored_reads':   anchored,
        'qualifying_reads': qualifying,
        'pct_qualifying':   round(100 * qualifying / anchored, 1) if anchored else 0.0,
    })

total_anchored   = len(df)
total_qualifying = int(is_qualifying.sum())
rows.append({
    'chr_end':          'ALL',
    'anchored_reads':   total_anchored,
    'qualifying_reads': total_qualifying,
    'pct_qualifying':   round(100 * total_qualifying / total_anchored, 1) if total_anchored else 0.0,
})

summary = pd.DataFrame(rows)

print(f'Writing {output_tsv}...')
with open(output_tsv, 'w') as f:
    f.write(f'# sample: {base_name}\n')
    if total_reads is not None:
        f.write(f'# reads after filter_reads.py: {total_reads}\n')
    f.write(f'# anchored: read carries that end\'s anchor\n')
    f.write(f'# qualifying: Adapter_After_Telomere == True AND repeat_length >= 30\n')
    f.write(f'#   (identical to the read set graphed in {base_name}_500bp.png)\n')
    summary.to_csv(f, sep='\t', index=False)

print(summary.to_string(index=False))

low = summary[(summary['chr_end'] != 'ALL') & (summary['qualifying_reads'] < 10)]
if len(low):
    ends = ', '.join(f"{r.chr_end}({r.qualifying_reads})" for r in low.itertuples())
    print(f'\nWarning: {len(low)} chromosome end(s) have fewer than 10 qualifying '
          f'reads and are poorly characterised: {ends}')
