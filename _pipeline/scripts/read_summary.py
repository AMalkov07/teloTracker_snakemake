#!/usr/bin/env python3
"""
read_summary.py — Per-chromosome-end summary of anchored and qualifying reads.

Two counts per chromosome end:

  anchored    Reads carrying that end's anchor. Every row of
              <base>_post_y_prime_probe.tsv is an anchored read, so this is
              simply the row count per chr_end.

  qualifying  Reads with a telomere that was actually found and measured —
              Adapter_After_Telomere is True AND repeat_length >= 30 (only the
              repeat length when the adapter check was off: --no-adapter-filter). This is
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
    python read_summary.py <base_name> <output_dir> [--no-adapter-filter]

The header also gives, per strand, the share of reads with a telomere that have an adapter after it,
and a warning when the AC strand is so low that the reads look adapter-trimmed.

Writes:
    <output_dir>/<base_name>_read_summary.tsv
"""

import glob
import os
import sys

import pandas as pd

from telomere_filters import MIN_TELOMERE_REPEAT, adapter_found, qualifying_mask

# Usage: read_summary.py <base_name> <output_dir> [--no-adapter-filter]
base_name  = sys.argv[1]
output_dir = sys.argv[2]
require_adapter = '--no-adapter-filter' not in sys.argv[3:]   # the adapter check was off

input_tsv  = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')
output_tsv = os.path.join(output_dir, f'{base_name}_read_summary.tsv')

print("Starting read_summary.py")
print(f'Opening {input_tsv}...')

df = pd.read_csv(input_tsv, sep='\t')

# Same filter as single_sample_plots.py df_graph (telomere_filters.py)
is_qualifying = qualifying_mask(df, require_adapter)

# Reads that reach the telomere mostly carry the sequencing adapter after it: about 86 % of AC-strand
# telomere reads in MinKNOW's own output. Far fewer means the adapters were trimmed off (dorado's
# default trimming, or before a public upload), and those reads drop out of the qualifying set for
# that reason alone. Warning only: nothing is changed.
TRIMMED_AC_SHARE = 40.0
adapter_share = {}
if require_adapter and 'Repeat_Type' in df.columns:
    has_telomere = pd.to_numeric(df['repeat_length'], errors='coerce') >= MIN_TELOMERE_REPEAT
    for strand in ('AC', 'TG'):
        sel = has_telomere & (df['Repeat_Type'] == strand)
        if sel.sum():
            adapter_share[strand] = (round(100 * adapter_found(df.loc[sel, 'Adapter_After_Telomere']).mean(), 1),
                                     int(sel.sum()))
trimmed_warning = None
ac = adapter_share.get('AC')
if ac and ac[1] >= 20 and ac[0] < TRIMMED_AC_SHARE:
    trimmed_warning = (f'only {ac[0]}% of AC-strand reads with a telomere have an adapter after it '
                       f'(MinKNOW output: about 86%). The reads look adapter-trimmed (dorado default '
                       f'trimming, or public data), so most telomere reads are not counted as qualifying. '
                       f'Basecall with --trim none, or set adapter_check: false.')

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

# Reads removed before anchoring, so the anchored counts below can be read in context.
# Both removals happen upstream and were previously silent.
#   fold-back : one molecule read twice in inverted orientation (detect_foldback_reads.py)
#   multi-end : a > half-anchor hit to two DIFFERENT ends' anchors -- a ligation chimera.
#               filter_for_reads_with_anchors.py drops these via drop_duplicates(keep=False).
blast_dir = os.path.join(output_dir, 'blast')
n_foldback = None
fb_file = os.path.join(blast_dir, f'{base_name}_foldback_read_ids.txt')
if os.path.exists(fb_file):
    with open(fb_file) as f:
        n_foldback = sum(1 for ln in f if ln.strip())
n_multi_end = None
# glob, since this script is not told the anchor set; only trust a single match
all_matches = glob.glob(os.path.join(blast_dir, f'all_matches_{base_name}_blasted_*.tsv'))
if len(all_matches) == 1:
    am = pd.read_csv(all_matches[0], sep='\t', usecols=['read_id', 'anchor_name'])
    n_multi_end = int((am.groupby('read_id')['anchor_name'].nunique() > 1).sum())

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
    if n_foldback is not None:
        f.write(f'# removed as fold-back (hairpin) reads: {n_foldback}\n')
    if n_multi_end is not None:
        f.write(f'# removed as multi-end chimeras (anchors of two different ends): {n_multi_end}\n')
    if adapter_share:
        f.write('# adapter after the telomere (reads with a telomere): '
                + ', '.join(f'{k} {v[0]}% of {v[1]}' for k, v in adapter_share.items()) + '\n')
    if trimmed_warning:
        f.write(f'# WARNING: {trimmed_warning}\n')
    f.write(f'# anchored: read carries that end\'s anchor\n')
    if require_adapter:
        f.write(f'# qualifying: Adapter_After_Telomere == True AND repeat_length >= {MIN_TELOMERE_REPEAT}\n')
    else:
        f.write(f'# qualifying: repeat_length >= {MIN_TELOMERE_REPEAT} (adapter check off)\n')
    f.write(f'#   (identical to the read set graphed in {base_name}_500bp.png)\n')
    summary.to_csv(f, sep='\t', index=False)

print(summary.to_string(index=False))

if trimmed_warning:
    print(f'\nWARNING: {trimmed_warning}')

low = summary[(summary['chr_end'] != 'ALL') & (summary['qualifying_reads'] < 10)]
if len(low):
    ends = ', '.join(f"{r.chr_end}({r.qualifying_reads})" for r in low.itertuples())
    print(f'\nWarning: {len(low)} chromosome end(s) have fewer than 10 qualifying '
          f'reads and are poorly characterised: {ends}')
