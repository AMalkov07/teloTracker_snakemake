import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.ticker import FuncFormatter

# Usage: y_prime_analysis.py <base_name> <output_dir>
# Note: anchor_set is still needed to locate blast results; kept as argv[2] used internally
# Shell calls: python y_prime_analysis.py $BASE_NAME $OUTPUT_DIR
# We derive anchor_set from the top_matches file by scanning for it, or keep it explicit.
# Per the shell script, y_prime_analysis.py receives: base_name output_dir
# anchor_set is not passed here — we read the post_telo_trimming tsv which doesn't need it directly,
# but we need anchor_set to find top_matches. We infer anchor_set from available files.

base_name  = sys.argv[1]
output_dir = sys.argv[2]

print("Starting y_prime_analysis.py")

blast_dir      = os.path.join(output_dir, 'blast')
y_prime_dir    = os.path.join(output_dir, 'blast', 'y_prime_probe')
graphs_dir     = os.path.join(output_dir, 'graphs')
figures_dir    = os.path.join(output_dir, 'graphs', 'y_prime_figures')
stats_dir      = os.path.join(output_dir, 'graphs', 'stats_for_y_primes')

os.makedirs(figures_dir, exist_ok=True)
os.makedirs(stats_dir,   exist_ok=True)

telomere_repeat_results_input = os.path.join(output_dir, f'{base_name}_post_telo_trimming.tsv')
output_tsv_f_name             = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')


def read_direction(match_start_on_anchor, match_end_on_anchor):
    return 'forward' if match_start_on_anchor < match_end_on_anchor else 'reverse'

def repeat_type_of_read(l_end_chr, alignment_direction):
    if l_end_chr:
        return 'AC' if alignment_direction == 'forward' else 'TG'
    else:
        return 'TG' if alignment_direction == 'forward' else 'AC'

def match_length_calc(match_start_on_anchor, match_end_on_anchor):
    return abs(match_start_on_anchor - match_end_on_anchor)

def y_prime_change_calc(anchor_name, base_name):
    if '6991' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 7,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7021' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 8,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 6,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7093' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 1, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 2, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 1,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 4, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7154' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 0,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 5,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 1, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7172' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 8,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 3, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 3, 'chr16R_anchor': 1
        }
    elif '7174' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 9,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 7,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7250' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 1, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 10,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 6, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7302' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 6,
            'chr13L_anchor': 4, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7321' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 4,
            'chr5L_anchor': 1, 'chr5R_anchor': 5, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 5, 'chr8L_anchor': 5, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 1, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 4,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 7, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7323' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 9,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 1,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 4, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7324' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 14,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 0, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 5, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 1, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 10,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 13, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif '7372' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 7,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif 'KRLT3' in base_name or '7637' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 1,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 4,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 2, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 7,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }
    elif 'tlc1' in base_name:
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 1, 'chr2L_anchor': 0, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 0,
            'chr5L_anchor': 1, 'chr5R_anchor': 0, 'chr6L_anchor': 1, 'chr6R_anchor': 1,
            'chr7L_anchor': 0, 'chr7R_anchor': 0, 'chr8L_anchor': 0, 'chr8R_anchor': 1,
            'chr9L_anchor': 0, 'chr9R_anchor': 1, 'chr10L_anchor': 0, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 1, 'chr12L_anchor': 0, 'chr12R_anchor': 0,
            'chr13L_anchor': 0, 'chr13R_anchor': 0, 'chr14L_anchor': 0, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 0, 'chr16L_anchor': 0, 'chr16R_anchor': 1
        }
    else:  # Default to 6991 (WildType)
        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0,
            'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0,
            'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0,
            'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 7,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1,
            'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }

    return y_prime_numbers[anchor_name]


def assign_delta_y_prime_sign(y_primes_relative_to_ref):
    if y_primes_relative_to_ref > 0:   return '+'
    elif y_primes_relative_to_ref < 0: return '-'
    else:                               return 'same'

def assign_reference_y_prime_end_status(reference_y_primes):
    if reference_y_primes == 0:    return 'Reference Y Primes = 0'
    if reference_y_primes == 1:    return 'Reference Y Primes = 1'
    elif reference_y_primes >= 2:  return 'Reference Y Primes = 2+'

def assign_read_y_prime_end_status(y_prime_probe_count):
    if y_prime_probe_count == 0:    return 'Read Y Primes = 0'
    if y_prime_probe_count == 1:    return 'Read Y Primes = 1'
    elif y_prime_probe_count >= 2:  return 'Read Y Primes = 2+'


# ---- Load all per-chromosome Y' probe blast results ----
probe_files = [f for f in os.listdir(y_prime_dir) if f.endswith('blasted_probe.tsv')]
df = pd.DataFrame()
for f in probe_files:
    df_temp = pd.read_csv(os.path.join(y_prime_dir, f), sep='\t')
    df = pd.concat([df, df_temp])

# Save combined probe results
combined_probe_out = os.path.join(y_prime_dir, f'all_{base_name}_probe_matches.tsv')
df.to_csv(combined_probe_out, sep='\t')

column_type_dict = {
    'read_id': str, 'total_read_length': int, 'read_bp_used_for_match': int,
    'match_start_on_read': int, 'match_end_on_read': int, 'anchor_name': str,
    'total_anchor_length': int, 'match_start_on_anchor': int, 'match_end_on_anchor': int,
    'pident': float, 'bitscore': float, 'evalue': float
}
df = df.astype(column_type_dict)

df_filtered = df[df['read_bp_used_for_match'] > df['total_anchor_length'] / 2]
df_filtered = df_filtered.sort_values(['read_id', 'match_start_on_read', 'bitscore'], ascending=False)
df_filtered['match_length'] = df_filtered.apply(
    lambda r: match_length_calc(r['match_start_on_anchor'], r['match_end_on_anchor']), axis=1)

chr_l_end_list = [f'chr{n}L_anchor' for n in range(1, 17)]
df_filtered['l_end_chr']          = df_filtered['anchor_name'].apply(lambda x: x in chr_l_end_list)
df_filtered['alignment_direction'] = df_filtered.apply(
    lambda r: read_direction(r['match_start_on_anchor'], r['match_end_on_anchor']), axis=1)
df_filtered['repeat_type']        = df_filtered.apply(
    lambda r: repeat_type_of_read(r['l_end_chr'], r['alignment_direction']), axis=1)

read_y_prime_counts      = df_filtered['read_id'].value_counts().reset_index()
read_y_prime_counts.columns = ['read_id', 'y_prime_probe_count']

df_repeat_results = pd.read_csv(telomere_repeat_results_input, sep='\t', index_col=0)
df_repeat_results['read_id'] = df_repeat_results['read_id'].astype(str)

df_combined = pd.merge(df_repeat_results, read_y_prime_counts, on='read_id', how='outer')

df_combined['y_prime_probe_count'] = df_combined['y_prime_probe_count'].fillna(0)
df_combined = df_combined.dropna(subset=['anchor_name'])

df_combined['reference_y_primes']           = df_combined['anchor_name'].apply(lambda a: y_prime_change_calc(a, base_name))
df_combined['y_primes_relative_to_ref']     = df_combined['y_prime_probe_count'] - df_combined['reference_y_primes']
df_combined['delta_y_prime_sign']           = df_combined['y_primes_relative_to_ref'].apply(assign_delta_y_prime_sign)
df_combined['reference_y_prime_end_status'] = df_combined['reference_y_primes'].apply(assign_reference_y_prime_end_status)
df_combined['read_y_prime_end_status']      = df_combined['y_prime_probe_count'].apply(assign_read_y_prime_end_status)
df_combined['y_prime_type_delta_y_prime_sign'] = df_combined['delta_y_prime_sign'] + ' ' + df_combined['reference_y_prime_end_status']

print(df_combined)

chr_end_list = ['1L','1R','2L','2R','3L','3R','4L','4R','5L','5R','6L','6R','7L','7R','8L','8R',
                '9L','9R','10L','10R','11L','11R','12L','12R','13L','13R','14L','14R','15L','15R','16L','16R']
chr_end_sorter = {c: i for i, c in enumerate(chr_end_list)}

df_combined['chr_end']      = df_combined['anchor_name'].apply(lambda a: a.strip('chr _anchor'))
df_combined['chr_end_rank'] = df_combined['chr_end'].map(chr_end_sorter)
df_combined.sort_values('chr_end_rank', ascending=True, inplace=True)

df_combined.to_csv(output_tsv_f_name, sep='\t', index=False)
print(f'Written: {output_tsv_f_name}')


# ---- Good telomere reads filter   NOTE: 'repeat_length' >= 30 can be adjusted as needed (e.g. 'repeat length' > 40) ----
df_good_reads = df_combined.loc[(df_combined['repeat_length'] >= 30) & (df_combined['Adapter_After_Telomere'] == True)]
print('Reads in df_good_reads:')
print(df_good_reads)
print(df_good_reads['y_prime_probe_count'].value_counts())
print(df_good_reads['y_primes_relative_to_ref'].value_counts())
print(df_good_reads['anchor_name'].value_counts())

# ---- Stats file output ----
stats_neutral  = df_good_reads[df_good_reads['delta_y_prime_sign'] == 'same']['y_primes_relative_to_ref'].describe()
stats_positive = df_good_reads[df_good_reads['delta_y_prime_sign'] == '+']['y_primes_relative_to_ref'].describe()
stats_negative = df_good_reads[df_good_reads['delta_y_prime_sign'] == '-']['y_primes_relative_to_ref'].describe()

output_stats_f_name = os.path.join(stats_dir, f'{base_name}_stats_y_prime.txt')
print(f'\nWriting {output_stats_f_name}')
with open(output_stats_f_name, 'w') as f:
    f.write(f'Neutral Y prime ends:\n{stats_neutral}\n'
            f'Positive Y prime ends:\n{stats_positive}\n'
            f'Negative Y prime ends:\n{stats_negative}\n')


# ---------------------------------------------------------------------------
# Y prime plots
# ---------------------------------------------------------------------------

def y_prime_count_violin_strip_plot(dataframe, plot_scale=(22, 9)):

    sns.set(style="whitegrid")

    chr_order = [f'{n}{arm}' for n in range(1, 17) for arm in ('L', 'R')]
    chr_list  = [c for c in chr_order if c in dataframe['chr_end'].unique()]

    if not chr_list:
        print(f'Warning: no chr_end values found in data, skipping y_prime_count_violin_strip_plot')
        return

    full_color_palette = ['#641E16', '#512E5F', '#4A235A', '#154360', '#2E86C1', '#0E6251', '#0B5345', '#145A32',
                          '#7D6608', '#7E5109', '#D35400', '#626567', '#4D5656', '#212F3C', '#17202A', '#F93409',
                          '#FFED24', '#CD6155', '#D7BDE2', '#BB8FCE', '#A9CCE3', '#AED6F1', '#A3E4D7', '#45B39D',
                          '#229954', '#F7DC6F', '#F5B041', '#E59866', '#D7DBDD', '#95A5A6', '#5D6D7E', '#ABB2B9',
                          '#FFC2B4', '#FFFEB4']
    chr_color_dict = dict(zip(chr_order, full_color_palette))
    chr_color_code = [chr_color_dict[c] for c in chr_list]
    sns.set_palette(chr_color_code)

    fig, ax = plt.subplots(figsize=plot_scale)
    sns.violinplot(x='chr_end', hue='chr_end', order=chr_list, y='y_primes_relative_to_ref', data=dataframe, gridsize=1000, cut=0, palette=chr_color_code,legend=False, ax=ax)
    sns.stripplot(x='chr_end', order=chr_list, y='y_primes_relative_to_ref', data=dataframe, linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="#FF2400", ax=ax)
    
    total_reads = len(dataframe)
    ax.set_title(f'{base_name} Delta Y Primes (N = {total_reads} Reads, Ave. Y Prime to ref = {dataframe["y_primes_relative_to_ref"].mean():.3f})',
                 fontweight="bold", fontsize=20, color='k', pad=15)
    plt.xlabel("Chromosome End", fontweight="bold", fontsize=20)
    plt.ylabel("Delta Y Primes", fontweight="bold", fontsize=30)
    plt.xticks(fontweight="bold", fontsize=13)
    plt.yticks(fontweight="bold", fontsize=20)

    """
    for i, anchor in enumerate(chr_list):
        df_anchor = dataframe[dataframe['chr_end'] == anchor]
        num_reads = df_anchor['y_primes_relative_to_ref'].count()
        plt.text(i, ax.get_ylim()[1] - 0.4, f'{num_reads * 100 / total_reads:.0f}%',
                 ha='center', fontsize=10, fontweight="bold")
        ref_y = df_anchor['reference_y_primes'].iloc[0] if len(df_anchor) > 0 else 'N/A'
        plt.text(i, ax.get_ylim()[0] + 0.2, f'Ref={ref_y}',
                 ha='center', fontsize=10, fontweight="bold")
    """

    save_file_name = os.path.join(figures_dir, f'{base_name}_delta_y_primes.png')
    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=300, format="png")
    plt.close()


def delta_sign_violin_strip_plot(dataframe, plot_scale=(16, 9)):

    sns.set(style="whitegrid")

    fig, ax = plt.subplots(figsize=plot_scale)
    sns.violinplot(x='delta_y_prime_sign', hue='delta_y_prime_sign', hue_order=["-", "same", "+"], order=["-", "same", "+"], y='y_primes_relative_to_ref',
                   data=dataframe, gridsize=1000, cut=0, palette=['tab:red', 'tab:blue', 'tab:green'], legend=False, ax=ax)
    ax.set(ylim=(-10, 20))
    sns.stripplot(x='delta_y_prime_sign', order=["-", "same", "+"], y='y_primes_relative_to_ref', data=dataframe,
                  linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="tab:gray", ax=ax)

    total_reads    = len(dataframe)
    total_same     = len(dataframe[dataframe['delta_y_prime_sign'] == "same"])
    total_positive = len(dataframe[dataframe['delta_y_prime_sign'] == "+"])
    total_negative = len(dataframe[dataframe['delta_y_prime_sign'] == "-"])

    ax.set_title(f'{base_name} Delta Y Primes (N = {total_reads} Reads, -:{total_negative}, =:{total_same}, & +:{total_positive})',
                 fontweight="bold", fontsize=20, color='k', pad=15)
    plt.xlabel("Sign of Delta Y Primes", fontweight="bold", fontsize=20)
    plt.ylabel("Number of Delta Y Primes", fontweight="bold", fontsize=30)
    plt.xticks(fontweight="bold", fontsize=13)
    plt.yticks(fontweight="bold", fontsize=20)

    save_file_name = os.path.join(figures_dir, f'{base_name}_sign_delta_y_primes.png')
    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=300, format="png")
    plt.close()


# ---- Call plots on good telomere reads ----
df_y_prime = df_good_reads.dropna(subset=['y_primes_relative_to_ref', 'chr_end'])

y_prime_count_violin_strip_plot(df_y_prime)

delta_sign_violin_strip_plot(df_y_prime)