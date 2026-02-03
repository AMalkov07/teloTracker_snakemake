import pandas as pd
import os
import sys
from repeatmasker_utils import load_and_aggregate_repeatmasker_results

print("Starting make_x_element_ends_pairs_repeatmasker_tsv.py")

# Parse command line arguments
if len(sys.argv) != 7:
    print(f"Usage: {sys.argv[0]} <repeatmasker_dir> <strain> <y_prime_probe_file> <output_all> <output_good> <output_gained>")
    print(f"Got {len(sys.argv)} arguments: {sys.argv}")
    sys.exit(1)

repeatmasker_results_dir = sys.argv[1]
strain = sys.argv[2]
y_prime_probe_file = sys.argv[3]
output_all = sys.argv[4]
output_good = sys.argv[5]
output_gained = sys.argv[6]

print(f'Reading from directory: {repeatmasker_results_dir}')
print(f'Strain: {strain}')
print(f'Y prime probe file: {y_prime_probe_file}')


# Chr_end extractor function for X element ends repeatmasker files
def extract_chr_end_x_element(filename):
    """Extract chr_end from X element repeatmasker filename."""
    # File is like: results/BASE/paired_x_element_ends_repeatmasker_results/BASE_chrXY_and_IDZ_x_element_ends_repeatmasker_results.ssv
    basename = os.path.basename(filename)
    # Remove suffix and extract chr_end
    filename_parts = basename.replace('_x_element_ends_repeatmasker_results_corrected.ssv', '')
    # Split by underscore and get the chr part (e.g., chr7L from BASE_chr7L_and_ID3)
    parts = filename_parts.split('_')
    # Find the part that starts with 'chr'
    chr_end_of_file = None
    for part in parts:
        if part.startswith('chr'):
            chr_end_of_file = part
            break

    if chr_end_of_file is None:
        print(f'Warning: Could not extract chr_end from {basename}')

    return chr_end_of_file


# Load and aggregate all RepeatMasker results using shared utility
df_all_repeatmakser_results = load_and_aggregate_repeatmasker_results(
    repeatmasker_results_dir,
    'x_element_ends_repeatmasker_results.ssv',
    chr_end_extractor=extract_chr_end_x_element
)

# Rename the chr_end column to match original script's naming
if 'chr_end' in df_all_repeatmakser_results.columns:
    df_all_repeatmakser_results.rename(columns={'chr_end': 'original_chr_end_anchor'}, inplace=True)

df_all_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)

print(df_all_repeatmakser_results)

# Save all repeatmasker results
df_all_repeatmakser_results.to_csv(output_all, sep='\t')

################## Testing these cutoff values for a good match ##################
df_good_repeatmakser_results = df_all_repeatmakser_results[df_all_repeatmakser_results['sub_match'] == False]

df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['SW_score'] >= 500]
df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['divergence_percent'] <= 2]

# Read y_prime_probe file
df_filter = pd.read_csv(y_prime_probe_file, sep='\t')
df_filter = df_filter.dropna(subset=["repeat_length"])
df_filter = df_filter[df_filter['Adapter_After_Telomere'] == True]

good_reads = df_filter['read_id'].to_list()
df_good_repeatmakser_results['good_read'] = df_good_repeatmakser_results['read_id'].apply(lambda x: x in good_reads)
print(df_good_repeatmakser_results['good_read'].value_counts())

df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['good_read'] == True]

df_good_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
df_good_repeatmakser_results.to_csv(output_good, sep='\t')

###############

df_filter = df_filter[df_filter['delta_y_prime_sign'] == '+']

reads_with_gain_of_y_prime = df_filter['read_id'].to_list()

df_good_repeatmakser_results['gained_y'] = df_good_repeatmakser_results['read_id'].apply(lambda x: x in reads_with_gain_of_y_prime)
print(df_good_repeatmakser_results['gained_y'].value_counts())

df_good_gained_y_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['gained_y'] == True]

df_good_gained_y_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
df_good_gained_y_repeatmakser_results.to_csv(output_gained, sep='\t')

print(f"Finished. Outputs saved to:")
print(f"  - {output_all}")
print(f"  - {output_good}")
print(f"  - {output_gained}")
