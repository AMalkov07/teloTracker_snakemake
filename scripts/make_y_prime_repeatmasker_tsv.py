import pandas as pd
import os
import sys
from repeatmasker_utils import load_and_aggregate_repeatmasker_results

print("Starting make_y_prime_repeatmasker_tsv.py")

# Parse command line arguments
if len(sys.argv) != 6:
    print(f"Usage: {sys.argv[0]} <repeatmasker_dir> <y_prime_probe_file> <output_all_repeatmasker> <output_good_end> <output_gained_y>")
    print(f"Got {len(sys.argv)} arguments: {sys.argv}")
    sys.exit(1)

repeatmasker_dir = sys.argv[1]
y_prime_probe_file = sys.argv[2]
output_all_repeatmasker = sys.argv[3]
output_good_end = sys.argv[4]
output_gained_y = sys.argv[5]

print(f'Reading from directory: {repeatmasker_dir}')
print(f'Y prime probe file: {y_prime_probe_file}')
print(f'Output files: {output_all_repeatmasker}, {output_good_end}, {output_gained_y}')


# Chr_end extractor function for Y' repeatmasker files
def extract_chr_end_y_prime(filename):
    """Extract chr_end from Y prime repeatmasker filename."""
    # File is like: results/BASE/read_repeatmasker_results/BASE_ANCHOR_chrXY_repeatmasker_results_corrected.ssv
    basename = os.path.basename(filename)
    # Remove the suffix
    filename_parts = basename.replace('_repeatmasker_results_corrected.ssv', '')
    # Get the last part which should be chrXY (e.g., chr1L, chr2R)
    chr_end = filename_parts.split('_')[-1]
    return chr_end


# Load and aggregate all RepeatMasker results using shared utility
df_all_repeatmakser_results = load_and_aggregate_repeatmasker_results(
    repeatmasker_dir,
    'repeatmasker_results.ssv',
    chr_end_extractor=extract_chr_end_y_prime
)

print(df_all_repeatmakser_results)

# Save all repeatmasker results
df_all_repeatmakser_results.to_csv(output_all_repeatmasker, sep='\t')

###############

# Read the y_prime_probe file
df_filter = pd.read_csv(y_prime_probe_file, sep='\t')

df_filter = df_filter.dropna(subset=["repeat_length"])
df_filter = df_filter[df_filter['Adapter_After_Telomere'] == True]

df_strand = df_filter[['read_id', 'Repeat_Type']]
df_all_repeatmakser_results = df_all_repeatmakser_results.merge(df_strand, how='left', on='read_id')

def determine_telomere_side(chr_end, strand):
    if ('L' in chr_end and strand == 'AC') or ('R' in chr_end and strand == 'TG'):
        telo_side = 'beginning'
    else:
        telo_side = 'end'
    return telo_side

df_all_repeatmakser_results['telomere_side'] = df_all_repeatmakser_results.apply(
    lambda row: determine_telomere_side(row["chr_end"], row["Repeat_Type"]), axis=1)

reads_with_good_ends = df_filter['read_id'].to_list()

print(f'Reads with good ends: {len(reads_with_good_ends)}')

######

df_all_repeatmakser_results['good_ends'] = df_all_repeatmakser_results['read_id'].apply(
    lambda x: x in reads_with_good_ends)

print(df_all_repeatmakser_results['good_ends'].value_counts())

df_good_end_y_repeatmakser_results = df_all_repeatmakser_results[
    df_all_repeatmakser_results['good_ends'] == True]
df_good_end_y_repeatmakser_results.to_csv(output_good_end, sep='\t')

######

df_filter = df_filter[df_filter['delta_y_prime_sign'] == '+']

reads_with_gain_of_y_prime = df_filter['read_id'].to_list()

print(f'Reads with gain of Y prime: {len(reads_with_gain_of_y_prime)}')

df_all_repeatmakser_results['gained_y'] = df_all_repeatmakser_results['read_id'].apply(
    lambda x: x in reads_with_gain_of_y_prime)

print(df_all_repeatmakser_results['gained_y'].value_counts())

df_gained_y_repeatmakser_results = df_all_repeatmakser_results[
    df_all_repeatmakser_results['gained_y'] == True]
df_gained_y_repeatmakser_results.to_csv(output_gained_y, sep='\t')

print(f"Finished. Outputs saved to:")
print(f"  - {output_all_repeatmasker}")
print(f"  - {output_good_end}")
print(f"  - {output_gained_y}")
