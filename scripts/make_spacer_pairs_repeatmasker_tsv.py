import pandas as pd
import os
import sys
import math
from repeatmasker_utils import load_and_aggregate_repeatmasker_results


# Parse command-line arguments
repeatmasker_results_dir = sys.argv[1]
strain = sys.argv[2]
anchors_and_distances_bed = sys.argv[3]
final_features_bed = sys.argv[4]
y_prime_probe_tsv = sys.argv[5]
output_all_tsv = sys.argv[6]
output_good_tsv = sys.argv[7]
output_good_gained_tsv = sys.argv[8]

print("Starting make_spacer_pairs_repeatmasker_tsv.py")

print(f'Opening spacer repeatmasker results from {repeatmasker_results_dir}...')

input_anchor_distances_bed_file = anchors_and_distances_bed
bed_file_of_features = final_features_bed


# Chr_end extractor function for spacer repeatmasker files
def extract_chr_end_spacer(filename):
    """Extract chr_end from spacer repeatmasker filename."""
    # Hardcoded parsing based on original script (lines 55-58)
    chr_end_of_file = filename.split("_repeatmasker_results_corrected.ssv")[0]
    chr_end_of_file = chr_end_of_file.split("rejection_")[1]
    chr_end_of_file = chr_end_of_file.split("_spacer")[0]
    chr_end_of_file = chr_end_of_file.split("_and_ID")[0]
    return chr_end_of_file


# Load and aggregate all RepeatMasker results using shared utility
df_all_repeatmakser_results = load_and_aggregate_repeatmasker_results(
    repeatmasker_results_dir,
    'results.ssv',  # Note: original script used broader pattern
    chr_end_extractor=extract_chr_end_spacer
)

# Rename the chr_end column to match original script's naming
if 'chr_end' in df_all_repeatmakser_results.columns:
    df_all_repeatmakser_results.rename(columns={'chr_end': 'original_chr_end_anchor'}, inplace=True)

print(df_all_repeatmakser_results)


# Read the BED file into a DataFrame
bed_columns = ["Chr", "Start", "End", "Name", "Strand", "Length"]
df_bed = pd.read_csv(bed_file_of_features, sep='\t', names=bed_columns)


distances_to_anchor_dict = {}

for chr_end in df_all_repeatmakser_results['original_chr_end_anchor'].unique():

    chopped_segment_length = 250

    distance_to_anchor_start = df_bed[df_bed["Name"] == f'{chr_end}_space_between_anchor']['Length'].iloc[0]
    distance_to_anchor_end = distance_to_anchor_start + df_bed[df_bed["Name"] == f'{chr_end}_anchor']['Length'].iloc[0]

    section_of_anchor_start = (distance_to_anchor_start / chopped_segment_length)
    section_of_anchor_end = (distance_to_anchor_end / chopped_segment_length)

    distances_to_anchor_dict[chr_end] = (section_of_anchor_start, section_of_anchor_end)

def label_anchor_on_chopped_segments(chr_end, section_number, row):

    section_of_anchor_start, section_of_anchor_end = distances_to_anchor_dict[chr_end]

    # Checks low to high if the section number is in the range of the anchor
    if section_number < math.floor(section_of_anchor_start):
        anchor_label = 'Before Anchor'
    elif section_number == math.floor(section_of_anchor_start):
        if section_number == section_of_anchor_start:
            anchor_label = 'Partial Before Anchor'
        else:
            anchor_label = 'Is Anchor'
    elif section_number <= section_of_anchor_end:
        anchor_label = 'Is Anchor'
    elif section_number < math.ceil(section_of_anchor_end):
        anchor_label = 'Partial After Anchor'
    else:
        anchor_label = 'After Anchor'

    return anchor_label


df_all_repeatmakser_results['anchor_label'] = df_all_repeatmakser_results.apply(
    lambda row: label_anchor_on_chopped_segments(row['original_chr_end_anchor'], row['section_number'], row), axis=1)


df_all_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)

df_all_repeatmakser_results.to_csv(output_all_tsv, sep='\t')

################## Testing these cutoff values for a good match ##################
df_good_repeatmakser_results = df_all_repeatmakser_results[df_all_repeatmakser_results['sub_match'] == False]

df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['SW_score'] >= 500]
df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['divergence_percent'] <= 2]

df_filter = pd.read_csv(y_prime_probe_tsv, sep='\t')
df_filter = df_filter.dropna(subset=["repeat_length"])
df_filter = df_filter[df_filter['Adapter_After_Telomere'] == True]

good_reads = df_filter['read_id'].to_list()
df_good_repeatmakser_results['good_read'] = df_good_repeatmakser_results['read_id'].apply(lambda x: x in good_reads)
print(df_good_repeatmakser_results['good_read'].value_counts())

df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['good_read'] == True]

df_good_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
df_good_repeatmakser_results.to_csv(output_good_tsv, sep='\t')

###############

df_filter = df_filter[df_filter['delta_y_prime_sign'] == '+']

reads_with_gain_of_y_prime = df_filter['read_id'].to_list()

df_good_repeatmakser_results['gained_y'] = df_good_repeatmakser_results['read_id'].apply(lambda x: x in reads_with_gain_of_y_prime)
print(df_good_repeatmakser_results['gained_y'].value_counts())

df_good_gained_y_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['gained_y'] == True]

df_good_gained_y_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
df_good_gained_y_repeatmakser_results.to_csv(output_good_gained_tsv, sep='\t')
