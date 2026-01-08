import pandas as pd
import os
import sys
import math


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

repeatmasker_results_dir_list = os.listdir(repeatmasker_results_dir)

all_results_files = [f'{repeatmasker_results_dir}/{f}' for f in repeatmasker_results_dir_list if f.startswith('dorado') and f.endswith('results.ssv')]

#print(len(all_results_files))

df_all_repeatmakser_results = pd.DataFrame()
for repeatmasker_results_file in all_results_files:
    
    # Need to remove lead spaces from the end_file that come when SW_score < 10,000 such that a space is in the start of the line
    corrected_repeatmasker_results_file = f"{repeatmasker_results_file.split('results.ssv')[0]}results_corrected.ssv"
    
    if os.path.isfile(corrected_repeatmasker_results_file) == 'fix_later':
        pass
    else:
        with open(repeatmasker_results_file, "r") as original_file:
            corrected_data = []
            for line_number,line in enumerate(original_file.readlines()):
                if line_number == 0:
                    corrected_data.append(line)
                else:
                    line = line.strip()
                    if line[-1] != "*":
                        line = f'{line} -'
                    corrected_data.append(line)
            
        with open(corrected_repeatmasker_results_file, "w") as corrected_file:
            for fixed_line in corrected_data:
                corrected_file.write(f'{fixed_line}\n')
    
    try:
        df_single_repeatmasker_results = pd.read_csv(corrected_repeatmasker_results_file, sep=r"\s+")
        chr_end_of_file = corrected_repeatmasker_results_file.split("_repeatmasker_results_corrected.ssv")[0]
        chr_end_of_file = chr_end_of_file.split("rejection_")[1]
        chr_end_of_file = chr_end_of_file.split("_spacer")[0]
        chr_end_of_file = chr_end_of_file.split("_and_ID")[0]
        df_single_repeatmasker_results["original_chr_end_anchor"] = chr_end_of_file
        df_single_repeatmasker_results = df_single_repeatmasker_results.dropna(axis=1, how='all')
        df_all_repeatmakser_results = pd.concat([df_all_repeatmakser_results, df_single_repeatmasker_results])
    except:
        print(f'Error in {corrected_repeatmasker_results_file}')
        pass

# Replace the "*" values with True and the "-" with False
df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '*', 'sub_match'] = True
df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '-', 'sub_match'] = False

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


