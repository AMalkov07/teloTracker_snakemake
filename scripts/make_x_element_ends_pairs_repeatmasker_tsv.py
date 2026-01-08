import pandas as pd
import os
import sys

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

# Not used in the script but kept for reference
input_anchor_distances_bed_file = f'{strain}_anchors_and_distances.bed'

repeatmasker_results_dir_list = os.listdir(repeatmasker_results_dir)

# Look for files ending in results.ssv (not just those starting with 'dorado')
all_results_files = [f'{repeatmasker_results_dir}/{f}' for f in repeatmasker_results_dir_list 
                     if f.endswith('x_element_ends_repeatmasker_results.ssv')]

print(f'Found {len(all_results_files)} x element repeatmasker result files')

df_all_repeatmakser_results = pd.DataFrame()
for repeatmasker_results_file in all_results_files:
    
    # Need to remove lead spaces from the file that come when SW_score < 10,000
    corrected_repeatmasker_results_file = f"{repeatmasker_results_file.split('results.ssv')[0]}results_corrected.ssv"
    
    if os.path.isfile(corrected_repeatmasker_results_file) == 'fix_later':
        pass
    else:
        with open(repeatmasker_results_file, "r") as original_file:
            corrected_data = []
            for line_number, line in enumerate(original_file.readlines()):
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
        
        # Extract chr_end from filename
        # File is like: results/BASE/paired_x_element_ends_repeatmasker_results/BASE_chrXY_and_IDZ_x_element_ends_repeatmasker_results.ssv
        filename = os.path.basename(corrected_repeatmasker_results_file)
        # Remove suffix and extract chr_end
        filename_parts = filename.replace('_x_element_ends_repeatmasker_results_corrected.ssv', '')
        # Split by underscore and get the chr part (e.g., chr7L from BASE_chr7L_and_ID3)
        parts = filename_parts.split('_')
        # Find the part that starts with 'chr'
        chr_end_of_file = None
        for part in parts:
            if part.startswith('chr'):
                chr_end_of_file = part
                break
        
        if chr_end_of_file is None:
            print(f'Warning: Could not extract chr_end from {filename}')
            continue
            
        df_single_repeatmasker_results["original_chr_end_anchor"] = chr_end_of_file
        df_single_repeatmasker_results = df_single_repeatmasker_results.dropna(axis=1, how='all')
        df_all_repeatmakser_results = pd.concat([df_all_repeatmakser_results, df_single_repeatmasker_results])
    except Exception as e:
        print(f'Error in {corrected_repeatmasker_results_file}: {e}')
        pass

# Replace the "*" values with True and the "-" with False
df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '*', 'sub_match'] = True
df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '-', 'sub_match'] = False

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