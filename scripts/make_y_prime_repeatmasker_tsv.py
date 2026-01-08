import pandas as pd
import os
import sys

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

# Get all repeatmasker result files
repeatmasker_results_dir_list = os.listdir(repeatmasker_dir)
# Look for files ending in results.ssv (not just those starting with 'dorado')
all_results_files = [f'{repeatmasker_dir}/{f}' for f in repeatmasker_results_dir_list 
                     if f.endswith('repeatmasker_results.ssv')]

print(f'Found {len(all_results_files)} repeatmasker result files')

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
        # File is like: results/BASE/read_repeatmasker_results/BASE_ANCHOR_chrXY_repeatmasker_results_corrected.ssv
        filename = os.path.basename(corrected_repeatmasker_results_file)
        # Remove the suffix
        filename_parts = filename.replace('_repeatmasker_results_corrected.ssv', '')
        # Get the last part which should be chrXY (e.g., chr1L, chr2R)
        chr_end = filename_parts.split('_')[-1]  # Gets the last underscore-separated part
        
        df_single_repeatmasker_results["chr_end"] = chr_end
        df_all_repeatmakser_results = pd.concat([df_all_repeatmakser_results, df_single_repeatmasker_results])
    except Exception as e:
        print(f'Error in {corrected_repeatmasker_results_file}: {e}')
        pass

# Replace the "*" values with True and the "-" with False
df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '*', 'sub_match'] = True
df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '-', 'sub_match'] = False

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
print(f"  - {output_gained_y}"
      )