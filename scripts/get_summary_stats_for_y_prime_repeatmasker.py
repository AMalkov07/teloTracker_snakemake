import pandas as pd
import os
import sys


print("Starting get_summary_stats_for_y_prime_repeatmasker.py")

file_name = sys.argv[1]

anchor_set = sys.argv[2]

strain_id = sys.argv[3]



print(f'Opening {file_name}...')

repeatmasker_results_dir = f'results/{file_name}/read_repeatmasker_results/'
repeatmasker_results_dir_list = os.listdir(f'results/{file_name}/read_repeatmasker_results/')

all_results_files = [f'{repeatmasker_results_dir}/{f}' for f in repeatmasker_results_dir_list if f.startswith('dorado') and f.endswith('results.ssv')]

print(f'Found {len(all_results_files)} Repeatmasker results files in {repeatmasker_results_dir}')

#strain_id = file_name.split("-")[1].split("_")[0]
print(f'Strain ID: {strain_id}')

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
        chr_end_of_file = chr_end_of_file.split("detect-only_")[1]
        df_single_repeatmasker_results["chr_end"] = chr_end_of_file
        df_all_repeatmakser_results = pd.concat([df_all_repeatmakser_results, df_single_repeatmasker_results])
    except:
        print(f'Error in {corrected_repeatmasker_results_file}')
        pass

# Replace the "*" values with True and the "-" with False
df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '*', 'sub_match'] = True
df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '-', 'sub_match'] = False

#print(df_all_repeatmakser_results)

# Split the Y' Group statuses to individual columns
df_all_repeatmakser_results['y_prime_size'] = df_all_repeatmakser_results['y_prime_group'].apply(lambda x: x.split('/')[0])
df_all_repeatmakser_results['y_prime_tandem_status'] = df_all_repeatmakser_results['y_prime_group'].apply(lambda x: x.split('/')[1])
df_all_repeatmakser_results['y_prime_color'] = df_all_repeatmakser_results['y_prime_group'].apply(lambda x: x.split('/')[2])

# Annotate if end has an ITS
def its_annotation(chr_end):
    # Select the chromosome ends with an ITS
    if chr_end in ['chr4R', 'chr5L', 'chr5R', 'chr6L', 'chr8L', 'chr12R', 'chr13L', 'chr14L', 'chr14R', 'chr16L']:
        if strain_id == '6991':
            if chr_end == 'chr4R':
                its_legnth = 381
            if chr_end == 'chr5L':
                its_legnth = 23
            if chr_end == 'chr5R':
                its_legnth = 15
            if chr_end == 'chr6L':
                its_legnth = 114
            if chr_end == 'chr8L':
                its_legnth = 157
            if chr_end == 'chr12R':
                its_legnth = 125
            if chr_end == 'chr13L':
                its_legnth = 51
            if chr_end == 'chr14L':
                its_legnth = 263
            if chr_end == 'chr14R':
                its_legnth = 72
            if chr_end == 'chr16L':
                its_legnth = 11
        elif strain_id == '7172':
            if chr_end == 'chr4R':
                its_legnth = 381
            if chr_end == 'chr5L':
                its_legnth = 25
            if chr_end == 'chr5R':
                its_legnth = 15
            if chr_end == 'chr6L':
                its_legnth = 111
            if chr_end == 'chr8L':
                its_legnth = 5
            if chr_end == 'chr12R':
                its_legnth = 125
            if chr_end == 'chr13L':
                its_legnth = 48
            if chr_end == 'chr14L':
                its_legnth = 263
            if chr_end == 'chr14R':
                its_legnth = 102
            if chr_end == 'chr16L':
                its_legnth = 6      
        elif strain_id == '7302':
            if chr_end == 'chr4R':
                its_legnth = 381
            if chr_end == 'chr5L':
                its_legnth = 25
            if chr_end == 'chr5R':
                its_legnth = 15
            if chr_end == 'chr6L':
                its_legnth = 111
            if chr_end == 'chr8L':
                its_legnth = 5
            if chr_end == 'chr12R':
                its_legnth = 125
            if chr_end == 'chr13L':
                its_legnth = 48
            if chr_end == 'chr14L':
                its_legnth = 263
            if chr_end == 'chr14R':
                its_legnth = 72
            if chr_end == 'chr16L':
                its_legnth = 6               
    else:
        its_legnth = 0
    return its_legnth
    
df_all_repeatmakser_results['ITS_length'] = df_all_repeatmakser_results['chr_end'].apply(lambda x: its_annotation(x))

df_all_repeatmakser_results.to_csv(f'results/{file_name}/{file_name}_repeatmasker.tsv', sep='\t')


# Get the previous master tsv results into a dataframe
filter_file = f'outputs/{file_name}_post_y_prime_probe.tsv'
df_filter = pd.read_csv(filter_file, sep='\t') 

df_filter = df_filter.dropna(subset=["repeat_length"])
df_filter = df_filter[df_filter['Adapter_After_Telomere'] == True]

##########################   Filter and analyze all the reads that gained Y's   ##########################

df_combined_repeatmakser_results = df_all_repeatmakser_results.merge(df_filter, left_on='read_id', right_on='read_id', how='inner', suffixes=(None, '_filter'))

df_combined_repeatmakser_results.to_csv(f'results/{file_name}/{file_name}_combined_repeatmasker.tsv', sep='\t')


##########################   Filter and analyze below only the reads that gained Y's   ##########################

df_gained_y_combined_repeatmakser_results = df_combined_repeatmakser_results[df_combined_repeatmakser_results['delta_y_prime_sign'] == '+']


###########   Compares the Y' counts for repeatmasker of full Y's and the Y' probe    ###########

# Get Y' counts for each read by repeatmasker
rm_y_prime_count_dict = df_gained_y_combined_repeatmakser_results['read_id'].value_counts().to_dict()
df_gained_y_combined_repeatmakser_results['repeatmasker_y_prime_count'] = df_gained_y_combined_repeatmakser_results['read_id'].apply(lambda x: rm_y_prime_count_dict[x] if x in rm_y_prime_count_dict else 0)

# Check if it matches the probe
df_gained_y_combined_repeatmakser_results['y_prime_count_delta'] = df_gained_y_combined_repeatmakser_results['repeatmasker_y_prime_count'] - df_gained_y_combined_repeatmakser_results['y_prime_probe_count']

# Remove sub matches that are not a full Y' (i.e. Y' is not gained)
df_gained_y_combined_repeatmakser_results = df_gained_y_combined_repeatmakser_results.dropna(subset=["sub_match"])

# Get Y's relative to the reference for repeatmasker's counts
df_gained_y_combined_repeatmakser_results['repeatmasker_delta_y_primes_to_ref'] = df_gained_y_combined_repeatmakser_results['repeatmasker_y_prime_count'] - df_gained_y_combined_repeatmakser_results['reference_y_primes']

##################################################################################################


###########   Gets the "multiplers" for the different Y' categories for normalization of counts    ###########

# Get multiplier for each unique Y' (Y' ID) repeated
def get_relative_y_prime_id(y_prime_id):
    # Split the Y' ID into the a list of the chromosome end with the Y' ID
    all_chr_in_id = y_prime_id.split("_")[2]
    chr_groups = all_chr_in_id.split(";")
    
    # Go through the chromosome ends and count the number of Y' ID's in each end
    # Add that to the running total of true Y's in the reference for the ID      
    y_prime_id_ref_occurances = 0
    for chr_group in chr_groups:
        y_prime_id_ref_occurances += len(chr_group.split(","))
        
    return y_prime_id_ref_occurances  
    
df_gained_y_combined_repeatmakser_results["y_prime_id_multiplier"] = df_gained_y_combined_repeatmakser_results['y_prime_id'].apply(lambda x: get_relative_y_prime_id(x))  

# Get one entry of every unique Y' ID to get the relative group info
df_relative_group_info = df_gained_y_combined_repeatmakser_results.drop_duplicates(subset=['y_prime_id'])

# Y prime size counters
y_prime_long_ref_occurances = 0
y_prime_short_ref_occurances = 0
# Y prime Tandeem status counters
y_prime_tandem_ref_occurances = 0
y_prime_solo_ref_occurances = 0
# Y prime color dictionary counter
y_prime_color_count_dict = {}

# Get the multiplier of each relative group info
for index, row in df_relative_group_info.iterrows():
    if row['y_prime_size'] == "Long":
        y_prime_long_ref_occurances += row['y_prime_id_multiplier']
    else: # Short
        y_prime_short_ref_occurances += row['y_prime_id_multiplier']
        
    if row['y_prime_tandem_status'] == "Tandem":
        y_prime_tandem_ref_occurances += row['y_prime_id_multiplier']
    else: # Solo
        y_prime_solo_ref_occurances += row['y_prime_id_multiplier']
        
    if row['y_prime_color'] in y_prime_color_count_dict.keys():
        y_prime_color_count_dict[row['y_prime_color']] += row['y_prime_id_multiplier']
    else:
        y_prime_color_count_dict[row['y_prime_color']] = row['y_prime_id_multiplier']

df_gained_y_combined_repeatmakser_results["y_prime_size_multiplier"] = df_gained_y_combined_repeatmakser_results['y_prime_size'].apply(
    lambda x: y_prime_long_ref_occurances if x == "Long" else y_prime_short_ref_occurances)

df_gained_y_combined_repeatmakser_results["y_prime_tandem_status_multiplier"] = df_gained_y_combined_repeatmakser_results['y_prime_tandem_status'].apply(
    lambda x: y_prime_tandem_ref_occurances if x == "Tandem" else y_prime_solo_ref_occurances)

df_gained_y_combined_repeatmakser_results["y_prime_color_multiplier"] = df_gained_y_combined_repeatmakser_results['y_prime_color'].apply(
    lambda x: y_prime_color_count_dict[x])

#############################################################################################################


# Saves to .tsv and prints the gained Y' dataframe with all the new columns 

df_gained_y_combined_repeatmakser_results.to_csv(f'results/{file_name}/{file_name}_gained_y_repeatmasker.tsv', sep='\t')

###########   Get Total Reads sequenced   ###########

# Get the previous master tsv results into a dataframe
all_reads_file = f'{file_name}/{file_name}_adapter_trimming_check.stats'
with open(all_reads_file, 'r') as f:
    first_line = f.readline()
    sequencing_run_total_reads = int(first_line.split('\t')[1])

# Get the best/top read anchor location for all reads tsv into a dataframe
achored_file = f'{file_name}/top_matches_{file_name}_blasted_{anchor_set}.tsv'
df_anchored = pd.read_csv(achored_file, sep='\t') 

total_anchored_reads = len(df_anchored["anchor_name"])


# Get total Telomere Reads sequenced
df_combined_repeatmakser_results_reads = df_combined_repeatmakser_results.drop_duplicates(subset=['read_id'])
telomere_read_counts = len(df_filter)
telomere_read_with_y_prime_counts = len(df_combined_repeatmakser_results_reads)
telomere_read_percentage_of_total = f"{((telomere_read_counts / sequencing_run_total_reads) * 100):.2f}"

# Get total Y' Gained Reads sequenced
df_gained_y_combined_repeatmakser_results_reads = df_gained_y_combined_repeatmakser_results.drop_duplicates(subset=['read_id'])
y_prime_gained_read_counts = len(df_gained_y_combined_repeatmakser_results_reads)
#y_prime_gained_percentage_of_total = f"{((y_prime_gained_read_counts / sequencing_run_total_reads) * 100):.2f}"
y_prime_gained_percentage_of_telomere = f"{((y_prime_gained_read_counts / telomere_read_counts) * 100):.2f}"

# Get the total Y's Gained
total_y_primes_gained = df_gained_y_combined_repeatmakser_results_reads['repeatmasker_delta_y_primes_to_ref'].sum()
y_primes_per_gaining_read = f"{(total_y_primes_gained / y_prime_gained_read_counts):.2f}"

dataframe_info = {'Total Reads':[sequencing_run_total_reads],
                    'Anchored Reads':[total_anchored_reads],
                    'Telomere Reads':[telomere_read_counts],
                    'Telomere Reads % of Total':[telomere_read_percentage_of_total], 
                    "Y' Gained Reads":[y_prime_gained_read_counts],
                    "Y' Gained Reads % of Telomere":[y_prime_gained_percentage_of_telomere],
                    "Total Y's Gained":[total_y_primes_gained],
                    "Y's Gained per Gaining Read":[y_primes_per_gaining_read]}


# Make dataframe for printing
#df_sequencing_run_info = sequencing_run_total_reads.to_frame(name='Total Reads in Sequencing Run')
# Reset the index and drop the index column
#df_sequencing_run_info.reset_index(drop=True, inplace=True)
#df_sequencing_run_info['Telomere Reads'] = y_prime_gained_read_counts
#df_sequencing_run_info['Telomere Reads Percentage of Total Reads'] = telomere_read_percentage_of_total
#df_sequencing_run_info["Y' Gained Reads"] = y_prime_gained_read_counts
#df_sequencing_run_info["Y' Gained Reads Percentage of Total Reads"] = y_prime_gained_percentage_of_total
#df_sequencing_run_info["Y' Gained Reads Percentage of Telomere Reads"] = y_prime_gained_percentage_of_telomere

df_sequencing_run_info = pd.DataFrame(dataframe_info)


###########   Gets the number of reads info for all ends, the reads that gained Y's, the frequency, and the percentage   ###########

# Get the total number of anchored reads for each chromosome end
df_anchored['Anchor_Name'] = df_anchored['anchor_name'].apply(lambda x: x.split('_')[0])
anchored_reads_for_chr_end = df_anchored["Anchor_Name"].value_counts()


# Get total reads/counts telomere reads for each chromosome end
df_filter['Chromosome_end'] = df_filter['chr_end'].apply(lambda x: f'chr{x}')
telomere_reads_for_chr_end = df_filter['Chromosome_end'].value_counts()

# Get the frequency of the read reaching the telomere sequence for each chromosome end
frequency_of_reaching_telomere_for_chr_end = (telomere_reads_for_chr_end.div(anchored_reads_for_chr_end, fill_value=0) * 100).round(2)

# Get total reads/counts sequenced with a gain of Y' for each chromosome end
y_prime_gained_reads_for_chr_end = df_gained_y_combined_repeatmakser_results_reads['chr_end'].value_counts()

# Get the frequency of Y' gained by the total reads sequenced for each chromosome end
frequency_of_gaining_y_primes_for_chr_end = (y_prime_gained_reads_for_chr_end.div(telomere_reads_for_chr_end, fill_value=0) * 100).round(2)

# Percentage the chromosome end makes up the reads with a gain of Y'
percentage_y_prime_gained_reads_for_chr_end = ((y_prime_gained_reads_for_chr_end / y_prime_gained_reads_for_chr_end.sum()) * 100).round(2)

# Make dataframe for printing
df_reads_sequenced_info = anchored_reads_for_chr_end.to_frame(name='Anchored Reads')
df_reads_sequenced_info['Telomere Reads'] = telomere_reads_for_chr_end
df_reads_sequenced_info['Freq. of Reaching Telomere'] = frequency_of_reaching_telomere_for_chr_end
df_reads_sequenced_info['Gain of Y Primes Reads'] = y_prime_gained_reads_for_chr_end
df_reads_sequenced_info['Freq. of Gaining Y Primes'] = frequency_of_gaining_y_primes_for_chr_end
df_reads_sequenced_info['Percentage of Gained Y Prime Reads'] = percentage_y_prime_gained_reads_for_chr_end  

####################################################################################################################################


###########   Gets the number of reads info for all 0, 1, 2+ ends   ###########

# Get total reads/counts sequenced with a gain of Y' for each reference Y' Group (0, 1, and 2+)
y_prime_gained_reads_ref_y_prime = df_gained_y_combined_repeatmakser_results_reads['reference_y_prime_end_status'].value_counts()

# Percentage the chromosome end makes up the reads with a gain of Y'
percentage_y_prime_gained_reads_ref_y_prime = ((y_prime_gained_reads_ref_y_prime / y_prime_gained_reads_ref_y_prime.sum()) * 100).round(2)

# Make dataframe for printing
df_ref_y_prime_group_info = y_prime_gained_reads_ref_y_prime.to_frame(name='Total Reference Y Group Counts')
df_ref_y_prime_group_info['Percentage Reference Y Group'] = percentage_y_prime_gained_reads_ref_y_prime

################################################################################


###########   Gets the delta Y's of the number of identified Y's (repeatmakser - probe) in the Y' gained reads   ###########

# Gets delta of the counts of Y' gained by repeatmakser - probe 
df_delta_y_prime_counts = df_gained_y_combined_repeatmakser_results_reads["y_prime_count_delta"].value_counts()

# Get the percentage of Y' gained deltas
percentage_y_prime_gained_reads_for_chr_end = ((df_delta_y_prime_counts / df_delta_y_prime_counts.sum()) * 100).round(2)

# Make dataframe for printing
df_delta_y_prime_info = df_delta_y_prime_counts.to_frame(name="Y' Delta Counts")
df_delta_y_prime_info["Frequency of Y' Deltas"] = percentage_y_prime_gained_reads_for_chr_end 

#############################################################################################################################


###########   Gets the Y' Info Statistics in the Y' gained reads   ###########

# Get counts of each unique Y' (Y' ID) for reads that gained Y's
### y_prime_id_dict = df_gained_y_combined_repeatmakser_results['y_prime_id'].value_counts().to_dict()

# Unique Y' ID
unique_y_prime_counts = df_gained_y_combined_repeatmakser_results['y_prime_id'].value_counts()

unique_y_prime_percentages = ((unique_y_prime_counts / unique_y_prime_counts.sum()) * 100).round(2)

# Make dataframe for printing
df_unique_y_prime_info = unique_y_prime_counts.to_frame(name="Unique Y' Counts")
df_unique_y_prime_info["Frequency of Y' Deltas"] = percentage_y_prime_gained_reads_for_chr_end

# Normalize the counts by the multiplier (Number of occurances of unique Y' ID in the reference genome)
unique_y_prime_multiplier_dict = dict(zip(df_gained_y_combined_repeatmakser_results['y_prime_id'], df_gained_y_combined_repeatmakser_results['y_prime_id_multiplier']))
unique_y_prime_normalized_counts = pd.Series({k: (unique_y_prime_counts.to_dict()[k] * unique_y_prime_multiplier_dict[k]) for k in unique_y_prime_multiplier_dict.keys()}, 
                                                name = 'Normalized y_prime_id')
#print(unique_y_prime_normalized_counts)

unique_y_prime_normalized_percentages = ((unique_y_prime_normalized_counts / unique_y_prime_normalized_counts.sum()) * 100).round(2).to_frame(name='Frequency')    
#print(unique_y_prime_normalized_percentages)

# Color
color_y_prime_counts = df_gained_y_combined_repeatmakser_results['y_prime_color'].value_counts()
color_y_prime_percentages = ((color_y_prime_counts / color_y_prime_counts.sum()) * 100).round(2).to_frame(name='Frequency')

# Normalize the counts by the multiplier (Number of occurances of unique Y' ID in the reference genome)
color_y_prime_multiplier_dict = dict(zip(df_gained_y_combined_repeatmakser_results['y_prime_color'], df_gained_y_combined_repeatmakser_results['y_prime_color_multiplier']))
color_y_prime_normalized_counts = unique_y_prime_counts.map(color_y_prime_multiplier_dict) * unique_y_prime_counts
#print(unique_y_prime_normalized_counts)



unique_y_prime_normalized_percentages = ((unique_y_prime_normalized_counts / unique_y_prime_normalized_counts.sum()) * 100).round(2).to_frame(name='Frequency')    
#print(unique_y_prime_normalized_percentages)   


# Size
size_y_prime_counts = df_gained_y_combined_repeatmakser_results['y_prime_size'].value_counts()
size_y_prime_percentages = ((size_y_prime_counts / size_y_prime_counts.sum()) * 100).round(2).to_frame(name='Frequency')  
#print(size_y_prime_percentages)


# Tandem Status
tandem_y_prime_counts = df_gained_y_combined_repeatmakser_results['y_prime_tandem_status'].value_counts()
tandem_y_prime_percentages = ((tandem_y_prime_counts / tandem_y_prime_counts.sum()) * 100).round(2).to_frame(name='Frequency')
#print(tandem_y_prime_percentages)    
    

with open(f'repeatmasker_summary_stats/{file_name}_repeatmasker_stats.txt', "w") as f:
    f.write(f'\t\t\tAnalysis for {file_name}\n\n')
    
    f.write(f'All ends\n')
    f.write(f'-----------------------------------------------------------------------------------------------------------------\n\n')
    
    f.write(f"Reads in Sequencing Run:\n")
    f.write(f"{df_sequencing_run_info.to_string(index=False)}\n\n")
    
    f.write(f'Counts of Telomere Reads by Chromosome End\n')
    f.write(f"{df_reads_sequenced_info.to_string(index=True)}\n\n")
    
    f.write(f"Counts of Reference Y' Groups in Reads with Y' gain:\n")
    f.write(f"{df_ref_y_prime_group_info}\n\n")      

    f.write(f"Delta Y's gained by (repeatmakser - probe) in Reads with Y' gain:\n")
    f.write(f"{df_delta_y_prime_info}\n\n\n")          
    
    
    f.write(f"% of each Unique Y' Gained in Reads with Y' gain:\n")
    f.write(f"{unique_y_prime_percentages}\n\n")
    
    f.write(f"% of each color of Y' Gained in Reads with Y' gain:\n")
    f.write(f"{color_y_prime_percentages}\n\n")        
    
    f.write(f"% of Solo/Tandem Y' Gained in Reads with Y' gain:\n")
    f.write(f"{tandem_y_prime_percentages}\n\n")    

    f.write(f"% of Long/Short Y' Gained in Reads with Y' gain:\n")
    f.write(f"{size_y_prime_percentages}\n\n")

print('Finshed!') 


