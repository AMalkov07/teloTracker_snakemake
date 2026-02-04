import pandas as pd
import sys

def calculate_y_prime_change(read_id, chr_end):

    y_prime_ref_order = y_prime_order_dict[chr_end]
    recombination_events = {}
    
    if read_id not in df_y_repeatmasker['read_id'].values:
        y_prime_match_order = []
        better_repeatmasker_y_prime_count = 0    
    else:
        df_read = df_y_repeatmasker[df_y_repeatmasker['read_id'] == read_id].reset_index(drop=True)
        telomere_side = df_read.at[0, 'telomere_side']

        y_prime_match_order = df_read['y_prime_id_and_match'].unique()
        better_repeatmasker_y_prime_count = len(y_prime_match_order)

        # Determine what the read's first Y' ID is
        if telomere_side == 'beginning':
            # The first Y' will have the largest start position on the read (it will be the last entry by default ordering)
            y_prime_match_order = y_prime_match_order[::-1]
        else: # telomere_side == 'end':
        # The first Y' will have the smallest start position on the read (it will be the first entry by default ordering)
            y_prime_match_order = y_prime_match_order



    for i, y_prime_and_match_id in enumerate(y_prime_match_order):
        #print(y_prime_and_match_id)

        current_y_in_read = y_prime_and_match_id.split("-")[0]

        # Catch exceptions
        if y_prime_ref_order == None or i >= len(y_prime_ref_order):
            expected_y_in_ref = None
        else:
            expected_y_in_ref = y_prime_ref_order[i]

        #print(y_prime_ref_order)
        #print(f'from {expected_y_in_ref} to {current_y_in_read}')

        if current_y_in_read == expected_y_in_ref:
            pass
        else:
            recombination_events[f'read_y_num_{i+1}'] = f'from {expected_y_in_ref} to {current_y_in_read}'

    # Include losses of Y's
    if y_prime_ref_order != None:
        if len(y_prime_match_order) < len(y_prime_ref_order):

            for i in range(len(y_prime_match_order), len(y_prime_ref_order)):
                expected_y_in_ref = y_prime_ref_order[i]
                recombination_events[f'ref_y_num_{i+1}'] = f'from {expected_y_in_ref} to {None}'

    # Check if it is the day0 reference order (no recombination)
    if len(recombination_events) == 0:
        recombination_status = "No Change"
    # Check the earliest place of recombination (1st Y' or elsewhere)
    else:
        if 'read_y_num_1' in recombination_events.keys() or 'ref_y_num_1' in recombination_events.keys():
            recombination_status = "1st Y' Change"
        else:
            recombination_status = "Y' Recombination"

    #print(recombination_events, recombination_status)

    return recombination_events, recombination_status, better_repeatmasker_y_prime_count

print("Starting get_stats_of_recombination.py")

input_base_name = sys.argv[1:]

df_combined = pd.DataFrame()

for base_name in input_base_name:

    print(f'Opening {base_name}...')

    #continue

    df_file = f'results/{base_name}/{base_name}_good_end_y_repeatmasker.tsv'

    strain_id = f'{base_name.split("dorado_")[1]}'
    strain_id = f'{strain_id.split("_day")[0]}'

    df_y_repeatmasker = pd.read_csv(df_file, sep='\t')

    df_input = pd.read_csv(f'results/outputs/{base_name}_post_y_prime_probe.tsv', sep='\t')

    output_file = f'results/{base_name}/{base_name}_y_prime_recombination.tsv'

    df_input = df_input.dropna(subset=["repeat_length"])
    df_input = df_input[df_input['Adapter_After_Telomere'] == True]


    if '6991' == strain_id: # WildType
        y_prime_order_dict = {
            "chr1L": None,
            "chr1R": None,
            "chr2L": ("ID4",),
            "chr2R": None,
            "chr3L": None,
            "chr3R": None,
            "chr4L": None,
            "chr4R": ("ID2", "ID2", "ID2", "ID2", "ID2", "ID2", "ID2"),
            "chr5L": ("ID6",),
            "chr5R": ("ID1",),
            "chr6L": ("ID4",),
            "chr6R": None,
            "chr7L": None,
            "chr7R": ("ID5",),
            "chr8L": ("ID1",),
            "chr8R": ("ID1",),
            "chr9L": ("ID6",),
            "chr9R": None,
            "chr10L": ("ID6",),
            "chr10R": None,
            "chr11L": None,
            "chr11R": None,
            "chr12L": ("ID1",),
            "chr12R": ("ID1", "ID2", "ID2", "ID2", "ID2", "ID2", "ID2"),
            "chr13L": ("ID1",),
            "chr13R": None,
            "chr14L": ("ID5", "ID2", "ID3", "ID3", "ID3"),
            "chr14R": ("ID6",),
            "chr15L": None,
            "chr15R": ("ID2",),
            "chr16L": ("ID5",),
            "chr16R": ("ID1",)
        }

    elif '7172' == strain_id: # mph1
        y_prime_order_dict = {
            "chr1L": None,
            "chr1R": None,
            "chr2L": ("ID4",),
            "chr2R": None,
            "chr3L": None,
            "chr3R": None,
            "chr4L": None,
            "chr4R": ("ID2", "ID2", "ID2", "ID2", "ID2", "ID2", "ID2"),
            "chr5L": ("ID6",),
            "chr5R": ("ID1",),
            "chr6L": ("ID4",),
            "chr6R": None,
            "chr7L": None,
            "chr7R": ("ID5",),
            "chr8L": ("ID1",),
            "chr8R": ("ID1",),
            "chr9L": ("ID6",),
            "chr9R": None,
            "chr10L": ("ID6",),
            "chr10R": None,
            "chr11L": None,
            "chr11R": None,
            "chr12L": ("ID1",),
            "chr12R": ("ID1", "ID2", "ID2", "ID2", "ID2", "ID2", "ID2", "ID2"),
            "chr13L": ("ID1",),
            "chr13R": None,
            "chr14L": ("ID5", "ID2", "ID3"),
            "chr14R": ("ID6",),
            "chr15L": None,
            "chr15R": ("ID2",),
            "chr16L": ("ID5", "ID2", "ID2"),
            "chr16R": ("ID1",)
        }

    elif '7302' == strain_id: # mph1
        y_prime_order_dict = {
            "chr1L": None,
            "chr1R": None,
            "chr2L": ("ID4",),
            "chr2R": None,
            "chr3L": None,
            "chr3R": None,
            "chr4L": None,
            "chr4R": ("ID2", "ID2", "ID2", "ID2", "ID2", "ID2", "ID2"),
            "chr5L": ("ID6",),
            "chr5R": ("ID1",),
            "chr6L": ("ID4",),
            "chr6R": None,
            "chr7L": None,
            "chr7R": ("ID5",),
            "chr8L": ("ID1",),
            "chr8R": ("ID1",),
            "chr9L": ("ID6",),
            "chr9R": None,
            "chr10L": ("ID6",),
            "chr10R": None,
            "chr11L": None,
            "chr11R": None,
            "chr12L": ("ID1",),
            "chr12R": ("ID1", "ID2", "ID2", "ID2", "ID2", "ID2"),
            "chr13L": ("ID7", "ID2", "ID7", "ID2"),
            "chr13R": None,
            "chr14L": ("ID5", "ID2", "ID3", "ID3", "ID3"),
            "chr14R": ("ID6",),
            "chr15L": None,
            "chr15R": ("ID2",),
            "chr16L": ("ID5",),
            "chr16R": ("ID1",)
        }

    else:
        raise ValueError(f'Unknown strain_id: {strain_id}')

    df_y_repeatmasker['y_prime_id_and_color'] = df_y_repeatmasker['y_prime_group'].apply(lambda x: x.split('/')[2])
    df_y_repeatmasker['y_prime_id'] = df_y_repeatmasker['y_prime_id_and_color'].apply(lambda x: x.split('_')[0])
    df_y_repeatmasker['y_prime_id_and_match'] = df_y_repeatmasker['y_prime_id'].astype(str) + '-' + df_y_repeatmasker['match_id'].astype(str)

    all_reads = df_input['read_id'].unique()

    y_prime_recomb_events_dict = {}
    y_prime_recomb_status_dict = {}
    better_repeatmasker_y_prime_count_dict = {}
    for read_id in all_reads:
        chr_end = f"chr{df_input[df_input['read_id'] == read_id].iloc[0]['chr_end']}"
        y_prime_recombinations, recombination_status, better_repeatmasker_y_prime_count = calculate_y_prime_change(read_id, chr_end)
        y_prime_recomb_events_dict[read_id] = y_prime_recombinations
        y_prime_recomb_status_dict[read_id] = recombination_status
        better_repeatmasker_y_prime_count_dict[read_id] = better_repeatmasker_y_prime_count


    df_input['y_prime_recombination_events'] = df_input['read_id'].apply(lambda x: y_prime_recomb_events_dict[x])
    df_input['y_prime_recombination_status'] = df_input['read_id'].apply(lambda x: y_prime_recomb_status_dict[x])
    
    #######################
    df_input['better_repeatmasker_y_prime_count'] = df_input['read_id'].apply(lambda x: better_repeatmasker_y_prime_count_dict[x])
    df_input['better_repeatmasker_y_primes_relative_to_ref'] = df_input['better_repeatmasker_y_prime_count'] - df_input['reference_y_primes']
    #######################


    def calc_delta_group(better_repeatmasker_y_primes_relative_to_ref):
        if better_repeatmasker_y_primes_relative_to_ref == 0:
            delta_group = 'Same Number'
        elif better_repeatmasker_y_primes_relative_to_ref < 0:
            delta_group = 'Loss'
        elif better_repeatmasker_y_primes_relative_to_ref == 1:
            delta_group = 'Gain 1'
        else:
            delta_group = 'Gain Multiple'

        return delta_group

    df_input['y_prime_delta_group'] = df_input['better_repeatmasker_y_primes_relative_to_ref'].apply(lambda x: calc_delta_group(x))


#    for chr_end, group in df_input.groupby('chr_end'):
#        print(f"Chromosome End: {chr_end}, with {len(group)} reads")
#        print(group['y_prime_recombination_status'].value_counts())
#        print('\nPercentages:')
#        test = (group['y_prime_recombination_status'].value_counts(normalize=True)*100).to_string()
#        print(test)
#        percentage = (group['y_prime_recombination_status'].value_counts(normalize=True)*100).apply(lambda x: f"{x:.2f}%")
#        print(percentage)
#        print('\n\n')

    print(f'Toal Reads: {len(df_input)}\n\n')
    
    df_input.to_csv(output_file, sep='\t', index=False)
    
    continue
    
    #print(df_input['y_prime_recombination_events'])
    #print(df_input['y_prime_recombination_status'].value_counts())
    #percentage = (df_input['y_prime_recombination_status'].value_counts(normalize=True)*100).apply(lambda x: f"{x:.2f}%")
    #print(percentage)

    #print('\n')
    #print(df_input['y_prime_delta_group'].value_counts())
    #print(df_input['reference_y_prime_end_status'].value_counts().apply(lambda x: f"{x} reads").to_string())
    #print((df_input['reference_y_prime_end_status'].value_counts(normalize=True)*100).apply(lambda x: f"{x:.2f} %").to_string())
    df_output = df_input.groupby(by=['reference_y_prime_end_status'])
    #print('\n')

    #print(df_output['y_prime_delta_group'].value_counts().apply(lambda x: f"{x} reads").to_string())
    #print((df_output['y_prime_delta_group'].value_counts(normalize=True)*100).apply(lambda x: f"{x:.2f} %").to_string())
    df_output = df_input.groupby(by=['reference_y_prime_end_status', 'y_prime_delta_group'])
    #print('\n')
    
    
    
    #print(df_output['y_prime_recombination_status'].value_counts().apply(lambda x: f"{x} reads").to_string())
    #print((df_output['y_prime_recombination_status'].value_counts(normalize=True)*100).apply(lambda x: f"{x:.2f} %").to_string())

    df_combined = pd.concat([df_combined, df_input], ignore_index=True)
    print(f'Finished {base_name}.\n\n')

    continue
##############################
    print('\n\n\n')

    print(f'Toal Reads: {len(df_input)}')
    print('\n')
    print(df_input['reference_y_primes'].value_counts().apply(lambda x: f"{x} reads").to_string())
    print((df_input['reference_y_primes'].value_counts(normalize=True)*100).apply(lambda x: f"{x:.2f} %").to_string())
    df_output = df_input.groupby(by=['reference_y_primes'])
    print()

    print(df_output['y_prime_delta_group'].value_counts().apply(lambda x: f"{x} reads").to_string())
    print((df_output['y_prime_delta_group'].value_counts(normalize=True)*100).apply(lambda x: f"{x:.2f} %").to_string())
    df_output = df_input.groupby(by=['reference_y_primes', 'y_prime_delta_group'])
    print()
    
    print(df_output['y_prime_recombination_status'].value_counts().apply(lambda x: f"{x} reads").to_string())
    print((df_output['y_prime_recombination_status'].value_counts(normalize=True)*100).apply(lambda x: f"{x:.2f} %").to_string())

    df_y_rec = df_input[
        (df_input['reference_y_prime_end_status'] == 'Reference Y Primes = 2+') & 
        (df_input['y_prime_delta_group'] == 'Loss') &
        (df_input['y_prime_recombination_status'] == "Y' Recombination")
                ]['y_prime_recombination_events'].to_list()

    for value in df_y_rec:
        print(value)

    print('\n####################\n')

    df_pre_y_rec = df_input[
        (df_input['reference_y_prime_end_status'] == 'Reference Y Primes = 2+') & 
        (df_input['y_prime_delta_group'] == 'Loss') &
        (df_input['y_prime_recombination_status'] == "1st Y' Change")
                ]['y_prime_recombination_events'].to_list()
    
    for value in df_pre_y_rec:
        print(value)

#df_combined.to_csv('wt_day3_and_day4_y_prime_recombination.tsv', sep='\t', index=False)

df_final = pd.read_csv(f'results/{sys.argv[1]}/{sys.argv[1]}_y_prime_recombination.tsv', sep='\t')

total_reads = len(df_final)
print(f'Total Reads: {total_reads}\n\n')

# 1. y_prime_recombination_status
counts = df_final['y_prime_recombination_status'].value_counts()
percentages = (counts / total_reads * 100).apply(lambda x: f"{x:.2f}%")
print("y_prime_recombination_status Counts:")
print(counts)
print("\nPercentages (of total reads):")
print(percentages)
print("\n")

# 2. reference_y_prime_end_status within each recombination status
grouped1 = df_final.groupby('y_prime_recombination_status')['reference_y_prime_end_status'].value_counts()
percentages = (grouped1 / total_reads * 100).apply(lambda x: f"{x:.2f}%")
print("reference_y_prime_end_status by y_prime_recombination_status:")
print(grouped1)
print("\nPercentages (of total reads):")
print(percentages)
print("\n")

# 3. y_prime_delta_group within each recombination + end status combo
grouped2 = df_final.groupby(['y_prime_recombination_status', 'reference_y_prime_end_status'])['y_prime_delta_group'].value_counts()
percentages = (grouped2 / total_reads * 100).apply(lambda x: f"{x:.2f}%")
print("y_prime_delta_group by recombination status + reference end status:")
print(grouped2)
print("\nPercentages (of total reads):")
print(percentages)
print("\n")