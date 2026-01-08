import pandas as pd
import os
import sys


input_base_name= sys.argv[1:]

print("Starting make_x_element_repeatmasker_tsv.py")


for base_name in input_base_name:
    
    print(f'Opening {base_name}...') 
    
    repeatmasker_results_dir = f'/Shared/malkova_lab/Ivan/nanopore_sequencing/telomere_analysis/working_new_telomere_analysis/{base_name}/x_element/x_element_repeatmasker_results/'
    repeatmasker_results_dir_list = os.listdir(f'/Shared/malkova_lab/Ivan/nanopore_sequencing/telomere_analysis/working_new_telomere_analysis/{base_name}/x_element/x_element_repeatmasker_results/')

    all_results_files = [f'{repeatmasker_results_dir}/{f}' for f in repeatmasker_results_dir_list if f.startswith('guppy-') and f.endswith('results.ssv')]

    print(len(all_results_files))

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
            chr_end_of_file = corrected_repeatmasker_results_file.split("_whole_x_regions_repeatmasker_results_corrected.ssv")[0]
            chr_end_of_file = chr_end_of_file.split("detect-only_")[1]
            df_single_repeatmasker_results["chr_end"] = chr_end_of_file
            df_single_repeatmasker_results = df_single_repeatmasker_results.dropna(axis=1, how='all')
            df_all_repeatmakser_results = pd.concat([df_all_repeatmakser_results, df_single_repeatmasker_results])
        except:
            print(f'Error in {corrected_repeatmasker_results_file}')
            pass

    # Replace the "*" values with True and the "-" with False
    df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '*', 'sub_match'] = True
    df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '-', 'sub_match'] = False
    

    print(df_all_repeatmakser_results)
    
    df_all_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)

    df_all_repeatmakser_results.to_csv(f'{base_name}_x_element_repeatmasker.tsv' , sep= '\t')

    ################## Testing these cutoff values for a good match ##################
    df_good_repeatmakser_results = df_all_repeatmakser_results[df_all_repeatmakser_results['sub_match'] == False]
    
    df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['SW_score'] >= 500]
    df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['divergence_percent'] <= 2]
    
    output_filter_file = f'outputs/{base_name}_post_y_prime_probe.tsv'
    df_filter = pd.read_csv(output_filter_file, sep='\t') 
    df_filter = df_filter.dropna(subset=["repeat_length"])
    df_filter = df_filter[df_filter['Adapter_After_Telomere'] == True]

    good_reads = df_filter['read_id'].to_list()
    df_good_repeatmakser_results['good_read'] = df_good_repeatmakser_results['read_id'].apply(lambda x: x in good_reads)
    print(df_good_repeatmakser_results['good_read'].value_counts())
    
    df_good_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['good_read'] == True]

    df_good_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
    df_good_repeatmakser_results.to_csv(f'{base_name}_good_x_element_repeatmasker.tsv' , sep= '\t')
    
    ###############

    df_filter = df_filter[df_filter['delta_y_prime_sign'] == '+']
    
    reads_with_gain_of_y_prime = df_filter['read_id'].to_list()

    df_good_repeatmakser_results['gained_y'] = df_good_repeatmakser_results['read_id'].apply(lambda x: x in reads_with_gain_of_y_prime)
    print(df_all_repeatmakser_results['gained_y'].value_counts())
    
    df_good_gained_y_repeatmakser_results = df_good_repeatmakser_results[df_good_repeatmakser_results['gained_y'] == True]
    
    df_good_gained_y_repeatmakser_results.sort_values(by=['original_chr_end_anchor', 'read_id', 'match_start_on_read'], inplace=True)
    df_good_gained_y_repeatmakser_results.to_csv(f'{base_name}_x_element_repeatmasker_gained_y.tsv' , sep= '\t')

    ###############
    
    df_filter_0_to_more = df_filter[df_filter['reference_y_primes'] == 0]
    
    reads_0_to_more = df_filter_0_to_more['read_id'].to_list()
    
    df_good_gained_y_repeatmakser_results['more_than_0'] = df_good_gained_y_repeatmakser_results['read_id'].apply(lambda x: x in reads_0_to_more)
    
    print(df_good_gained_y_repeatmakser_results['more_than_0'].value_counts())
    
    df_more_than_0_y_repeatmakser_results = df_good_gained_y_repeatmakser_results[df_good_gained_y_repeatmakser_results['more_than_0'] == True]
    df_more_than_0_y_repeatmakser_results.to_csv(f'{base_name}_x_element_repeatmasker_gained_y_from_0.tsv' , sep= '\t')