import pandas as pd
import os
import sys


input_base_name = sys.argv[1:]

print("Starting make_y_prime_repeatmasker_tsv.py")


for file_name in input_base_name:
    
    print(f'Opening {file_name}...') 
    
    repeatmasker_results_dir = f'results/{file_name}/read_repeatmasker_results/'
    repeatmasker_results_dir_list = os.listdir(f'results/{file_name}/read_repeatmasker_results/')

    all_results_files = [f'{repeatmasker_results_dir}/{f}' for f in repeatmasker_results_dir_list if f.startswith('dorado') and f.endswith('results.ssv')]

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
            chr_end_of_file = corrected_repeatmasker_results_file.split("_repeatmasker_results_corrected.ssv")[0]
            chr_end_of_file = chr_end_of_file.split("rejection_")[1]
            df_single_repeatmasker_results["chr_end"] = chr_end_of_file
            df_all_repeatmakser_results = pd.concat([df_all_repeatmakser_results, df_single_repeatmasker_results])
        except:
            print(f'Error in {corrected_repeatmasker_results_file}')
            pass

    # Replace the "*" values with True and the "-" with False
    df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '*', 'sub_match'] = True
    df_all_repeatmakser_results.loc[df_all_repeatmakser_results['sub_match'] == '-', 'sub_match'] = False
    #df_all_repeatmakser_results['sub_match'] = df_all_repeatmakser_results['sub_match'].apply(lambda x: x == '*')

    print(df_all_repeatmakser_results)

    df_all_repeatmakser_results.to_csv(f'results/{file_name}_repeatmasker.tsv' , sep= '\t')

    ###############

    output_filter_file = f'results/outputs/{file_name}_post_y_prime_probe.tsv'
    df_filter = pd.read_csv(output_filter_file, sep='\t') 
    
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
    
    df_all_repeatmakser_results['telomere_side'] = df_all_repeatmakser_results.apply(lambda row: determine_telomere_side(row["chr_end"], row["Repeat_Type"]), axis=1)
    
    reads_with_good_ends = df_filter['read_id'].to_list()
    
    print(len(reads_with_good_ends))
    
    ######
    
    df_all_repeatmakser_results['good_ends'] = df_all_repeatmakser_results['read_id'].apply(lambda x: x in reads_with_good_ends)
    
    print(df_all_repeatmakser_results['good_ends'].value_counts())
    
    df_good_end_y_repeatmakser_results = df_all_repeatmakser_results[df_all_repeatmakser_results['good_ends'] == True]
    df_good_end_y_repeatmakser_results.to_csv(f'results/{file_name}_good_end_y_repeatmasker.tsv' , sep= '\t')
    
    
    ######
    
    df_filter = df_filter[df_filter['delta_y_prime_sign'] == '+']
    
    reads_with_gain_of_y_prime = df_filter['read_id'].to_list()
    
    print(len(reads_with_gain_of_y_prime))

    df_all_repeatmakser_results['gained_y'] = df_all_repeatmakser_results['read_id'].apply(lambda x: x in reads_with_gain_of_y_prime)
    
    print(df_all_repeatmakser_results['gained_y'].value_counts())
    
    df_gained_y_repeatmakser_results = df_all_repeatmakser_results[df_all_repeatmakser_results['gained_y'] == True]
    df_gained_y_repeatmakser_results.to_csv(f'results/{file_name}_gained_y_repeatmasker.tsv' , sep= '\t')
