import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os


input_base_name = sys.argv[1:]

print("Starting new_plot_switch_distances.py")


def calculate_index_of_switch(chr_end_matching_order_list, anchored_chr_end):

      conservate_switch_index = 0
      aggressive_switch_index = 0

      score = 0
      max_score = 0
      first_switch = True
      for i, matching_chr_end in enumerate(chr_end_matching_order_list):
            if matching_chr_end != anchored_chr_end:
                  if i == 0:
                        switch_to_chr_end = matching_chr_end
                  
                  score += 1
                  if score >= max_score:
                        aggressive_switch_index = i
            else:
                  if i == 0:
                        conservate_switch_index = 0
                        aggressive_switch_index = 0
                        switch_to_chr_end = anchored_chr_end
                        break
                  
                  else:
                        score -= 1

                        if first_switch == True:
                              first_switch = False
                              conservate_switch_index = i
      
      return switch_to_chr_end, conservate_switch_index, aggressive_switch_index

      


df_master_spacer_switches = pd.DataFrame()

for base_name in input_base_name:

      print(f'Starting {base_name}')

      continue

      read_id_list = []
      strain_list = []
      genotype_list = []
      day_list = []
      repeat_list = []
      
      strain = base_name.split('_')[1]
      
      if strain == '6991':
            genotype = 'wt'
      
      elif strain == '7172' or strain == '7302':
            genotype = 'mph1'
      
      day = base_name.split('_')[2]
      
      if 'repeat' not in base_name:
            repeat = 'original'
      else:
            repeat = base_name.split('_')[3]


      spacer_switchpoint_file = f'{base_name}_paired_good_spacer_repeatmasker.tsv'
      
      x_ends_switchpoint_file = f'{base_name}_good_x_element_ends_paired_repeatmasker.tsv'
      
      y_prime_repeatmasker_file = f'{base_name}_combined_repeatmasker.tsv'        
      
      recombination_info_file = f'{base_name}_y_prime_recombination.tsv'   


      print(f'Opening {spacer_switchpoint_file}...')
      
      df_paired_spacer_data = pd.read_csv(spacer_switchpoint_file, sep='\t')
      
      df_paired_x_ends_data = pd.read_csv(x_ends_switchpoint_file, sep='\t')
      
      df_y_prime_data = pd.read_csv(y_prime_repeatmasker_file, sep='\t')
      
      #############################   REMOVE THE PROBLEM ENDS   #############################
      #df_paired_spacer_data = df_paired_spacer_data[df_paired_spacer_data['original_chr_end_anchor'] != 'chr1L'].reset_index(drop=True)
      #df_paired_spacer_data = df_paired_spacer_data[df_paired_spacer_data['original_chr_end_anchor'] != 'chr1R'].reset_index(drop=True)
      
      reads_of_y_prime_switches = df_paired_spacer_data['read_id'].unique()
      print(len(reads_of_y_prime_switches))
      
      df_paired_x_ends_data = df_paired_x_ends_data[df_paired_x_ends_data['original_chr_end_anchor'] != 'chr1L'].reset_index(drop=True)
      df_paired_x_ends_data = df_paired_x_ends_data[df_paired_x_ends_data['original_chr_end_anchor'] != 'chr1R'].reset_index(drop=True)
      
      reads_of_y_prime_switches_x_ends = df_paired_x_ends_data['read_id'].unique()
      print(len(reads_of_y_prime_switches_x_ends))
      
      
      reads_of_y_prime_y_gain = df_y_prime_data['read_id'].unique()
      
      #df_y_prime_data = df_y_prime_data[df_y_prime_data['chr_end'] != 'chr1L'].reset_index(drop=True)
      #df_y_prime_data = df_y_prime_data[df_y_prime_data['chr_end'] != 'chr1R'].reset_index(drop=True)
      
      ###########################################################################################
      
      
      print(len(reads_of_y_prime_y_gain))
      
      
      ### For 0 --> More than 0 Y' events
      df_filter_for_0_to_many = df_y_prime_data[
                                                (df_y_prime_data['reference_y_primes'] == 0)
                                                & (df_y_prime_data['y_prime_id'] != None)]
      
      reads_of_y_prime_0_to_many = df_filter_for_0_to_many['read_id'].unique()
      #print(f"0 to Many Y' Reads: {len(reads_of_y_prime_0_to_many)}")
      
      # Must be in y_prime_0_to_many
      #reads_of_y_prime_switches = np.intersect1d(reads_of_y_prime_switches, reads_of_y_prime_0_to_many)



      ### For > 0 --> Reecombination Y' events
      df_filter_for_0_to_many = df_y_prime_data[
                                                (df_y_prime_data['reference_y_primes'] > 0)
                                                & (df_y_prime_data['y_prime_id'] != None)]
      
      reads_of_y_prime_0_to_many = df_filter_for_0_to_many['read_id'].unique()
      #print(f"0 to Many Y' Reads: {len(reads_of_y_prime_0_to_many)}")
      
      # Must be in y_prime_0_to_many
      #reads_of_y_prime_switches = np.intersect1d(reads_of_y_prime_switches, reads_of_y_prime_0_to_many)
      
      
      
      

      ### For any first Y' recombination event
      df_recombination_info = pd.read_csv(recombination_info_file, sep='\t')
      
      #print(df_recombination_info['y_prime_recombination_status'].value_counts())
      
      # Add filter of 1L and 1R
      df_filter_for_first_recombinations = df_recombination_info[
            (df_recombination_info['y_prime_recombination_status'] == "1st Y' Change")]
            #& (df_recombination_info['anchor_name'] != 'chr1L_anchor')
            #& (df_recombination_info['anchor_name'] != 'chr1R_anchor')]
      
      
      read_with_first_recombinations = df_filter_for_first_recombinations['read_id'].unique()
      
      print(f"Reads with 1st Y' recombinations: {len(read_with_first_recombinations)}")
      
      df_0_to_more_probe = df_filter_for_first_recombinations[
            df_filter_for_first_recombinations['y_prime_probe_count'] == 0]
      
      df_0_to_more_RM = df_filter_for_first_recombinations[
            df_filter_for_first_recombinations['better_repeatmasker_y_prime_count'] == 0]
      
      df_not_0_to_more_probe = df_filter_for_first_recombinations[
            df_filter_for_first_recombinations['y_prime_probe_count'] > 0]
      df_not_0_to_more_RM = df_filter_for_first_recombinations[
            df_filter_for_first_recombinations['better_repeatmasker_y_prime_count'] > 0]
      
      print(f"0 to More Reads: {len(df_0_to_more_probe)}")
      print(f"Not 0 to More Reads: {len(df_not_0_to_more_probe)}")
      
      print(f"0 to More Reads RM: {len(df_0_to_more_RM)}")
      print(f"Not 0 to More Reads RM: {len(df_not_0_to_more_RM)}")
      
      


      conservate_switch_distance_list = []
      aggressive_switch_distance_list = []
      spacer_switching_chr_end_pair_list = []
      
      has_x_end_list = []
      has_its_in_donor_list = []
      
      for read_id in reads_of_y_prime_switches:
            
            try:
                  df_read_paired_spacer_data = df_paired_spacer_data[df_paired_spacer_data['read_id'] == read_id].reset_index(drop=True)          
                  
                  anchored_chr_end = df_read_paired_spacer_data['original_chr_end_anchor'].iloc[0]
                  
                  strand = df_read_paired_spacer_data[df_read_paired_spacer_data['anchor_label'] == 'Is Anchor']['strand'].mode().iloc[0]
                  
                  df_read_paired_spacer_data['mathched_chr_end'] = df_read_paired_spacer_data['chr_end_tract'].apply(lambda x: x.split('_')[0])
                  
                  chr_end_matching_order_list = df_read_paired_spacer_data['mathched_chr_end'].to_list()
                  
                  
                  #print('\n\n\n')
                  #print(f'My end {anchored_chr_end}, {strand}')
                  #print(chr_end_matching_order_list)
                  
                  
                  if 'L' in anchored_chr_end:
                        if strand == '+':
                              # telomere_side is the top of the list
                              switch_to_chr_end, conservate_switch_index, aggressive_switch_index = calculate_index_of_switch(chr_end_matching_order_list, anchored_chr_end)
                        
                              if conservate_switch_index == 0:
                                    conservate_switch_distance = 0
                              else:
                                    conservate_min_switch_distance = df_read_paired_spacer_data['match_end_on_read'][conservate_switch_index]
                                    conservate_max_switch_distance = df_read_paired_spacer_data['match_start_on_read'][conservate_switch_index + 1]
                                    conservate_switch_distance = (conservate_max_switch_distance + conservate_min_switch_distance)/2 - df_read_paired_spacer_data['match_start_on_read'][0]
                              
                              if aggressive_switch_index == 0:
                                    aggressive_switch_distance = 0
                              else:
                                    aggressive_min_switch_distance = df_read_paired_spacer_data['match_end_on_read'][aggressive_switch_index]
                                    aggressive_max_switch_distance = df_read_paired_spacer_data['match_start_on_read'][aggressive_switch_index + 1]
                                    aggressive_switch_distance = (aggressive_max_switch_distance + aggressive_min_switch_distance)/2 - df_read_paired_spacer_data['match_start_on_read'][0]    
                                    
                        else: # strand = 'C'
                              # telomere_side is the bottom of the list
                              chr_end_matching_order_list.reverse()
                              switch_to_chr_end, conservate_switch_index, aggressive_switch_index = calculate_index_of_switch(chr_end_matching_order_list, anchored_chr_end)
                              
                              
                              list_index_length = len(chr_end_matching_order_list) - 1
                        
                              if conservate_switch_index == 0:
                                    conservate_switch_distance = 0
                              else:
                                    conservate_min_switch_distance = df_read_paired_spacer_data['match_start_on_read'][list_index_length - (conservate_switch_index)]
                                    conservate_max_switch_distance = df_read_paired_spacer_data['match_end_on_read'][list_index_length - (conservate_switch_index + 1)]
                                    conservate_switch_distance = df_read_paired_spacer_data['match_start_on_read'][list_index_length] - (conservate_max_switch_distance + conservate_min_switch_distance)/2

                              if aggressive_switch_index == 0:
                                    aggressive_switch_distance = 0
                              else:
                                    aggressive_min_switch_distance = df_read_paired_spacer_data['match_start_on_read'][list_index_length - conservate_switch_index]
                                    aggressive_max_switch_distance = df_read_paired_spacer_data['match_end_on_read'][list_index_length - (conservate_switch_index + 1)]
                                    aggressive_switch_distance = df_read_paired_spacer_data['match_start_on_read'][list_index_length] - (aggressive_max_switch_distance + aggressive_min_switch_distance)/2 
                                    
                  else: #R in anchored_chr_end
                        if strand == '+':
                              # telomere_side is the bottom of the list
                              chr_end_matching_order_list.reverse()
                              switch_to_chr_end, conservate_switch_index, aggressive_switch_index = calculate_index_of_switch(chr_end_matching_order_list, anchored_chr_end)
                              
                              
                              list_index_length = len(chr_end_matching_order_list) - 1
                        
                              if conservate_switch_index == 0:
                                    conservate_switch_distance = 0
                              else:
                                    conservate_min_switch_distance = df_read_paired_spacer_data['match_start_on_read'][list_index_length - (conservate_switch_index)]
                                    conservate_max_switch_distance = df_read_paired_spacer_data['match_end_on_read'][list_index_length - (conservate_switch_index + 1)]
                                    conservate_switch_distance = df_read_paired_spacer_data['match_start_on_read'][list_index_length] - (conservate_max_switch_distance + conservate_min_switch_distance)/2

                              if aggressive_switch_index == 0:
                                    aggressive_switch_distance = 0
                              else:
                                    aggressive_min_switch_distance = df_read_paired_spacer_data['match_start_on_read'][list_index_length - conservate_switch_index]
                                    aggressive_max_switch_distance = df_read_paired_spacer_data['match_end_on_read'][list_index_length - (conservate_switch_index + 1)]
                                    aggressive_switch_distance = df_read_paired_spacer_data['match_start_on_read'][list_index_length] - (aggressive_max_switch_distance + aggressive_min_switch_distance)/2 
                                    
                        else: # strand = 'C'
                              # telomere_side is the top of the list
                              switch_to_chr_end, conservate_switch_index, aggressive_switch_index = calculate_index_of_switch(chr_end_matching_order_list, anchored_chr_end)
                              
                              if conservate_switch_index == 0:
                                    conservate_switch_distance = 0
                              else:
                                    conservate_min_switch_distance = df_read_paired_spacer_data['match_end_on_read'][conservate_switch_index]
                                    conservate_max_switch_distance = df_read_paired_spacer_data['match_start_on_read'][conservate_switch_index + 1]
                                    conservate_switch_distance = (conservate_max_switch_distance + conservate_min_switch_distance)/2 - df_read_paired_spacer_data['match_start_on_read'][0]
                              
                              if aggressive_switch_index == 0:
                                    aggressive_switch_distance = 0
                              else:
                                    aggressive_min_switch_distance = df_read_paired_spacer_data['match_end_on_read'][aggressive_switch_index]
                                    aggressive_max_switch_distance = df_read_paired_spacer_data['match_start_on_read'][aggressive_switch_index + 1]
                                    aggressive_switch_distance = (aggressive_max_switch_distance + aggressive_min_switch_distance)/2 - df_read_paired_spacer_data['match_start_on_read'][0]
                        
                  #print(anchored_chr_end, strand, switch_to_chr_end, conservate_switch_distance, aggressive_switch_distance)
                  
                  
                  if anchored_chr_end == switch_to_chr_end:
                        spacer_switching_chr_end_pair = 'No Switch'
                  else:
                        spacer_switching_chr_end_pair = f'Anchored {anchored_chr_end} switch to {switch_to_chr_end}'
                              
                              
                  if read_id in reads_of_y_prime_switches_x_ends:
                                    
                        df_read_paired_x_ends_data = df_paired_x_ends_data[df_paired_x_ends_data['read_id'] == read_id].reset_index(drop=True)
                        
                        x_element_end = df_read_paired_x_ends_data['x_element_ends'].iloc[0].split('_')[0]
                        
                        #print(x_element_end)
                        
                        if x_element_end == anchored_chr_end:
                              
                              x_element_read_outcome = "No Switch in X element"
                              
                        else:
                              x_element_read_outcome = "Switch in X element"
                  
                        #df_paired_x_ends_data['leftover_on_read']
                        
                  else:
                        x_element_read_outcome = 'No Results'

                  #print(read_id)

                  ends_with_its = ['chr4R', 'chr5L', 'chr5R', 'chr6L', 'chr8L', 'chr12R', 'chr13L', 'chr14L', 'chr14R', 'chr16L']
                  
                  try:
                        donor_y_prime_chr_ends = df_y_prime_data[df_y_prime_data['read_id'] == read_id]['y_prime_id'].iloc[0]
                        
                        for end_with_its in ends_with_its:
                              if end_with_its in donor_y_prime_chr_ends:
                                    its_in_donor = True
                                    break
                              else:
                                    its_in_donor = False
                  except IndexError:
                        its_in_donor = 'No Result'
                        print(read_id)
                        print(read_id in reads_of_y_prime_y_gain)
                        print(its_in_donor)
                  
                  read_id_list.append(read_id)
                  strain_list.append(strain)
                  genotype_list.append(genotype)
                  day_list.append(day)
                  repeat_list.append(repeat)
                  conservate_switch_distance_list.append(conservate_switch_distance)
                  aggressive_switch_distance_list.append(aggressive_switch_distance)
                  spacer_switching_chr_end_pair_list.append(spacer_switching_chr_end_pair)
                  has_x_end_list.append(x_element_read_outcome)
                  has_its_in_donor_list.append(its_in_donor)
            except:
                  print(f'Error with read_id {read_id} in strain {strain}, day {day}, repeat {repeat}')
                  # print the error message
                  print(sys.exc_info()[0])
                  
                  read_id_list.append(read_id)
                  strain_list.append(None)
                  genotype_list.append(None)
                  day_list.append(None)
                  repeat_list.append(None)
                  conservate_switch_distance_list.append(None)
                  aggressive_switch_distance_list.append(None)
                  spacer_switching_chr_end_pair_list.append(None)
                  has_x_end_list.append(None)
                  has_its_in_donor_list.append(None)
      
      strain_spacer_switches_dict = {
      'read_id': read_id_list,
      'strain': strain_list,
      'genotype': genotype_list,
      'day': day_list,
      'repeat': repeat_list,
      'conservate_switch_distance': conservate_switch_distance_list,
      'aggressive_switch_distance': aggressive_switch_distance_list,
      'anchor_and_spacer_switch_to_chr_end_pair': spacer_switching_chr_end_pair_list,
      'has_x_end': has_x_end_list,
      'has_its_in_donor': has_its_in_donor_list
      }
      
      df_strain_spacer_switches = pd.DataFrame(strain_spacer_switches_dict, index=None)
            
      df_master_spacer_switches = pd.concat([df_master_spacer_switches, df_strain_spacer_switches])



########################################################################

#df_master_spacer_switches.to_csv('all_spacer_switches-all_changes.tsv', index=False, sep='\t')
df_master_spacer_switches = pd.read_csv('new_output_combined_y_prime_recombination_with_spacer_switches.tsv', sep='\t')

df_master_spacer_switches['group'] = df_master_spacer_switches['source_file'].apply(
      lambda x: x.split('_')[1])


df_stats_wt = df_master_spacer_switches[df_master_spacer_switches['group'] == '6991']

total_wt_reads = len(df_stats_wt)

# Print the total number of reads and percentage of each recombination status
print(f'Total WT Reads: {total_wt_reads}')
print(df_stats_wt['y_prime_recombination_status'].value_counts())
print(df_stats_wt['y_prime_recombination_status'].value_counts() / total_wt_reads * 100)
print('\n')

df_y_prime_recombination = df_stats_wt[df_stats_wt['y_prime_recombination_status'] == "Y' Recombination"]
print(f'Total WT Reads with Y\' Recombination: {len(df_y_prime_recombination)}')
print(df_y_prime_recombination['y_prime_delta_group'].value_counts())
print(df_y_prime_recombination['y_prime_delta_group'].value_counts() / total_wt_reads * 100)
df_y_prime_recombination_loss = df_y_prime_recombination[
      df_y_prime_recombination['y_prime_delta_group'] == "Loss"]
print(f'Total WT Reads with Y\' Recombination Loss: {len(df_y_prime_recombination_loss)}')
print(df_y_prime_recombination_loss['chr_end'].value_counts())
print(df_y_prime_recombination_loss['chr_end'].value_counts() / total_wt_reads * 100)
print('\n')

df_pre_y_prime_0_start = df_stats_wt[(df_stats_wt['y_prime_recombination_status'] == "1st Y' Change")
                             & (df_stats_wt['reference_y_primes'] == 0)]
print(f'Total WT Reads with 1st Y\' Change starting 0: {len(df_pre_y_prime_0_start)}')
print(df_pre_y_prime_0_start['y_prime_delta_group'].value_counts())
print(df_pre_y_prime_0_start['y_prime_delta_group'].value_counts() / total_wt_reads * 100)
print('\n')

df_pre_y_prime_1m_start = df_stats_wt[(df_stats_wt['y_prime_recombination_status'] == "1st Y' Change")
                             & (df_stats_wt['reference_y_primes'] > 0)]
print(f'Total WT Reads with 1st Y\' Change starting 1+: {len(df_pre_y_prime_1m_start)}')
print(df_pre_y_prime_1m_start['y_prime_delta_group'].value_counts())
print(df_pre_y_prime_1m_start['y_prime_delta_group'].value_counts() / total_wt_reads * 100)
print('\n')


# Check if group is '7172' or '7302'
df_stats_mph1 = df_master_spacer_switches[df_master_spacer_switches['group'].isin(['7172', '7302'])]


total_mph1_reads = len(df_stats_mph1)

# Print the total number of reads and percentage of each recombination status
print(f'Total mph1 Reads: {total_mph1_reads}')
print(df_stats_mph1['y_prime_recombination_status'].value_counts())
print(df_stats_mph1['y_prime_recombination_status'].value_counts() / total_mph1_reads * 100)
print('\n')

df_y_prime_recombination = df_stats_mph1[df_stats_mph1['y_prime_recombination_status'] == "Y' Recombination"]
print(f'Total mph1 Reads with Y\' Recombination: {len(df_y_prime_recombination)}')
print(df_y_prime_recombination['y_prime_delta_group'].value_counts())
print(df_y_prime_recombination['y_prime_delta_group'].value_counts() / total_mph1_reads * 100)
df_y_prime_recombination_loss = df_y_prime_recombination[
      df_y_prime_recombination['y_prime_delta_group'] == "Loss"]
print(f'Total mph1 Reads with Y\' Recombination Loss: {len(df_y_prime_recombination_loss)}')
print(df_y_prime_recombination_loss['chr_end'].value_counts())
print(df_y_prime_recombination_loss['chr_end'].value_counts() / total_mph1_reads * 100)

print('\n')

df_pre_y_prime_0_start = df_stats_mph1[(df_stats_mph1['y_prime_recombination_status'] == "1st Y' Change")
                             & (df_stats_mph1['reference_y_primes'] == 0)]
print(f'Total mph1 Reads with 1st Y\' Change starting 0: {len(df_pre_y_prime_0_start)}')
print(df_pre_y_prime_0_start['y_prime_delta_group'].value_counts())
print(df_pre_y_prime_0_start['y_prime_delta_group'].value_counts() / total_mph1_reads * 100)
print('\n')

df_pre_y_prime_1m_start = df_stats_mph1[(df_stats_mph1['y_prime_recombination_status'] == "1st Y' Change")
                             & (df_stats_mph1['reference_y_primes'] > 0)]
print(f'Total mph1 Reads with 1st Y\' Change starting 1+: {len(df_pre_y_prime_1m_start)}')
print(df_pre_y_prime_1m_start['y_prime_delta_group'].value_counts())
print(df_pre_y_prime_1m_start['y_prime_delta_group'].value_counts() / total_mph1_reads * 100)
print('\n')




print(df_master_spacer_switches)










print(f'\nAll Ends\n')
print(df_master_spacer_switches['genotype'].value_counts())
print('\n')
#print(df_master_spacer_switches['anchor_and_spacer_switch_to_chr_end_pair'].value_counts())
#print('\n\n\n')

df_wt = df_master_spacer_switches[df_master_spacer_switches['genotype'] == 'wt']
#print(f'\nWT Ends\n')
#print(df_wt['genotype'].value_counts())
#print('\n')
#print(df_wt['anchor_and_spacer_switch_to_chr_end_pair'].value_counts())
#print('\n')
#print(df_wt[df_wt['anchor_and_spacer_switch_to_chr_end_pair'] == 'Anchored chr2R switch to chr13L']['aggressive_switch_distance'].value_counts())
#print('\n')
#print(df_wt['has_x_end'].value_counts())
#print('\n')
#print('\n\n\n')

df_mph = df_master_spacer_switches[df_master_spacer_switches['genotype'] == 'mph1']
#print(f'\nMPH Ends\n')
#print(df_mph['genotype'].value_counts())
#print('\n')
#print(df_mph['anchor_and_spacer_switch_to_chr_end_pair'].value_counts())
#print('\n')
#print(df_mph[df_mph['anchor_and_spacer_switch_to_chr_end_pair'] == 'Anchored chr2R switch to chr13L']['aggressive_switch_distance'].value_counts())
#print('\n')
#print(df_mph['has_x_end'].value_counts())
#print('\n\n\n')




wt_num = len(df_wt)

# Switch in Spacer
wt_in_spacer_switch_num = len(df_wt[df_wt['anchor_and_spacer_switch_to_chr_end_pair'] != 'No Switch'])
df_wt_not_switch_in_spacer_num = df_wt[df_wt['anchor_and_spacer_switch_to_chr_end_pair'] == 'No Switch']

# Switch in X element
wt_in_x_num = len(df_wt_not_switch_in_spacer_num[df_wt_not_switch_in_spacer_num['has_x_end'] == 'Switch in X element'])
df_wt_no_switch_found = df_wt_not_switch_in_spacer_num[df_wt_not_switch_in_spacer_num['has_x_end'] != 'Switch in X element']

# No Switch Found
wt_no_switch_found_num = len(df_wt_no_switch_found)

### ITS in Donor
wt_no_switch_found_with_its_in_donor_num = len(df_wt_no_switch_found[df_wt_no_switch_found['has_its_in_donor'] == True])

### Unknown
wt_no_switch_found_no_its_in_donor_num = len(df_wt_no_switch_found[df_wt_no_switch_found['has_its_in_donor'] != True])

print(f"WT Switch 1st Y' Events: {wt_num}")
print(f'WT Switch in Spacer Number: {wt_in_spacer_switch_num}, {((wt_in_spacer_switch_num/wt_num) * 100):.2f}%')
print(f'WT Switch in X Element Number: {wt_in_x_num}, {((wt_in_x_num/wt_num) * 100):.2f}%')
print(f'WT No Switch Found with Donor has ITS: {wt_no_switch_found_with_its_in_donor_num}, {((wt_no_switch_found_with_its_in_donor_num/wt_num) * 100):.2f}%')
print(f'WT No Switch Found - Unknown: {wt_no_switch_found_no_its_in_donor_num}, {((wt_no_switch_found_no_its_in_donor_num/wt_num) * 100):.2f}%')
print(f'WT Total No Switch Found: {wt_no_switch_found_num}, {((wt_no_switch_found_num/wt_num) * 100):.2f}%')

df_wt = df_wt[df_wt['anchor_and_spacer_switch_to_chr_end_pair'] != 'No Switch']


print('\n')


mph_num = len(df_mph)

# Switch in Spacer
mph_in_spacer_switch_num = len(df_mph[df_mph['anchor_and_spacer_switch_to_chr_end_pair'] != 'No Switch'])
df_mph_not_switch_in_spacer_num = df_mph[df_mph['anchor_and_spacer_switch_to_chr_end_pair'] == 'No Switch']

# Switch in X element
mph_in_x_num = len(df_mph_not_switch_in_spacer_num[df_mph_not_switch_in_spacer_num['has_x_end'] == 'Switch in X element'])
df_mph_no_switch_found = df_mph_not_switch_in_spacer_num[df_mph_not_switch_in_spacer_num['has_x_end'] != 'Switch in X element']

# No Switch Found
mph_no_switch_found_num = len(df_mph_no_switch_found)

### ITS in Donor
mph_no_switch_found_with_its_in_donor_num = len(df_mph_no_switch_found[df_mph_no_switch_found['has_its_in_donor'] == True])

### Unknown
mph_no_switch_found_no_its_in_donor_num = len(df_mph_no_switch_found[df_mph_no_switch_found['has_its_in_donor'] != True])

print(f"mph Switch 1st Y' Events: {mph_num}")
print(f'mph Switch in Spacer Number: {mph_in_spacer_switch_num}, {((mph_in_spacer_switch_num/mph_num) * 100):.2f}%')
print(f'mph Switch in X Element Number: {mph_in_x_num}, {((mph_in_x_num/mph_num) * 100):.2f}%')
print(f'mph No Switch Found with Donor has ITS: {mph_no_switch_found_with_its_in_donor_num}, {((mph_no_switch_found_with_its_in_donor_num/mph_num) * 100):.2f}%')
print(f'mph No Switch Found - Unknown: {mph_no_switch_found_no_its_in_donor_num}, {((mph_no_switch_found_no_its_in_donor_num/mph_num) * 100):.2f}%')
print(f'mph Total No Switch Found: {mph_no_switch_found_num}, {((mph_no_switch_found_num/mph_num) * 100):.2f}%')









"""
#stop

for switch_type in ['conservate_switch_distance', 'aggressive_switch_distance']:


      df_master_spacer_switches['switch_distance_rev'] = df_master_spacer_switches[switch_type] * -1

      df_no_0 = df_master_spacer_switches.loc[df_master_spacer_switches[switch_type] > 0]
      #print(df_no_0['genotype'].value_counts())

      figure_directory = f'spacer_outputs/{switch_type}_switch_location_figures'
      os.makedirs(figure_directory, exist_ok=True)

      sns.histplot(data=df_master_spacer_switches, x='switch_distance_rev', hue='genotype', cumulative=True, stat="density", common_norm=False, kde=True, binwidth=50)

      save_file_name = f'all_switches_cumulative_density_curve_of_{switch_type}_switch_location.png'
      plt.title(f"{save_file_name.removesuffix('.png')}")
      print(f'Graphing {figure_directory}/{save_file_name}...')
      plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
      plt.show()


      sns.histplot(data=df_master_spacer_switches, x='switch_distance_rev', hue='genotype', cumulative=False, stat="percent", common_norm=False, kde=True, binwidth=250)

      save_file_name = f'all_switches_frequency_curve_of_{switch_type}_switch_location.png'
      plt.title(f"{save_file_name.removesuffix('.png')}")
      print(f'Graphing {figure_directory}/{save_file_name}...')
      plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
      plt.show()
      



      sns.histplot(data=df_no_0, x='switch_distance_rev', hue='genotype', cumulative=True, stat="density", common_norm=False, kde=True, binwidth=50)

      save_file_name = f'all_switches_cumulative_density_curve_no_0_of_{switch_type}_switch_location.png'
      plt.title(f"{save_file_name.removesuffix('.png')}")
      print(f'Graphing {figure_directory}/{save_file_name}...')
      plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
      plt.show()


      sns.histplot(data=df_no_0, x='switch_distance_rev', hue='genotype', cumulative=False, stat="percent", common_norm=False, kde=True, binwidth=250)

      save_file_name = f'all_switches_frequency_curve_no_0_of_{switch_type}_switch_location.png'
      plt.title(f"{save_file_name.removesuffix('.png')}")
      print(f'Graphing {figure_directory}/{save_file_name}...')
      plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
      plt.show()


      #sns.barplot(data=)


"""

