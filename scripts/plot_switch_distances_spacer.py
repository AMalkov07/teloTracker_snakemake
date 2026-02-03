import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys


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

      strain_list = []
      genotype_list = []
      day_list = []
      repeat_list = []
      
      strain = base_name.split('_')[0]
      strain = strain.split('-')[1]
      
      if strain == '6991':
            genotype = 'wt'
      
      elif strain == '7172' or strain == '7302':
            genotype = 'mph1'
      
      day = base_name.split('_')[1]
      
      if 'repeat' not in base_name.split('_')[2]:
            repeat = 'original'
      else:
            repeat = base_name.split('_')[2]


      spacer_switchpoint_file = f'{base_name}_paired_good_gained_spacer_repeatmasker.tsv'

      print(f'Opening {spacer_switchpoint_file}...')
      
      df_paired_spacer_data = pd.read_csv(spacer_switchpoint_file, sep='\t')
      
      #############################   REMOVE THE PROBLEM ENDS   #############################
      df_paired_spacer_data = df_paired_spacer_data[df_paired_spacer_data['anchor_label'] != 'chr1L'].reset_index(drop=True)
      df_paired_spacer_data = df_paired_spacer_data[df_paired_spacer_data['anchor_label'] != 'chr1R'].reset_index(drop=True)
      
      reads_of_0_to_1_switches = df_paired_spacer_data['read_id'].unique()


      conservate_switch_distance_list = []
      aggressive_switch_distance_list = []
      spacer_switching_chr_end_pair_list = []
      
      for read_id in reads_of_0_to_1_switches:
            
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
                  
            print(anchored_chr_end, strand, switch_to_chr_end, conservate_switch_distance, aggressive_switch_distance)
            
            if anchored_chr_end == switch_to_chr_end:
                  spacer_switching_chr_end_pair = 'No Switch'
            else:
                  spacer_switching_chr_end_pair = f'Anchored {anchored_chr_end} switch to {switch_to_chr_end}'
                        
            strain_list.append(strain)
            genotype_list.append(genotype)
            day_list.append(day)
            repeat_list.append(repeat)
            conservate_switch_distance_list.append(conservate_switch_distance)
            aggressive_switch_distance_list.append(aggressive_switch_distance)
            spacer_switching_chr_end_pair_list.append(spacer_switching_chr_end_pair)
      
      strain_spacer_switches_dict = {
      'strain': strain_list,
      'genotype': genotype_list,
      'day': day_list,
      'repeat': repeat_list,
      'conservate_switch_distance': conservate_switch_distance_list,
      'aggressive_switch_distance': aggressive_switch_distance_list,
      'anchor_and_spacer_switch_to_chr_end_pair': spacer_switching_chr_end_pair_list
      }
      
      df_strain_spacer_switches = pd.DataFrame(strain_spacer_switches_dict, index=None)
            
      df_master_spacer_switches = pd.concat([df_master_spacer_switches, df_strain_spacer_switches])
            
                      



#################################### Showing only day4 ! ####################################

df_master_spacer_switches = df_master_spacer_switches[df_master_spacer_switches['day'] == 'day4']


print(f'\nAll Ends\n')
print(df_master_spacer_switches['genotype'].value_counts())
print('\n')
print(df_master_spacer_switches['anchor_and_spacer_switch_to_chr_end_pair'].value_counts())
print('\n\n\n')

df_wt = df_master_spacer_switches[df_master_spacer_switches['genotype'] == 'wt']
print(f'\nWT Ends\n')
print(df_wt['genotype'].value_counts())
print('\n')
print(df_wt['anchor_and_spacer_switch_to_chr_end_pair'].value_counts())
print('\n')
print(df_wt[df_wt['anchor_and_spacer_switch_to_chr_end_pair'] == 'Anchored chr2R switch to chr13L']['aggressive_switch_distance'].value_counts())
print('\n\n\n')

df_mph = df_master_spacer_switches[df_master_spacer_switches['genotype'] == 'mph1']
print(f'\nMPH Ends\n')
print(df_mph['genotype'].value_counts())
print('\n')
print(df_mph['anchor_and_spacer_switch_to_chr_end_pair'].value_counts())
print('\n')
print(df_mph[df_mph['anchor_and_spacer_switch_to_chr_end_pair'] == 'Anchored chr2R switch to chr13L']['aggressive_switch_distance'].value_counts())
print('\n\n\n')


for switch_type in ['conservate_switch_distance', 'aggressive_switch_distance']:

      #print(df_master_spacer_switches['genotype'].value_counts())

      df_master_spacer_switches['switch_distance_rev'] = df_master_spacer_switches[switch_type] * -1

      df_no_0 = df_master_spacer_switches.loc[df_master_spacer_switches[switch_type] > 0]
      #print(df_no_0['genotype'].value_counts())

      figure_directory = f'spacer_outputs/{switch_type}_switch_location_figures'

      sns.histplot(data=df_master_spacer_switches, x='switch_distance_rev', hue='genotype', cumulative=True, stat="density", common_norm=False, kde=True, binwidth=50)

      save_file_name = f'cumulative_density_curve_of_{switch_type}_switch_location.png'
      plt.title(f"{save_file_name.removesuffix('.png')}")
      print(f'Graphing {figure_directory}/{save_file_name}...')
      #plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
      plt.show()


      sns.histplot(data=df_master_spacer_switches, x='switch_distance_rev', hue='genotype', cumulative=False, stat="percent", common_norm=False, kde=True, binwidth=250)

      save_file_name = f'frequency_curve_of_{switch_type}_switch_location.png'
      plt.title(f"{save_file_name.removesuffix('.png')}")
      print(f'Graphing {figure_directory}/{save_file_name}...')
      #plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
      plt.show()


      sns.histplot(data=df_no_0, x='switch_distance_rev', hue='genotype', cumulative=True, stat="density", common_norm=False, kde=True, binwidth=50)

      save_file_name = f'cumulative_density_curve_no_0_of_{switch_type}_switch_location.png'
      plt.title(f"{save_file_name.removesuffix('.png')}")
      print(f'Graphing {figure_directory}/{save_file_name}...')
      #plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
      plt.show()


      sns.histplot(data=df_no_0, x='switch_distance_rev', hue='genotype', cumulative=False, stat="percent", common_norm=False, kde=True, binwidth=250)

      save_file_name = f'frequency_curve_no_0_of_{switch_type}_switch_location.png'
      plt.title(f"{save_file_name.removesuffix('.png')}")
      print(f'Graphing {figure_directory}/{save_file_name}...')
      #plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
      plt.show()
