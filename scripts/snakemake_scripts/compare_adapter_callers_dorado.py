import pandas as pd
import os
import sys


# Base name of the file/run wanting to analyze
#base_name='dorado-6991_day0_PromethION_6991_day0_PromethION_TeloTag_yes_rejection-detect-only'
#base_name = f'{sys.argv[1]}'
anchor_set=f'{sys.argv[7]}'    # "new_all_best_anchors"  # f'{sys.argv[2]}' # Default = new_all_best_anchors

print("Starting compare_adapter_callers.py")

#dorado_summary_f_name=f'{base_name}/{base_name}_sequencing_summary.tsv'
dorado_summary_f_name=f'{sys.argv[1]}'
#porechop_results_f_name=f'{base_name}/{base_name}_porechopped.tsv'
porechop_results_f_name=f'{sys.argv[2]}'
#anchor_blast_f_name=f'{base_name}/{base_name}_blasted_{anchor_set}.tsv'
anchor_blast_f_name=f'{sys.argv[3]}'
#best_anchor_telo_f_name=f'{base_name}/top_matches_{base_name}_blasted_{anchor_set}.tsv'
best_anchor_telo_f_name=f'{sys.argv[4]}'


# Output Summary
output_stats_f_name=f'{sys.argv[5]}'
output_table_f_name=f'{sys.argv[6]}'


def check_for_adpt_or_tag(adapter_bool, tag_bool) -> bool:
    if adapter_bool == True or tag_bool == True:
        either_bool = True
    else:
        either_bool = False
    return either_bool

def check_for_adpt_and_tag(adapter_bool, tag_bool) -> bool:
    if adapter_bool == True and tag_bool == True:
        both_bool = True
    else:
        both_bool = False
    return both_bool



# Load TSV file into a Pandas DataFrame with custom column headers
df_all_anchor = pd.read_csv(anchor_blast_f_name, sep='\t')
#names=['read_id', 'total_read_length', 'read_bp_used_for_match', 'match_start_on_read', 'match_end_on_read', 'anchor_name', 'total_anchor_length', 'match_start_on_anchor', 'match_end_on_anchor', 'pident', 'bitscore', 'evalue']

# Filter for only sequencing reads with more than 2500 basepairs in read_bp_used_for_match
df_all_anchor = df_all_anchor[df_all_anchor['total_read_length'] > 2500]
# Filter for only reads with more than 2500 bp (half) matching the anchor
df_all_anchor = df_all_anchor[df_all_anchor['read_bp_used_for_match'] > 2500]
# Keep only one row for each unique entry in anchor_name (chromosome) per read_id (read)
df_all_anchor = df_all_anchor.drop_duplicates(subset=['read_id', 'anchor_name'], keep="first")
# Gets best match (top remaining) for the anchor
df_all_anchor = df_all_anchor.drop_duplicates(subset=['read_id'], keep="first")

# Get total number of reads with anchors
total_anchored_qscore_pass_reads = len(df_all_anchor)

df_dorado = pd.read_csv(dorado_summary_f_name, sep='\t')

# Get Total Numbers of Reads

total_passed_reads = len(df_dorado)

df_dorado_2500 = df_dorado[df_dorado['sequence_length_template'] >= 2500]
total_2500bp_passed_reads = len(df_dorado_2500)

df_dorado_5000 = df_dorado[df_dorado['sequence_length_template'] >= 5000]
total_5000bp_passed_reads = len(df_dorado_5000)

df_dorado_10000 = df_dorado[df_dorado['sequence_length_template'] >= 10000]
total_10000bp_passed_reads = len(df_dorado_10000)

df_dorado_20000 = df_dorado[df_dorado['sequence_length_template'] >= 20000]
total_20000bp_passed_reads = len(df_dorado_20000)


# Get Total Number of Anchored Reads

df_all_anchor_dorado = pd.merge(df_all_anchor, df_dorado, on='read_id', how='inner')

total_anchored_qscore_pass_reads = len(df_all_anchor_dorado)

df_porechop = pd.read_csv(porechop_results_f_name, sep = '\t')

df_all_anchor_dorado_porechop = pd.merge(df_all_anchor_dorado, df_porechop, on='read_id', how='inner')


columns_to_include = ['read_id',
                      'Adapter_After_Telomere', 'Tag_After_Telomere', 'Repeat_Type'] # Telomere_Sequence


df_all_anchor_dorado_porechop = df_all_anchor_dorado_porechop[columns_to_include]

df_best_anchor = pd.read_csv(best_anchor_telo_f_name, sep='\t')

total_pass_to_best_anchor_reads = len(df_best_anchor)


df_best_anchor_dorado_porechop = pd.merge(df_all_anchor_dorado_porechop, df_best_anchor, on='read_id', how='inner')

#######################
# Put False for both dorado adapters as I am not checking them
df_best_anchor_dorado_porechop['dorado_telomere_adapter'] = False
df_best_anchor_dorado_porechop['dorado_centromere_adapter'] = False
#######################
df_best_anchor_dorado_porechop["both_adapter_and_tag"] = df_best_anchor_dorado_porechop.apply(lambda row: check_for_adpt_and_tag(row["Adapter_After_Telomere"], row["Tag_After_Telomere"]), axis=1)
df_best_anchor_dorado_porechop["either_adapter_or_tag"] = df_best_anchor_dorado_porechop.apply(lambda row: check_for_adpt_or_tag(row["Adapter_After_Telomere"], row["Tag_After_Telomere"]), axis=1)

columns_to_include = ['read_id', 'Repeat_Type', 'Adapter_After_Telomere', 'Tag_After_Telomere',
                    'both_adapter_and_tag', 'either_adapter_or_tag', 'anchor_name', 'dorado_telomere_adapter', 'dorado_centromere_adapter'] # Telomere_Sequence

df_final = df_best_anchor_dorado_porechop[columns_to_include]

#print(df_final['Repeat_Type'].value_counts())
#print(df_final['anchor_name'].value_counts())

print('\ndf_final\n')
print(df_final)
df_final.to_csv(output_table_f_name, sep = '\t')



df_final_dorado_adpt = df_final[df_final['dorado_telomere_adapter'] == True]
total_dorado_adapter_reads = len(df_final_dorado_adpt)
#print('df_final_dorado_adpt')
#print(df_final_dorado_adpt['Repeat_Type'].value_counts())
#print(df_final_dorado_adpt['anchor_name'].value_counts())
#print('\n')

df_final_porechop_adpt = df_final[df_final['Adapter_After_Telomere'] == True]
total_porechop_adapter_reads = len(df_final_porechop_adpt)
#print('df_final_porechop_adpt')
#print(df_final_porechop_adpt['Repeat_Type'].value_counts())
#print(df_final_porechop_adpt['anchor_name'].value_counts())
#print('\n')

df_final_either_adpt = df_final[(df_final['dorado_telomere_adapter'] == True) | (df_final['Adapter_After_Telomere'] == True)]
total_either_adapter_reads = len(df_final_either_adpt)
print('df_final_either_adpt')
print(df_final_either_adpt['Repeat_Type'].value_counts())
print(df_final_either_adpt['anchor_name'].value_counts())
print('\n')

df_final_both_adpt = df_final_porechop_adpt[df_final_porechop_adpt['dorado_telomere_adapter'] == True]
total_both_adapter_reads = len(df_final_both_adpt)
#print('df_final_both_adpt')
#print(df_final_both_adpt['Repeat_Type'].value_counts())
#print(df_final_both_adpt['anchor_name'].value_counts())
#print('\n')

df_final_tag = df_final[df_final['Tag_After_Telomere'] == True]
total_tag_reads = len(df_final_tag)
#print('df_final_tag')
#print(df_final_tag['Repeat_Type'].value_counts())
#print(df_final_tag['anchor_name'].value_counts())
#print('\n')

df_final_tag_and_adpt = df_final_tag[df_final_tag['Adapter_After_Telomere'] == True]
total_tag_and_adpt_reads = len(df_final_tag_and_adpt)
#print('df_final_tag_and_adpt')
#print(df_final_tag_and_adpt['Repeat_Type'].value_counts())
#print(df_final_tag_and_adpt['anchor_name'].value_counts())
#print('\n')



df_final_tag_only = df_final_tag[df_final_tag['Adapter_After_Telomere'] == False]
total_tag_only = len(df_final_tag_only)

df_final_adpt = df_final_either_adpt
total_adpt_reads = len(df_final_adpt)

df_final_adpt_only = df_final_adpt[df_final_adpt['Tag_After_Telomere'] == False]
total_adpt_only = len(df_final_adpt_only)


lines=[
f'Total_passed_reads:\t{total_passed_reads}\tTotal_2.5kb_reads:\t{total_2500bp_passed_reads}_({round(total_2500bp_passed_reads/total_passed_reads*100,4)}%)\tTotal_5kb_reads:\t{total_5000bp_passed_reads}_({round(total_5000bp_passed_reads/total_passed_reads*100,4)}%)\tTotal_10kb_reads:\t{total_10000bp_passed_reads}_({round(total_10000bp_passed_reads/total_passed_reads*100,4)}%)\tTotal_20kb_reads:\t{total_20000bp_passed_reads}_({round(total_20000bp_passed_reads/total_passed_reads*100,4)}%)\n',
f'Total_qscore_pass_anchored_reads:\t{total_anchored_qscore_pass_reads}\t\t\tPercentage_of_5kb:\t{round(total_anchored_qscore_pass_reads/total_5000bp_passed_reads*100,4)}%\tPercentage_of_Total:\t{round(total_anchored_qscore_pass_reads/total_passed_reads*100,4)}%\n',
f'Total_qscore_pass_to_telomere_reads:\t{total_pass_to_best_anchor_reads}\t\t\tPercentage_of_Anchored:\t{round(total_pass_to_best_anchor_reads/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_pass_to_best_anchor_reads/total_passed_reads*100,4)}%\n',

f'Total_qscore_pass_dorado_adapter_reads:\t{total_dorado_adapter_reads}\tPercentage_of_to_Telomere:\t{round(total_dorado_adapter_reads/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_dorado_adapter_reads/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_dorado_adapter_reads/total_passed_reads*100,4)}%\n',
f'Total_qscore_pass_porechop_adapter_reads:\t{total_porechop_adapter_reads}\tPercentage_of_to_Telomere:\t{round(total_porechop_adapter_reads/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_porechop_adapter_reads/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_porechop_adapter_reads/total_passed_reads*100,4)}%\n',
f'Total_qscore_pass_either_called_adapter_reads:\t{total_either_adapter_reads}\tPercentage_of_to_Telomere:\t{round(total_either_adapter_reads/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_either_adapter_reads/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_either_adapter_reads/total_passed_reads*100,4)}%\n',
f'Total_qscore_pass_both_called_adapter_reads:\t{total_both_adapter_reads}\tPercentage_of_to_Telomere:\t{round(total_both_adapter_reads/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_both_adapter_reads/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_both_adapter_reads/total_passed_reads*100,4)}%\n',

f'Total_tag_reads:\t{total_tag_reads}\tPercentage_of_to_Telomere:\t{round(total_tag_reads/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_tag_reads/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_tag_reads/total_passed_reads*100,4)}%\n',
f'Total_tag_only_reads:\t{total_tag_only}\tPercentage_of_to_Telomere:\t{round(total_tag_only/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_tag_only/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_tag_only/total_passed_reads*100,4)}%\n',

f'Total_adapter_reads:\t{total_adpt_reads}\tPercentage_of_to_Telomere:\t{round(total_adpt_reads/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_adpt_reads/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_adpt_reads/total_passed_reads*100,4)}%\n',
f'Total_adapter_only_reads:\t{total_adpt_only}\tPercentage_of_to_Telomere:\t{round(total_adpt_only/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_adpt_only/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_adpt_only/total_passed_reads*100,4)}%\n',

f'Total_both_tag_and_adapter_reads:\t{total_tag_and_adpt_reads}\tPercentage_of_to_Telomere:\t{round(total_tag_and_adpt_reads/total_pass_to_best_anchor_reads*100,4)}%\tPercentage_of_Anchored:\t{round(total_tag_and_adpt_reads/total_anchored_qscore_pass_reads*100,4)}%\tPercentage_of_Total:\t{round(total_tag_and_adpt_reads/total_passed_reads*100,4)}%\n'
]

print(f'Writting adapter statistics to: {output_stats_f_name}')
with open (output_stats_f_name, 'w') as f:
    for line in lines:
        f.write(line)

        line=line.strip('\n')
        #print(line)


