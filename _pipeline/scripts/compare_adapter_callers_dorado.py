import pandas as pd
import os
import sys

# Usage: compare_adapter_callers_dorado.py <base_name> <anchor_set> <output_dir>
base_name  = sys.argv[1]
anchor_set = sys.argv[2]
output_dir = sys.argv[3]
input_fasta = sys.argv[4]

print("Starting compare_adapter_callers_dorado.py")

blast_dir    = os.path.join(output_dir, 'blast')
porechop_dir = os.path.join(output_dir, 'porechop')

porechop_results_f_name = os.path.join(porechop_dir, f'{base_name}_porechopped.tsv')
anchor_blast_f_name     = os.path.join(blast_dir, f'{base_name}_blasted_{anchor_set}.tsv')
best_anchor_telo_f_name = os.path.join(blast_dir, f'top_matches_{base_name}_blasted_{anchor_set}.tsv')

output_stats_f_name = os.path.join(output_dir, f'{base_name}_adapter_trimming_check.stats')
output_table_f_name = os.path.join(output_dir, f'{base_name}_adapter_trimming_check.tsv')

def check_for_adpt_or_tag(adapter_bool, tag_bool):
    return adapter_bool == True or tag_bool == True

def check_for_adpt_and_tag(adapter_bool, tag_bool):
    return adapter_bool == True and tag_bool == True

# Reads the number of lines in fasta file and divides by 2 to get the number of reads
def count_fastq_reads(fasta_file):
    with open(fasta_file, 'r') as f:
        num_lines = sum(1 for line in f)
    return num_lines // 2

total_reads = count_fastq_reads(input_fasta)

df_best_anchor = pd.read_csv(best_anchor_telo_f_name, sep='\t')
total_anchored_reads = len(df_best_anchor)

df_porechop = pd.read_csv(porechop_results_f_name, sep='\t')
df_best_anchor_porechop = pd.merge(df_porechop, df_best_anchor, on='read_id', how='inner')

# Dorado adapter columns set to False (not used in this pipeline version)
df_best_anchor_porechop['dorado_telomere_adapter']  = False
df_best_anchor_porechop['dorado_centromere_adapter'] = False

df_best_anchor_porechop['both_adapter_and_tag']  = df_best_anchor_porechop.apply(
    lambda r: check_for_adpt_and_tag(r['Adapter_After_Telomere'], r['Tag_After_Telomere']), axis=1)
df_best_anchor_porechop['either_adapter_or_tag'] = df_best_anchor_porechop.apply(
    lambda r: check_for_adpt_or_tag(r['Adapter_After_Telomere'], r['Tag_After_Telomere']), axis=1)

columns_to_include = ['read_id', 'Repeat_Type', 'Adapter_After_Telomere', 'Tag_After_Telomere',
                      'both_adapter_and_tag', 'either_adapter_or_tag', 'anchor_name',
                      'dorado_telomere_adapter', 'dorado_centromere_adapter']

df_final = df_best_anchor_porechop[columns_to_include]

print('\ndf_final\n')
print(df_final)
df_final.to_csv(output_table_f_name, sep='\t')

# Summary counts
df_porechop_adpt  = df_final[df_final['Adapter_After_Telomere'] == True]
df_tag            = df_final[df_final['Tag_After_Telomere'] == True]
df_tag_and_adpt   = df_tag[df_tag['Adapter_After_Telomere'] == True]
df_tag_only       = df_tag[df_tag['Adapter_After_Telomere'] == False]
df_either_adpt    = df_final[(df_final['dorado_telomere_adapter'] == True) | (df_final['Adapter_After_Telomere'] == True)]
df_adpt_only      = df_either_adpt[df_either_adpt['Tag_After_Telomere'] == False]

total_porechop_adapter_reads = len(df_porechop_adpt)
total_tag_reads              = len(df_tag)
total_tag_and_adpt_reads     = len(df_tag_and_adpt)
total_tag_only               = len(df_tag_only)
total_adpt_reads             = len(df_either_adpt)
total_adpt_only              = len(df_adpt_only)


lines = [
    f'Total_reads:\t{total_reads}\n',
    f'Total_anchored_reads:\t{total_anchored_reads}\n',
    f'Total_porechop_adapter_reads:\t{total_porechop_adapter_reads}\t'
    f'Total_tag_reads:\t{total_tag_reads}\t'
    f'Total_tag_only_reads:\t{total_tag_only}\t'
    f'Total_adapter_reads:\t{total_adpt_reads}\t'
    f'Total_adapter_only_reads:\t{total_adpt_only}\t'
    f'Total_both_tag_and_adapter_reads:\t{total_tag_and_adpt_reads}\t'
]

print(f'Writing adapter statistics to: {output_stats_f_name}')
with open(output_stats_f_name, 'w') as f:
    for line in lines:
        f.write(line)
