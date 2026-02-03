import pandas as pd
import pysam
import sys
import os

base_name=f'{sys.argv[1]}'

anchor_set=f'{sys.argv[2]}'    # "new_all_best_anchors"  # f'{sys.argv[2]}' # Default = new_all_best_anchors

input=f'{base_name}_blasted_{anchor_set}'
all_reads_f_name=f'results/{base_name}/{base_name}.fasta'
input_f_name_best=f'results/{base_name}/top_matches_{input}.tsv'

individual_outputs_dir=f'results/{base_name}/chr_anchor_included_individual_files/'

os.system(f'mkdir -p {individual_outputs_dir}')

print("Starting split_and_label_all_reads_include_anchor.py")

chr_list = ['chr10R',
'chr11R',
'chr12R',
'chr13R',
'chr14R',
'chr15R',
'chr16R',
'chr17R',
'chr1R',
'chr2R',
'chr3R',
'chr4R',
'chr5R',
'chr6R',
'chr7R',
'chr8R',
'chr9R',
'chr10L',
'chr11L',
'chr12L',
'chr13L',
'chr14L',
'chr15L',
'chr16L',
'chr17L',
'chr1L',
'chr2L',
'chr3L',
'chr4L',
'chr5L',
'chr6L',
'chr7L',
'chr8L',
'chr9L'
]

# Load TSV file into a Pandas DataFrame with custom column headers
df_best = pd.read_csv(input_f_name_best, sep='\t')

fasta_index = pysam.FastaFile(all_reads_f_name)

for chr in chr_list:
    
    anchor = f'{chr}_anchor'

    output_read_f_name=f'{individual_outputs_dir}/{input}_{anchor}_reads.fasta'
    print(f'Starting to make {output_read_f_name}...')


    df_new = df_best[df_best['anchor_name'] == anchor]
    df_new = df_new.sort_values(['read_length_past_anchor'], ascending=False)
    df_new.to_csv(f'{individual_outputs_dir}/{input}_{anchor}_reads.tsv', sep="\t",index=False)

    read_list_for_chr_end = zip(df_new['read_id'], df_new['repeat_type'])

    print(f'Grabbing reads in {anchor}...')
    # Open the output file for writing
    split_reads_list = []

    for read_name, repeat_type in read_list_for_chr_end:
        #print(read_name)
        sequence = fasta_index.fetch(read_name)
        sequence = sequence.strip('\n')
        row_index = df_best.loc[df_best['read_id'] == read_name].index[0]
        #print(row_index)
        qstart_index = df_best.at[row_index,'match_start_on_read'] - 1
        qend_index = df_best.at[row_index,'match_end_on_read'] - 1
        wanted_section_of_read = df_best.at[row_index,'wanted_section_of_read']
        if wanted_section_of_read == 'before_match_start_on_read':
            chopped_seq = sequence[:qend_index]
        elif wanted_section_of_read == 'after_match_end_on_read':
            chopped_seq = sequence[qstart_index:]
        else:
            continue
        #print(sequence)
        #print(len(chopped_seq))
        chopped_seq = chopped_seq.split('\n')[0]
        
        if len(chopped_seq) == 0:
            continue
        
        header = f'{read_name} {repeat_type} {chr}'
        #print(header)
        split_reads_list.append(f'>{header}\n{chopped_seq}\n')


    with open(f'{output_read_f_name}', "w") as f_out:
        f_out.writelines(split_reads_list)
