import pandas as pd
import pysam
import sys
import os

# Usage: split_and_label_all_reads_include_anchor.py <base_name> <anchor_set> <output_dir>
base_name  = sys.argv[1]
anchor_set = sys.argv[2]
output_dir = sys.argv[3]

blast_input       = f'{base_name}_blasted_{anchor_set}'
all_reads_f_name  = os.path.join(output_dir, f'{base_name}.fasta')
input_f_name_best = os.path.join(output_dir, 'blast', f'top_matches_{blast_input}.tsv')

individual_outputs_dir = os.path.join(output_dir, 'blast', 'chr_anchor_reads')
os.makedirs(individual_outputs_dir, exist_ok=True)

print("Starting split_and_label_all_reads_include_anchor.py")

chr_list = (
    [f'chr{n}R' for n in range(1, 17)] +
    [f'chr{n}L' for n in range(1, 17)]
)

df_best     = pd.read_csv(input_f_name_best, sep='\t')
fasta_index = pysam.FastaFile(all_reads_f_name)

for chr in chr_list:
    anchor = f'{chr}_anchor'
    output_read_f_name = os.path.join(individual_outputs_dir, f'{blast_input}_{anchor}_reads.fasta')
    print(f'Starting to make {output_read_f_name}...')

    df_new = df_best[df_best['anchor_name'] == anchor].sort_values(['read_length_past_anchor'], ascending=False)
    df_new.to_csv(os.path.join(individual_outputs_dir, f'{blast_input}_{anchor}_reads.tsv'), sep='\t', index=False)

    read_list_for_chr_end = zip(df_new['read_id'], df_new['repeat_type'])

    print(f'Grabbing reads in {anchor}...')
    split_reads_list = []

    for read_name, repeat_type in read_list_for_chr_end:
        sequence = fasta_index.fetch(read_name).strip('\n')
        row_index = df_best.loc[df_best['read_id'] == read_name].index[0]
        qstart_index = df_best.at[row_index, 'match_start_on_read'] - 1
        qend_index   = df_best.at[row_index, 'match_end_on_read'] - 1
        wanted_section = df_best.at[row_index, 'wanted_section_of_read']

        if wanted_section == 'before_match_start_on_read':
            chopped_seq = sequence[:qend_index]
        elif wanted_section == 'after_match_end_on_read':
            chopped_seq = sequence[qstart_index:]
        else:
            continue

        chopped_seq = chopped_seq.split('\n')[0]
        if len(chopped_seq) == 0:
            continue

        header = f'{read_name} {repeat_type} {chr}'
        split_reads_list.append(f'>{header}\n{chopped_seq}\n')

    with open(output_read_f_name, 'w') as f_out:
        f_out.writelines(split_reads_list)
