import pandas as pd
import pysam
import sys
import os

# --- MODIFIED: Capture explicit paths from Snakemake ---
# Order: python script.py {input.tsv} {input.fasta} {output.dir} {wildcard_base_name} {wildcard_anchor_set}
input_f_name_best = sys.argv[1]    # results/{BASE}/top_matches_...tsv
all_reads_f_name = sys.argv[2]     # path to the .fasta file
individual_outputs_dir = sys.argv[3] # results/{BASE}/chr_anchor_included_individual_files/
base_name = sys.argv[4]            # Need this for the internal file naming
anchor_set = sys.argv[5]           # Need this for the internal file naming

# This prefix is used for naming individual files inside the folder
file_prefix = f"{base_name}_blasted_{anchor_set}"

print(f"Starting split_and_label_all_reads_include_anchor.py")
print(f"Reading from: {input_f_name_best}")

# --- REMOVED: os.system(f'mkdir -p {individual_outputs_dir}') ---
# Snakemake handles directory creation automatically.

chr_list = ['chr10R', 'chr11R', 'chr12R', 'chr13R', 'chr14R', 'chr15R', 'chr16R', 'chr17R', 'chr1R', 
            'chr2R', 'chr3R', 'chr4R', 'chr5R', 'chr6R', 'chr7R', 'chr8R', 'chr9R', 'chr10L', 
            'chr11L', 'chr12L', 'chr13L', 'chr14L', 'chr15L', 'chr16L', 'chr17L', 'chr1L', 
            'chr2L', 'chr3L', 'chr4L', 'chr5L', 'chr6L', 'chr7L', 'chr8L', 'chr9L']

# Load TSV file
df_best = pd.read_csv(input_f_name_best, sep='\t')
fasta_index = pysam.FastaFile(all_reads_f_name)

for chr_name in chr_list:
    anchor = f'{chr_name}_anchor'
    
    # MODIFIED: Use the directory passed from Snakemake
    output_read_f_name = os.path.join(individual_outputs_dir, f'{file_prefix}_{anchor}_reads.fasta')
    output_tsv_f_name = os.path.join(individual_outputs_dir, f'{file_prefix}_{anchor}_reads.tsv')
    
    print(f'Starting to make {output_read_f_name}...')

    df_new = df_best[df_best['anchor_name'] == anchor]
    df_new = df_new.sort_values(['read_length_past_anchor'], ascending=False)
    
    # Save the individual TSV
    df_new.to_csv(output_tsv_f_name, sep="\t", index=False)

    read_list_for_chr_end = zip(df_new['read_id'], df_new['repeat_type'])

    split_reads_list = []
    for read_name, repeat_type in read_list_for_chr_end:
        sequence = fasta_index.fetch(read_name)
        sequence = sequence.strip('\n')
        
        row_index = df_best.loc[df_best['read_id'] == read_name].index[0]
        qstart_index = df_best.at[row_index,'match_start_on_read'] - 1
        qend_index = df_best.at[row_index,'match_end_on_read'] - 1
        wanted_section_of_read = df_best.at[row_index,'wanted_section_of_read']
        
        if wanted_section_of_read == 'before_match_start_on_read':
            chopped_seq = sequence[:qend_index]
        elif wanted_section_of_read == 'after_match_end_on_read':
            chopped_seq = sequence[qstart_index:]
        else:
            continue
            
        chopped_seq = chopped_seq.split('\n')[0]
        if len(chopped_seq) == 0:
            continue
        
        header = f'{read_name} {repeat_type} {chr_name}'
        split_reads_list.append(f'>{header}\n{chopped_seq}\n')

    with open(output_read_f_name, "w") as f_out:
        f_out.writelines(split_reads_list)

print("Done with split_and_label_all_reads_include_anchor.py")