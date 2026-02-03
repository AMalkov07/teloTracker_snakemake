# Make pairings for spacers

import pandas as pd
import os
import sys
from Bio import SeqIO


print("Starting make_pairings_from_y_primes.py")

base_name = sys.argv[1]

anchor_set = sys.argv[2]

strain_id = sys.argv[3]


print(f'Opening {base_name}...')

df_file = f'results/{base_name}_gained_y_repeatmasker.tsv'
input_chr_end_read_files_dir = f'results/{base_name}/chr_anchor_included_individual_files/'

output_pairing_fasta_file_dir = f'results/{base_name}/paired_by_y_prime_reads/'

os.system(f'mkdir -p {output_pairing_fasta_file_dir}')



df = pd.read_csv(df_file, sep='\t')

df_0_to_1 = df[(df['reference_y_primes'] == 0) & (df['repeatmasker_y_prime_count'] == 1)]


if '6991' == strain_id: # WildType
    y_prime_0_ends = ['chr1L', 'chr1R', 'chr2R', 'chr3L', 'chr3R', 'chr4L', 'chr6R',
                'chr7L','chr9R', 'chr10R', 'chr11L', 'chr11R', 'chr13R', 'chr15L']

elif '7172' == strain_id: # mph1
    y_prime_0_ends = ['chr1L', 'chr1R', 'chr2R', 'chr3L', 'chr3R', 'chr4L', 'chr6R',
                'chr7L','chr9R', 'chr10R', 'chr11L', 'chr11R', 'chr13R', 'chr15L']

elif '7302' == strain_id: # mph1
    y_prime_0_ends = ['chr1L', 'chr1R', 'chr2R', 'chr3L', 'chr3R', 'chr4L', 'chr6R',
                'chr7L','chr9R', 'chr10R', 'chr11L', 'chr11R', 'chr13R', 'chr15L']
    
else:
    raise ValueError(f'Unknown strain_id: {strain_id}')


for y_prime_0_end in y_prime_0_ends:
    df_y_prime_0_end = df_0_to_1[df_0_to_1['chr_end'] == y_prime_0_end]
    y_prime_0_end_and_color_list = df_y_prime_0_end['y_prime_color'].unique()

    for color_tag in y_prime_0_end_and_color_list:
        df_y_prime_0_end_and_color = df_y_prime_0_end[df_y_prime_0_end['y_prime_color'] == color_tag]
        y_prime_0_end_and_color_read_list = df_y_prime_0_end_and_color['read_id'].unique()

        y_prime_id = color_tag.split('_')[0]
        input_chr_end_read_file = f'{input_chr_end_read_files_dir}/{base_name}_blasted_{anchor_set}_{y_prime_0_end}_anchor_reads.fasta'
        output_pairing_fasta_file = f'{output_pairing_fasta_file_dir}/{y_prime_0_end}_and_{y_prime_id}.fasta'

        print(f'Creating file: {output_pairing_fasta_file}')

        with open(input_chr_end_read_file, 'r') as input_file, open(output_pairing_fasta_file, 'w') as output_file:
            # Parse the input file using SeqIO
            for record in SeqIO.parse(input_file, 'fasta'):
                header = record.description
                sequence = str(record.seq)

                chr_header = header.split(" ")[0]
                #print(chr_header)
                #print(y_prime_0_end_and_color_read_list)
                # Check if the header matches any string in the list
                if chr_header in y_prime_0_end_and_color_read_list:
                    # Write the record to the output file
                    SeqIO.write(record, output_file, 'fasta')


