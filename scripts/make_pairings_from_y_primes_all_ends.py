# Make pairings for spacers

import pandas as pd
import os
import sys
from Bio import SeqIO


print("Starting make_pairings_from_y_primes.py")

input_base_name = sys.argv[1:]

#base_name = sys.argv[1]
anchor_set = 'test_anchors' #sys.argv[2]

for base_name in input_base_name:


    print(f'Opening {base_name}...')

    df_file = f'results/{base_name}_good_end_y_repeatmasker.tsv'
    input_chr_end_read_files_dir = f'results/{base_name}/chr_anchor_included_individual_files/'

    output_pairing_fasta_file_dir = f'results/{base_name}/paired_by_y_prime_reads/'

    os.system(f'mkdir -p {output_pairing_fasta_file_dir}')

    strain_id = f'{base_name.split("dorado_")[1]}'
    strain_id = f'{strain_id.split("_day")[0]}'

    df = pd.read_csv(df_file, sep='\t')


    chr_end_first_y_prime_id_dict = {}

    if '6991' == strain_id: # WildType
        y_prime_0_ends = ['chr1L', 'chr1R', 'chr2R', 'chr3L', 'chr3R', 'chr4L', 'chr6R',
                    'chr7L','chr9R', 'chr10R', 'chr11L', 'chr11R', 'chr13R', 'chr15L']
        chr_end_first_y_prime_id_dict = {
            'chr1L': None, 'chr1R': None,
            'chr2L': 'ID4', 'chr2R': None,
            'chr3L': None, 'chr3R': None,
            'chr4L': None, 'chr4R': 'ID2',
            'chr5L': 'ID6', 'chr5R': 'ID1',
            'chr6L': 'ID4', 'chr6R': None,
            'chr7L': None, 'chr7R': 'ID5',
            'chr8L': 'ID1', 'chr8R': 'ID1',
            'chr9L': 'ID6', 'chr9R': None,
            'chr10L': 'ID6', 'chr10R': None,
            'chr11L': None, 'chr11R': None,
            'chr12L': 'ID1', 'chr12R': 'ID1',
            'chr13L': 'ID1', 'chr13R': None,        
            'chr14L': 'ID5', 'chr14R': 'ID6',
            'chr15L': None, 'chr15R': 'ID2',
            'chr16L': 'ID5', 'chr16R': 'ID1'
            }
        
    elif '7172' == strain_id: # mph1
        y_prime_0_ends = ['chr1L', 'chr1R', 'chr2R', 'chr3L', 'chr3R', 'chr4L', 'chr6R',
                    'chr7L','chr9R', 'chr10R', 'chr11L', 'chr11R', 'chr13R', 'chr15L']
        chr_end_first_y_prime_id_dict = {
            'chr1L': None, 'chr1R': None,
            'chr2L': 'ID4', 'chr2R': None,
            'chr3L': None, 'chr3R': None,
            'chr4L': None, 'chr4R': 'ID2',
            'chr5L': 'ID6', 'chr5R': 'ID1',
            'chr6L': 'ID4', 'chr6R': None,
            'chr7L': None, 'chr7R': 'ID5',
            'chr8L': 'ID1', 'chr8R': 'ID1',
            'chr9L': 'ID6', 'chr9R': None,
            'chr10L': 'ID6', 'chr10R': None,
            'chr11L': None, 'chr11R': None,
            'chr12L': 'ID1', 'chr12R': 'ID1',
            'chr13L': 'ID1', 'chr13R': None,        
            'chr14L': 'ID5', 'chr14R': 'ID6',
            'chr15L': None, 'chr15R': 'ID2',
            'chr16L': 'ID5', 'chr16R': 'ID1'
            }    
        
    elif '7302' == strain_id: # mph1
        y_prime_0_ends = ['chr1L', 'chr1R', 'chr2R', 'chr3L', 'chr3R', 'chr4L', 'chr6R',
                    'chr7L','chr9R', 'chr10R', 'chr11L', 'chr11R', 'chr13R', 'chr15L']
        chr_end_first_y_prime_id_dict = {
            'chr1L': None, 'chr1R': None,
            'chr2L': 'ID4', 'chr2R': None,
            'chr3L': None, 'chr3R': None,
            'chr4L': None, 'chr4R': 'ID2',
            'chr5L': 'ID6', 'chr5R': 'ID1',
            'chr6L': 'ID4', 'chr6R': None,
            'chr7L': None, 'chr7R': 'ID5',
            'chr8L': 'ID1', 'chr8R': 'ID1',
            'chr9L': 'ID6', 'chr9R': None,
            'chr10L': 'ID6', 'chr10R': None,
            'chr11L': None, 'chr11R': None,
            'chr12L': 'ID1', 'chr12R': 'ID1',
            'chr13L': 'ID7', 'chr13R': None,        
            'chr14L': 'ID5', 'chr14R': 'ID6',
            'chr15L': None, 'chr15R': 'ID2',
            'chr16L': 'ID5', 'chr16R': 'ID1'
            }
        
    else:
        raise ValueError(f'Unknown strain_id: {strain_id}')

    def calculate_y_prime_change(df_read):
        #print(df_read)
        chr_end = df_read.at[0, 'chr_end']
        ref_first_y_prime_id = chr_end_first_y_prime_id_dict[chr_end]
        
        telomere_side = df_read.at[0, 'telomere_side']
        # Determine what the read's first Y' ID is
        if telomere_side == 'beginning':
            # The first Y' will have the largest start position on the read (it will be the last entry by default ordering)
            read_first_y_prime_id_and_color = df_read.at[(len(df_read)-1), 'y_prime_id_and_color']
            #print(read_first_y_prime_id_and_color)
            read_first_y_prime_id = read_first_y_prime_id_and_color.split('_')[0]
        else: # telomere_side == 'end':
            # The first Y' will have the smallest start position on the read (it will be the first entry by default ordering)
            read_first_y_prime_id_and_color = df_read.at[0, 'y_prime_id_and_color']
            #print(read_first_y_prime_id_and_color)            
            read_first_y_prime_id = read_first_y_prime_id_and_color.split('_')[0]
        
        if read_first_y_prime_id == ref_first_y_prime_id:
            y_prime_change = None
        else:
            y_prime_change = f'{chr_end}_and_{read_first_y_prime_id}'
        
        return y_prime_change


    df['y_prime_id_and_color'] = df['y_prime_group'].apply(lambda x: x.split('/')[2])

    all_reads = df['read_id'].unique()

    y_prime_change_dict = {}
    for read_id in all_reads:
        df_read = df[df['read_id'] == read_id].reset_index(drop=True)
        y_prime_change = calculate_y_prime_change(df_read)
        y_prime_change_dict[read_id] = y_prime_change

    df['y_prime_change'] = df['read_id'].apply(lambda x: y_prime_change_dict[x])

    df = df.dropna(subset=['y_prime_change'])

    y_prime_change_pairings = df['y_prime_change'].unique()
    for pairing in y_prime_change_pairings:
        df_pair = df[df['y_prime_change'] == pairing]
        pairing_read_list = df_pair['read_id'].unique()

        chr_end = pairing.split('_')[0]
        input_chr_end_read_file = f'{input_chr_end_read_files_dir}/{base_name}_blasted_{anchor_set}_{chr_end}_anchor_reads.fasta'
        output_pairing_fasta_file = f'{output_pairing_fasta_file_dir}/{pairing}.fasta'

        print(f'Creating file: {output_pairing_fasta_file}')

        with open(input_chr_end_read_file, 'r') as input_file, open(output_pairing_fasta_file, 'w') as output_file:
            # Parse the input file using SeqIO
            for record in SeqIO.parse(input_file, 'fasta'):
                header = record.description
                sequence = str(record.seq)

                chr_header = header.split(" ")[0]
                # Check if the header matches any string in the list
                if chr_header in pairing_read_list:
                    # Write the record to the output file
                    SeqIO.write(record, output_file, 'fasta')



