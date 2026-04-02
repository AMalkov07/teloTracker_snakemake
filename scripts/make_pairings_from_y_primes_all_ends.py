# Make pairings for spacers

import pandas as pd
import os
import re
import sys
from Bio import SeqIO


def build_chr_end_first_y_prime_id_dict(y_prime_lib_fasta):
    """Build mapping of chr_end -> first Y prime ID from library FASTA.

    Parses headers like:
        >Y_Prime_chr4R1,2,3,6,7;chr12R6,7#Long/Tandem/ID2_Red
    For each chr_end, finds the Y prime entry that contains position 1.
    Chr_ends with no Y primes map to None.
    """
    pos1_id = {}  # chr_end -> ID for position 1
    all_chr_ends_with_yprimes = set()

    with open(y_prime_lib_fasta) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            header = line.strip().lstrip('>')
            if '#' not in header:
                continue
            prefix, metadata = header.split('#', 1)
            parts = metadata.split('/')
            if len(parts) < 3:
                continue
            id_color = parts[2]
            y_id = id_color.split('_', 1)[0]

            origin_body = prefix.replace('Y_Prime_', '')
            for group in origin_body.split(';'):
                gm = re.match(r'(chr\d+[LR])([\d,]+)', group)
                if gm:
                    chr_end = gm.group(1)
                    positions = gm.group(2).split(',')
                    all_chr_ends_with_yprimes.add(chr_end)
                    if '1' in positions:
                        pos1_id[chr_end] = y_id

    # Build full dict for all 32 chr_ends
    result = {}
    for n in range(1, 17):
        for s in ('L', 'R'):
            ce = f'chr{n}{s}'
            result[ce] = pos1_id.get(ce, None)

    # Also derive y_prime_0_ends
    all_possible = set(f'chr{n}{s}' for n in range(1, 17) for s in ('L', 'R'))
    y_prime_0_ends = sorted(all_possible - all_chr_ends_with_yprimes)

    return result, y_prime_0_ends


print("Starting make_pairings_from_y_primes.py")

input_base_name = sys.argv[1:]

#base_name = sys.argv[1]
anchor_set = 'test_anchors' #sys.argv[2]

for base_name in input_base_name:


    print(f'Opening {base_name}...')

    df_file = f'results/{base_name}/{base_name}_good_end_y_repeatmasker.tsv'
    input_chr_end_read_files_dir = f'results/{base_name}/chr_anchor_included_individual_files/'

    output_pairing_fasta_file_dir = f'results/{base_name}/paired_by_y_prime_reads/'

    os.system(f'mkdir -p {output_pairing_fasta_file_dir}')

    strain_id = f'{base_name.split("dorado_")[1]}'
    strain_id = f'{strain_id.split("_day")[0]}'

    df = pd.read_csv(df_file, sep='\t')

    # Build chr_end -> first Y prime ID mapping dynamically from library FASTA
    y_prime_lib_fasta = f'references/extracted_yprimes_{strain_id}.fasta'
    if not os.path.exists(y_prime_lib_fasta):
        raise FileNotFoundError(f'Y prime library not found: {y_prime_lib_fasta}')

    chr_end_first_y_prime_id_dict, y_prime_0_ends = build_chr_end_first_y_prime_id_dict(y_prime_lib_fasta)
    print(f'  Built first Y prime ID dict for {sum(1 for v in chr_end_first_y_prime_id_dict.values() if v is not None)} chr_ends')

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



