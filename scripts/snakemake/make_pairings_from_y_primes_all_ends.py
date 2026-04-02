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
    pos1_id = {}
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
                    if '1' in positions:
                        pos1_id[chr_end] = y_id

    result = {}
    for n in range(1, 17):
        for s in ('L', 'R'):
            ce = f'chr{n}{s}'
            result[ce] = pos1_id.get(ce, None)
    return result


print("Starting make_pairings_from_y_primes_all_ends.py")

# Parse command line arguments
if len(sys.argv) < 7:
    print(f"Usage: {sys.argv[0]} <good_end_y_file> <input_chr_reads_dir> <output_pairings_dir> <strain_id> <anchor_set> <base_name> <y_prime_lib_fasta>")
    print(f"Got {len(sys.argv)} arguments: {sys.argv}")
    sys.exit(1)

good_end_y_file = sys.argv[1]
input_chr_end_read_files_dir = sys.argv[2]
output_pairing_fasta_file_dir = sys.argv[3]
strain_id = sys.argv[4]
anchor_set = sys.argv[5]
base_name = sys.argv[6]
y_prime_lib_fasta = sys.argv[7]

print(f'Opening {good_end_y_file}...')
print(f'Input directory: {input_chr_end_read_files_dir}')
print(f'Output directory: {output_pairing_fasta_file_dir}')
print(f'Strain: {strain_id}, Anchor set: {anchor_set}')
print(f'Base name: {base_name}')

# Read input data
df = pd.read_csv(good_end_y_file, sep='\t')

# Build chr_end -> first Y prime ID mapping dynamically from library FASTA
chr_end_first_y_prime_id_dict = build_chr_end_first_y_prime_id_dict(y_prime_lib_fasta)
print(f'Built first Y prime ID dict for {sum(1 for v in chr_end_first_y_prime_id_dict.values() if v is not None)} chr_ends')

def calculate_y_prime_change(df_read, chr_end_first_y_prime_id_dict):
    chr_end = df_read.at[0, 'chr_end']
    ref_first_y_prime_id = chr_end_first_y_prime_id_dict[chr_end]
    
    telomere_side = df_read.at[0, 'telomere_side']
    # Determine what the read's first Y' ID is
    if telomere_side == 'beginning':
        # The first Y' will have the largest start position on the read
        read_first_y_prime_id_and_color = df_read.at[(len(df_read)-1), 'y_prime_id_and_color']
        read_first_y_prime_id = read_first_y_prime_id_and_color.split('_')[0]
    else:  # telomere_side == 'end':
        # The first Y' will have the smallest start position on the read
        read_first_y_prime_id_and_color = df_read.at[0, 'y_prime_id_and_color']
        read_first_y_prime_id = read_first_y_prime_id_and_color.split('_')[0]
    
    if read_first_y_prime_id == ref_first_y_prime_id:
        y_prime_change = None
    else:
        y_prime_change = f'{chr_end}_and_{read_first_y_prime_id}'
    
    return y_prime_change

# Process Y prime data
df['y_prime_id_and_color'] = df['y_prime_group'].apply(lambda x: x.split('/')[2])

all_reads = df['read_id'].unique()

y_prime_change_dict = {}
for read_id in all_reads:
    df_read = df[df['read_id'] == read_id].reset_index(drop=True)
    y_prime_change = calculate_y_prime_change(df_read, chr_end_first_y_prime_id_dict)
    y_prime_change_dict[read_id] = y_prime_change

df['y_prime_change'] = df['read_id'].apply(lambda x: y_prime_change_dict[x])

df = df.dropna(subset=['y_prime_change'])

y_prime_change_pairings = df['y_prime_change'].unique()
print(f'Found {len(y_prime_change_pairings)} unique pairings')

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
            
            chr_header = header.split(" ")[0]
            # Check if the header matches any string in the list
            if chr_header in pairing_read_list:
                # Write the record to the output file
                SeqIO.write(record, output_file, 'fasta')

print(f"Finished. Created {len(y_prime_change_pairings)} pairing files in {output_pairing_fasta_file_dir}")