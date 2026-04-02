import pandas as pd
import re
import sys


def build_y_prime_order_dict(y_prime_lib_fasta, labeled_bed_path):
    """Build y_prime_order_dict dynamically from Y prime library FASTA and BED file.

    For each chr_end, determines the ordered tuple of Y prime group IDs
    (anchor-proximal to telomere-proximal) by:
    1. Parsing library FASTA headers to map (chr_end, position) -> ID
    2. Parsing the BED file for Y_Prime features at each chr_end
    3. Looking up each Y prime's group ID

    Returns dict: chr_end -> tuple of ID strings, or None if no Y primes.
    """
    # Step 1: Parse library FASTA headers to build (chr_end, position) -> ID mapping
    pos_to_id = {}  # (chr_end, position_str) -> ID
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
            y_id = id_color.split('_', 1)[0]  # e.g., 'ID2'

            # Parse locations from prefix: Y_Prime_chr4R1,2,3,6,7;chr12R6,7
            origin_body = prefix.replace('Y_Prime_', '')
            for group in origin_body.split(';'):
                gm = re.match(r'(chr\d+[LR])([\d,]+)', group)
                if gm:
                    chr_end = gm.group(1)
                    positions = gm.group(2).split(',')
                    for pos in positions:
                        pos_to_id[(chr_end, pos)] = y_id

    # Step 2: Parse BED file for Y_Prime features
    chr_end_yprimes = {}  # chr_end -> [(position_num, ID)]
    all_chr_ends = set()
    with open(labeled_bed_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            feat_name = fields[3]
            m = re.match(r'(chr\d+[LR])_', feat_name)
            if m:
                all_chr_ends.add(m.group(1))
            # Match Y_Prime features: e.g., chr4R_Y_Prime_3
            m = re.match(r'(chr\d+[LR])_Y_Prime_(\d+)', feat_name)
            if not m:
                continue
            chr_end = m.group(1)
            pos_num = int(m.group(2))
            pos_str = str(pos_num)

            y_id = pos_to_id.get((chr_end, pos_str), feat_name)
            if chr_end not in chr_end_yprimes:
                chr_end_yprimes[chr_end] = []
            chr_end_yprimes[chr_end].append((pos_num, y_id))

    # Step 3: Build the order dict
    all_possible = [f'chr{n}{s}' for n in range(1, 17) for s in ('L', 'R')]
    all_chr_ends.update(all_possible)

    y_prime_order_dict = {}
    for ce in sorted(all_chr_ends):
        if ce not in chr_end_yprimes or not chr_end_yprimes[ce]:
            y_prime_order_dict[ce] = None
        else:
            sorted_yps = sorted(chr_end_yprimes[ce], key=lambda x: x[0])
            y_prime_order_dict[ce] = tuple(y_id for _, y_id in sorted_yps)

    return y_prime_order_dict


def calculate_y_prime_change(read_id, chr_end, df_y_repeatmasker, y_prime_order_dict):
    y_prime_ref_order = y_prime_order_dict[chr_end]
    recombination_events = {}
    
    if read_id not in df_y_repeatmasker['read_id'].values:
        y_prime_match_order = []
        better_repeatmasker_y_prime_count = 0    
    else:
        df_read = df_y_repeatmasker[df_y_repeatmasker['read_id'] == read_id].reset_index(drop=True)
        telomere_side = df_read.at[0, 'telomere_side']

        y_prime_match_order = df_read['y_prime_id_and_match'].unique()
        better_repeatmasker_y_prime_count = len(y_prime_match_order)

        # Determine what the read's first Y' ID is
        if telomere_side == 'beginning':
            # The first Y' will have the largest start position on the read
            y_prime_match_order = y_prime_match_order[::-1]
        # else: telomere_side == 'end', order stays as is

    for i, y_prime_and_match_id in enumerate(y_prime_match_order):
        current_y_in_read = y_prime_and_match_id.split("-")[0]

        # Catch exceptions
        if y_prime_ref_order == None or i >= len(y_prime_ref_order):
            expected_y_in_ref = None
        else:
            expected_y_in_ref = y_prime_ref_order[i]

        if current_y_in_read == expected_y_in_ref:
            pass
        else:
            recombination_events[f'read_y_num_{i+1}'] = f'from {expected_y_in_ref} to {current_y_in_read}'

    # Include losses of Y's
    if y_prime_ref_order != None:
        if len(y_prime_match_order) < len(y_prime_ref_order):
            for i in range(len(y_prime_match_order), len(y_prime_ref_order)):
                expected_y_in_ref = y_prime_ref_order[i]
                recombination_events[f'ref_y_num_{i+1}'] = f'from {expected_y_in_ref} to {None}'

    # Check if it is the day0 reference order (no recombination)
    if len(recombination_events) == 0:
        recombination_status = "No Change"
    # Check the earliest place of recombination (1st Y' or elsewhere)
    else:
        if 'read_y_num_1' in recombination_events.keys() or 'ref_y_num_1' in recombination_events.keys():
            recombination_status = "1st Y' Change"
        else:
            recombination_status = "Y' Recombination"

    return recombination_events, recombination_status, better_repeatmasker_y_prime_count

print("Starting get_stats_of_recombination.py")

# Parse command line arguments
good_end_y_file = sys.argv[1]
y_prime_probe_file = sys.argv[2]
strain_id = sys.argv[3]
output_file = sys.argv[4]
y_prime_lib_fasta = sys.argv[5]
labeled_bed_path = sys.argv[6]

print(f'Opening {good_end_y_file}...')

# Build Y prime order dict dynamically from library FASTA + BED
y_prime_order_dict = build_y_prime_order_dict(y_prime_lib_fasta, labeled_bed_path)
print(f"Built Y prime order dict for {sum(1 for v in y_prime_order_dict.values() if v is not None)} chr_ends with Y primes")

# Read input files
df_y_repeatmasker = pd.read_csv(good_end_y_file, sep='\t')
df_input = pd.read_csv(y_prime_probe_file, sep='\t')

# Filter input data
df_input = df_input.dropna(subset=["repeat_length"])
df_input = df_input[df_input['Adapter_After_Telomere'] == True]

# Process Y prime data
df_y_repeatmasker['y_prime_id_and_color'] = df_y_repeatmasker['y_prime_group'].apply(lambda x: x.split('/')[2])
df_y_repeatmasker['y_prime_id'] = df_y_repeatmasker['y_prime_id_and_color'].apply(lambda x: x.split('_')[0])
df_y_repeatmasker['y_prime_id_and_match'] = (df_y_repeatmasker['y_prime_id'].astype(str) + '-' + 
                                               df_y_repeatmasker['match_id'].astype(str))

all_reads = df_input['read_id'].unique()

y_prime_recomb_events_dict = {}
y_prime_recomb_status_dict = {}
better_repeatmasker_y_prime_count_dict = {}

for read_id in all_reads:
    chr_end = f"chr{df_input[df_input['read_id'] == read_id].iloc[0]['chr_end']}"
    y_prime_recombinations, recombination_status, better_repeatmasker_y_prime_count = calculate_y_prime_change(
        read_id, chr_end, df_y_repeatmasker, y_prime_order_dict)
    y_prime_recomb_events_dict[read_id] = y_prime_recombinations
    y_prime_recomb_status_dict[read_id] = recombination_status
    better_repeatmasker_y_prime_count_dict[read_id] = better_repeatmasker_y_prime_count

df_input['y_prime_recombination_events'] = df_input['read_id'].apply(lambda x: y_prime_recomb_events_dict[x])
df_input['y_prime_recombination_status'] = df_input['read_id'].apply(lambda x: y_prime_recomb_status_dict[x])

#######################
df_input['better_repeatmasker_y_prime_count'] = df_input['read_id'].apply(lambda x: better_repeatmasker_y_prime_count_dict[x])
df_input['better_repeatmasker_y_primes_relative_to_ref'] = (df_input['better_repeatmasker_y_prime_count'] - 
                                                              df_input['reference_y_primes'])
#######################

def calc_delta_group(better_repeatmasker_y_primes_relative_to_ref):
    if better_repeatmasker_y_primes_relative_to_ref == 0:
        delta_group = 'Same Number'
    elif better_repeatmasker_y_primes_relative_to_ref < 0:
        delta_group = 'Loss'
    elif better_repeatmasker_y_primes_relative_to_ref == 1:
        delta_group = 'Gain 1'
    else:
        delta_group = 'Gain Multiple'
    return delta_group

df_input['y_prime_delta_group'] = df_input['better_repeatmasker_y_primes_relative_to_ref'].apply(lambda x: calc_delta_group(x))

print(f'Total Reads: {len(df_input)}\n')

# Save output
df_input.to_csv(output_file, sep='\t', index=False)

# Print summary statistics
total_reads = len(df_input)
print(f'Total Reads: {total_reads}\n')

# 1. y_prime_recombination_status
counts = df_input['y_prime_recombination_status'].value_counts()
percentages = (counts / total_reads * 100).apply(lambda x: f"{x:.2f}%")
print("y_prime_recombination_status Counts:")
print(counts)
print("\nPercentages (of total reads):")
print(percentages)
print("\n")

# 2. reference_y_prime_end_status within each recombination status
grouped1 = df_input.groupby('y_prime_recombination_status')['reference_y_prime_end_status'].value_counts()
percentages = (grouped1 / total_reads * 100).apply(lambda x: f"{x:.2f}%")
print("reference_y_prime_end_status by y_prime_recombination_status:")
print(grouped1)
print("\nPercentages (of total reads):")
print(percentages)
print("\n")

# 3. y_prime_delta_group within each recombination + end status combo
grouped2 = df_input.groupby(['y_prime_recombination_status', 'reference_y_prime_end_status'])['y_prime_delta_group'].value_counts()
percentages = (grouped2 / total_reads * 100).apply(lambda x: f"{x:.2f}%")
print("y_prime_delta_group by recombination status + reference end status:")
print(grouped2)
print("\nPercentages (of total reads):")
print(percentages)
print("\n")

print(f"Finished. Output saved to: {output_file}")