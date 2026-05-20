import pandas as pd
import os
import sys
import re
import pysam

# Usage: fine_telomere_trimming.py <base_name> <anchor_set> <output_dir>
base_name  = sys.argv[1]
anchor_set = sys.argv[2]
output_dir = sys.argv[3]

print("Starting fine_telomere_trimming.py")

blast_dir    = os.path.join(output_dir, 'blast')
porechop_dir = os.path.join(output_dir, 'porechop')

input_reads_f_name      = os.path.join(porechop_dir, f'{base_name}_trimmed.fasta')
tag_and_adapter_info    = os.path.join(output_dir, f'{base_name}_adapter_trimming_check.tsv')
best_anchor_telo_f_name = os.path.join(blast_dir, f'top_matches_{base_name}_blasted_{anchor_set}.tsv')
output_tsv_f_name       = os.path.join(output_dir, f'{base_name}_post_telo_trimming.tsv')

repeat_trim_dir = os.path.join(output_dir, 'repeat_trim_files')
os.makedirs(repeat_trim_dir, exist_ok=True)

output_both_no_telo_repeat_read_f_name = os.path.join(repeat_trim_dir, f'{base_name}_final_no_telo_repeat_reads_best.fasta')
output_both_read_f_name                = os.path.join(repeat_trim_dir, f'{base_name}_final_telomere_repeat_reads_best.fasta')
output_both_trim_f_name                = os.path.join(repeat_trim_dir, f'{base_name}_final_telomere_repeat_outside_trim_reads_best.tsv')


def trim_for_telomere_repeat(read_sequence, repeat_type, minimum_repeat=15,
                              maximum_distance_for_seed_start=60,
                              small_telomere_gap_distance=3, large_telomere_gap_distance=20,
                              long_repeat_extension_threshold=20):

    read_sequence_length = len(read_sequence)

    if repeat_type == 'AC':
        telomere_repeat = r"([C]{1,3}A)+"
    else:
        telomere_repeat = r"([G]{1,3}T)+"
        read_sequence = read_sequence[::-1]

    telomere_start_checking      = True
    first_telomere_span          = True
    distance_from_start_of_read  = 0
    working_telomere_sequence    = read_sequence

    while telomere_start_checking:
        telomere_search = re.search(telomere_repeat, working_telomere_sequence)

        current_telomere_span        = telomere_search.span()
        current_telomere_span_length = current_telomere_span[1] - current_telomere_span[0]
        running_telomere_span        = [current_telomere_span[0] + distance_from_start_of_read,
                                        current_telomere_span[1] + distance_from_start_of_read]
        running_telomere_span_length = running_telomere_span[1] - running_telomere_span[0]

        if first_telomere_span:
            if running_telomere_span[0] > maximum_distance_for_seed_start:
                return None
            elif running_telomere_span_length >= minimum_repeat:
                start_of_telomere          = running_telomere_span[0]
                working_end_of_telomere    = running_telomere_span[1]
                first_telomere_span        = False
                distance_from_start_of_read = running_telomere_span[1]
                working_telomere_sequence  = read_sequence[distance_from_start_of_read:]
            else:
                distance_from_start_of_read = running_telomere_span[1]
                working_telomere_sequence  = read_sequence[distance_from_start_of_read:]
                continue
        else:
            gap = running_telomere_span[0] - working_end_of_telomere
            if gap <= small_telomere_gap_distance:
                working_end_of_telomere    = running_telomere_span[1]
                distance_from_start_of_read = running_telomere_span[1]
                working_telomere_sequence  = read_sequence[distance_from_start_of_read:]
                continue
            else:
                if current_telomere_span_length < long_repeat_extension_threshold:
                    end_of_telomere = working_end_of_telomere
                    break
                elif current_telomere_span_length >= long_repeat_extension_threshold and gap < large_telomere_gap_distance:
                    working_end_of_telomere    = running_telomere_span[1]
                    distance_from_start_of_read = running_telomere_span[1]
                    working_telomere_sequence  = read_sequence[distance_from_start_of_read:]
                    continue
                else:
                    end_of_telomere = working_end_of_telomere
                    break

    final_telomere_distance_to_telo_repeat = start_of_telomere
    final_telomere_repeat_length           = end_of_telomere - start_of_telomere

    if repeat_type == 'AC':
        final_telomere_repeat_sequence = read_sequence[start_of_telomere:end_of_telomere]
        final_telomere_outside_trim    = [read_sequence[:start_of_telomere],
                                          read_sequence[end_of_telomere:(end_of_telomere + 10)]]
    else:
        read_sequence  = read_sequence[::-1]
        start_of_telo  = read_sequence_length - end_of_telomere
        end_of_telo    = read_sequence_length - start_of_telomere
        final_telomere_repeat_sequence = read_sequence[start_of_telo:end_of_telo]
        final_telomere_outside_trim    = [read_sequence[(start_of_telo - 10):start_of_telo],
                                          read_sequence[end_of_telo:]]

    return (final_telomere_repeat_sequence, final_telomere_repeat_length,
            final_telomere_outside_trim, final_telomere_distance_to_telo_repeat)


# Index the trimmed fasta
print(f'Making file and Indexing Reads...')
os.system(f'samtools faidx {input_reads_f_name}')
fasta_index = pysam.FastaFile(input_reads_f_name)

df_info = pd.read_csv(tag_and_adapter_info, sep='\t', index_col=[0])
columns_to_include = ['read_id', 'Repeat_Type', 'Adapter_After_Telomere', 'Tag_After_Telomere',
                      'both_adapter_and_tag', 'either_adapter_or_tag', 'anchor_name']
df_info_filtered = df_info[columns_to_include].copy()
print(df_info_filtered)

chr_list = (
    [f'chr{n}R' for n in range(1, 17)] +
    [f'chr{n}L' for n in range(1, 17)]
)

running_both_no_telo_repeat_reads_list = []
running_both_to_telomere_reads_list    = []
running_outside_trimmed_reads_list     = []
add_to_tsv_list                        = []

for chr in chr_list:
    output_no_telo_AC = os.path.join(repeat_trim_dir, f'{chr}_telomere_no_telo_repeat_AC_reads_best.fasta')
    output_no_telo_TG = os.path.join(repeat_trim_dir, f'{chr}_telomere_no_telo_repeat_TG_reads_best.fasta')
    output_AC         = os.path.join(repeat_trim_dir, f'{chr}_telomere_trimmed_AC_reads_best.fasta')
    output_TG         = os.path.join(repeat_trim_dir, f'{chr}_telomere_trimmed_TG_reads_best.fasta')

    df_chr     = df_info_filtered[df_info_filtered['anchor_name'] == f'{chr}_anchor']
    read_list  = df_chr['read_id'].tolist()

    print(f'Starting to make {chr}...')
    print(f'Grabbing reads in {chr}...')

    no_telo_repeat_reads_AC = []
    no_telo_repeat_reads_TG = []
    trimmed_reads_AC        = []
    trimmed_reads_TG        = []
    outside_trimmed_reads   = []

    for read_name in read_list:
        sequence    = fasta_index.fetch(read_name).strip('\n')
        row_index   = df_info_filtered.loc[df_info_filtered['read_id'] == read_name].index[0]
        repeat_type = df_info_filtered.at[row_index, 'Repeat_Type']

        if repeat_type == 'AC':
            try:
                trimming_results = trim_for_telomere_repeat(sequence, 'AC')
            except AttributeError:
                continue
            if trimming_results is None:
                header = read_name + ' AC ' + chr + '\n'
                no_telo_repeat_reads_AC += [header, sequence + '\n']
                repeat_length = length_to_trim = None
            else:
                header = read_name + ' AC ' + chr + '\n'
                trimmed_read, repeat_length, outside_trim, length_to_trim = trimming_results
                trimmed_reads_AC      += [header, trimmed_read + '\n']
                outside_trimmed_reads += [header, str(outside_trim) + '\n', trimmed_read + '\n']

        elif repeat_type == 'TG':
            try:
                trimming_results = trim_for_telomere_repeat(sequence, 'TG')
            except AttributeError:
                continue
            if trimming_results is None:
                header = read_name + ' TG ' + chr + '\n'
                no_telo_repeat_reads_TG += [header, sequence + '\n']
                repeat_length = length_to_trim = None
            else:
                header = read_name + ' TG ' + chr + '\n'
                trimmed_read, repeat_length, outside_trim, length_to_trim = trimming_results
                trimmed_reads_TG      += [header, trimmed_read + '\n']
                outside_trimmed_reads += [header, str(outside_trim) + '\n', trimmed_read + '\n']
        else:
            continue

        add_to_tsv_list.append([read_name, repeat_length, length_to_trim])

    with open(output_AC, 'w')         as f: f.writelines(trimmed_reads_AC)
    with open(output_TG, 'w')         as f: f.writelines(trimmed_reads_TG)
    with open(output_no_telo_AC, 'w') as f: f.writelines(no_telo_repeat_reads_AC)
    with open(output_no_telo_TG, 'w') as f: f.writelines(no_telo_repeat_reads_TG)

    running_both_no_telo_repeat_reads_list = (no_telo_repeat_reads_AC + no_telo_repeat_reads_TG
                                               + running_both_no_telo_repeat_reads_list)
    running_both_to_telomere_reads_list    = (trimmed_reads_AC + trimmed_reads_TG
                                               + running_both_to_telomere_reads_list)
    running_outside_trimmed_reads_list     = outside_trimmed_reads + running_outside_trimmed_reads_list

with open(output_both_no_telo_repeat_read_f_name, 'w') as f: f.writelines(running_both_no_telo_repeat_reads_list)
with open(output_both_read_f_name, 'w')               as f: f.writelines(running_both_to_telomere_reads_list)
with open(output_both_trim_f_name, 'w')               as f: f.writelines(running_outside_trimmed_reads_list)

df_trimmed_repeat = pd.DataFrame(add_to_tsv_list, columns=['read_id', 'repeat_length', 'length_to_trim_end_of_read'])
df_trimmed_info   = pd.merge(df_trimmed_repeat, df_info_filtered, on='read_id', how='inner')

df_best_anchor_overhang          = pd.read_csv(best_anchor_telo_f_name, sep='\t')
df_best_anchor_overhang_selected = df_best_anchor_overhang[['read_id', 'read_length_past_anchor']].copy()

df_trimmed_and_overhang = pd.merge(df_trimmed_info, df_best_anchor_overhang_selected, on='read_id', how='inner')

def trim_read_lengths_past_anchor(read_delta, trim_length):
    return trim_length if trim_length is None else read_delta - trim_length

df_trimmed_and_overhang['trimmed_read_length_past_anchor'] = df_trimmed_and_overhang.apply(
    lambda r: trim_read_lengths_past_anchor(r['read_length_past_anchor'], r['length_to_trim_end_of_read']), axis=1)

columns_in_order = ['read_id', 'anchor_name', 'Repeat_Type', 'repeat_length',
                    'Adapter_After_Telomere', 'Tag_After_Telomere',
                    'both_adapter_and_tag', 'either_adapter_or_tag',
                    'trimmed_read_length_past_anchor', 'length_to_trim_end_of_read']

df_final = df_trimmed_and_overhang[columns_in_order].copy()
print(df_final)
df_final.to_csv(output_tsv_f_name, sep='\t')