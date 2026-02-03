import pandas as pd
import os
import sys
import re
import pysam


base_name=f'{sys.argv[1]}'

anchor_set=f'{sys.argv[2]}'    # "new_all_best_anchors"  # f'{sys.argv[2]}' # Default = new_all_best_anchors

print("Starting fine_telomere_trimming.py")


def trim_for_telomere_repeat(read_sequence,repeat_type,minimum_repeat=15,maximum_distance_for_seed_start=60,
                                small_telomere_gap_distance=3, large_telomere_gap_distance=20, long_repeat_extension_threshold=20):

    read_sequence_length = len(read_sequence)

    # Assign Telomere Repeat to look for
    if repeat_type == 'AC':
        telomere_repeat = r"([C]{1,3}A)+"
    else: # repeat_type == 'TG'
        # Reverse of repeat as we will reverse the string when searching
        telomere_repeat = r"([G]{1,3}T)+"
        # Need to reverse the TG read_name sequence for the search
        read_sequence = read_sequence[::-1]
    #print(read_sequence[0:50])
    #print(telomere_repeat)

    telomere_start_checking = True
    first_telomere_span = True
    distance_from_start_of_read = 0
    working_telomere_sequence = read_sequence

    # Continually gets span legnths and uses start_of_telomere and end_of_telomere as markers for total length
    while telomere_start_checking == True:
        telomere_search = re.search(telomere_repeat,working_telomere_sequence)

        # The basepair positions on the truncated read
        # Only used in checking for extension
        current_telomere_span = telomere_search.span()
        current_telomere_span_length = current_telomere_span[1] - current_telomere_span[0]

        # The basepair positions on the orignal read
        running_telomere_span = [(current_telomere_span[0]+distance_from_start_of_read) , (current_telomere_span[1]+distance_from_start_of_read)]
        running_telomere_span_length = running_telomere_span[1] - running_telomere_span[0]

        #print(running_telomere_span_length)
        if first_telomere_span == True:
            if running_telomere_span[0] > maximum_distance_for_seed_start: # Default is 60 basepairs
                return None
            # Default is 15 basepairs
            elif running_telomere_span_length >= minimum_repeat:
                start_of_telomere = running_telomere_span[0]
                working_end_of_telomere = running_telomere_span[1]
                first_telomere_span = False
                distance_from_start_of_read = (running_telomere_span[1])
                working_telomere_sequence = read_sequence[distance_from_start_of_read:]
            else:
                distance_from_start_of_read = (running_telomere_span[1])
                working_telomere_sequence = read_sequence[distance_from_start_of_read:]
                continue # The sequence is considered the start of the telomere
        else:
            gap_between_telomere_repeats = running_telomere_span[0] - working_end_of_telomere
            #print(f'Gap = {gap_between_telomere_repeats}')
            #print(running_telomere_span)

            # Continue extension of telomere repeat if gap is smaller than small_telomere_gap_distance (3 basepairs)
            if gap_between_telomere_repeats <= small_telomere_gap_distance:
                working_end_of_telomere = running_telomere_span[1]
                distance_from_start_of_read = (running_telomere_span[1])
                working_telomere_sequence = read_sequence[distance_from_start_of_read:]
                continue
            # Break extension if gap is larger than small_telomere_gap_distance and:
            # 1. current_telomere_span_length < long_repeat_extension_threshold
            # or 2. gap_between_telomere_repeats > large_telomere_gap_distance
            # If both false (the elif case), then extend still
            else:
                if current_telomere_span_length < long_repeat_extension_threshold:
                    end_of_telomere = working_end_of_telomere
                    break
                elif current_telomere_span_length >= long_repeat_extension_threshold and gap_between_telomere_repeats < large_telomere_gap_distance:
                    working_end_of_telomere = running_telomere_span[1]
                    distance_from_start_of_read = (running_telomere_span[1])
                    working_telomere_sequence = read_sequence[distance_from_start_of_read:]
                    continue
                else:
                    end_of_telomere = working_end_of_telomere
                    break

    # Get Telomere Repeat Info
    final_telomere_distance_to_telo_repeat = start_of_telomere
    final_telomere_repeat_length = end_of_telomere - start_of_telomere

    # Assign Telomere Repeat Sequences
    if repeat_type == 'AC':
        final_telomere_repeat_sequence = read_sequence[start_of_telomere:end_of_telomere]
        final_telomere_outside_trim = [read_sequence[:start_of_telomere], read_sequence[end_of_telomere:(end_of_telomere+10)]]
    else: # repeat_type == 'TG'
        # Need to correct for reverse the TG read_name sequence for the search
        read_sequence = read_sequence[::-1]
        start_of_telo = (read_sequence_length - end_of_telomere) # Changes the index value for the reversed sequence and makes new variable
        end_of_telo = (read_sequence_length - start_of_telomere) # Changes the index value for the reversed sequence and makes new variable
        # Assign Telomere Repeat
        final_telomere_repeat_sequence = read_sequence[start_of_telo:end_of_telo]
        final_telomere_outside_trim = [read_sequence[(start_of_telo-10):start_of_telo], read_sequence[end_of_telo:]]


    #print(f'Final Span = {final_telomere_repeat_span}')
    #print(f'Final Length = {final_telomere_repeat_length}')

    return (final_telomere_repeat_sequence, final_telomere_repeat_length, final_telomere_outside_trim, final_telomere_distance_to_telo_repeat)

def check_for_adpt_or_tag(adapter_bool, tag_bool) -> bool:
    if adapter_bool == True or tag_bool == True:
        either_bool = True
    else:
        either_bool = False
    return either_bool


input_reads_f_name=f'results/{base_name}/{base_name}_all_chromosome_anchored_reads_pre_trim.fasta'

tag_and_adapter_info_f_name=f'results/{base_name}/{base_name}_adapter_trimming_check.tsv'
best_anchor_telo_f_name=f'results/{base_name}/top_matches_{base_name}_blasted_{anchor_set}.tsv'

output_tsv_f_name=f'results/{base_name}/{base_name}_post_telo_trimming.tsv'

repeat_trim_dir = f'results/{base_name}/repeat_trim_files/'
output_both_no_telo_repeat_read_f_name=f'{repeat_trim_dir}/{base_name}_final_no_telo_repeat_reads_best.fasta'
output_both_read_f_name=f'{repeat_trim_dir}/{base_name}_final_telomere_repeat_reads_best.fasta'
output_both_trim_f_name=f'{repeat_trim_dir}/{base_name}_final_telomere_repeat_outside_trim_reads_best.tsv'


os.system(f'mkdir -p {repeat_trim_dir}')

# Load the FAIDX file
print(f'Making file and Indexing Reads...')
if os.path.exists(f'{input_reads_f_name}.fai'):

    ####### Remove later, if you need don't think your index file is messed up ##########
    os.system(f'samtools faidx {input_reads_f_name}')

    pass
else:
    os.system(f'samtools faidx {input_reads_f_name}')

fasta_index = pysam.FastaFile(input_reads_f_name)

df_info = pd.read_csv(f'{tag_and_adapter_info_f_name}', sep ='\t', index_col=[0])
# Names from blast = ['qseqid', 'qlen', 'length', 'qstart', 'qend', 'sseqid', 'slen', 'sstart', 'send', 'pident', 'bitscore', 'evalue', 'wanted_section', 'l_end', 'overhang_length', 'correct_overhang'])
# ['read_id', 'total_read_length', 'read_bp_used_for_match', 'match_start_on_read', 'match_end_on_read', 'anchor_name', 'total_anchor_length', 'match_start_on_anchor', 'match_end_on_anchor', 'pident', 'bitscore', 'evalue', 'wanted_section_of_read', 'l_end_chr', 'read_length_past_anchor', 'expected_wt_overhang_distance'])

columns_to_include = ['read_id', 'Repeat_Type', 'Adapter_After_Telomere', 'Tag_After_Telomere', 'both_adapter_and_tag', 'either_adapter_or_tag', 'anchor_name']

df_info_filtered = df_info[columns_to_include].copy()
print(df_info_filtered)


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


running_both_no_telo_repeat_reads_list = []
running_both_to_telomere_reads_list = []
running_outside_trimmed_reads_list = []
add_to_tsv_list = []

for chr in chr_list:


    output_no_telo_repeat_AC_read_f_name=f'{repeat_trim_dir}/{chr}_telomere_no_telo_repeat_AC_reads_best.fasta'
    output_no_telo_repeat_TG_read_f_name=f'{repeat_trim_dir}/{chr}_telomere_no_telo_repeat_TG_reads_best.fasta'

    output_AC_read_f_name=f'{repeat_trim_dir}/{chr}_telomere_trimmed_AC_reads_best.fasta'
    output_TG_read_f_name=f'{repeat_trim_dir}/{chr}_telomere_trimmed_TG_reads_best.fasta'


    df_info_filtered_chr = df_info_filtered[df_info_filtered['anchor_name'] == f'{chr}_anchor']

    chr_wanted_read_list = df_info_filtered_chr['read_id'].tolist()

    print(f'Starting to make {chr}...')

    print(f'Grabbing reads in {chr}...')
    # Open the output file for writing

    no_telo_repeat_reads_AC = []
    no_telo_repeat_reads_TG = []
    trimmed_reads_AC = []
    trimmed_reads_TG = []
    outside_trimmed_reads = []

    for read_name in chr_wanted_read_list:

        sequence = fasta_index.fetch(read_name)
        read_sequence = sequence.strip('\n')
        row_index = df_info_filtered.loc[df_info_filtered['read_id'] == read_name].index[0]
        repeat_type = df_info_filtered.at[row_index,'Repeat_Type']

        if repeat_type == 'AC':
            # Skips read if not at least 10 bp's of telomere repeat
            try:
                trimming_results = trim_for_telomere_repeat(read_sequence, "AC")
            except AttributeError:
                #print(f'Skipping a AC read {AttributeError}')
                continue

            # No telomere repeat read if no telomere repeat found in first 75 basepairs
            if trimming_results == None:
                #print(f'For an AC read : NO REPEAT FOUND AT END')
                header = read_name + " AC " + chr + "\n"
                no_telo_repeat_reads_AC.append(f'{header}')
                no_telo_repeat_reads_AC.append(f'{read_sequence}\n')
                repeat_length = None
                length_to_trim_end_of_read = None
            else:
                header = read_name + " AC " + chr + "\n"
                trimmed_reads_AC.append(f'{header}')
                outside_trimmed_reads.append(f'{header}')

                trimmed_read = trimming_results[0]
                repeat_length = trimming_results[1]
                outside_trim = trimming_results[2]
                length_to_trim_end_of_read = trimming_results[3] # Based on the length to telo repeat

                trimmed_reads_AC.append(f'{trimmed_read}\n')
                outside_trimmed_reads.append(f'{outside_trim}\n')
                outside_trimmed_reads.append(f'{trimmed_read}\n')

        elif repeat_type == 'TG':
            # Skips read if not at least 10 bp's of telomere repeat
            try:
                trimming_results = trim_for_telomere_repeat(read_sequence, "TG")
            except AttributeError:
                #print(f'Skipping a TG read {AttributeError}')
                continue

            # No telomere repeat read if no telomere repeat found in first 75 basepairs
            if trimming_results == None:
                #print(f'For a TG read : NO REPEAT FOUND AT END')
                header = read_name + " TG " + chr + "\n"
                no_telo_repeat_reads_TG.append(f'{header}')
                no_telo_repeat_reads_TG.append(f'{read_sequence}\n')
                repeat_length = None
                length_to_trim_end_of_read = None
            else:
                header = read_name + " TG " + chr + "\n"
                trimmed_reads_TG.append(f'{header}')
                outside_trimmed_reads.append(f'{header}')

                trimmed_read = trimming_results[0]
                repeat_length = trimming_results[1]
                outside_trim = trimming_results[2]
                length_to_trim_end_of_read = trimming_results[3] # Based on the length to telo repeat

                trimmed_reads_TG.append(f'{trimmed_read}\n')
                outside_trimmed_reads.append(f'{outside_trim}\n')
                outside_trimmed_reads.append(f'{trimmed_read}\n')
        else:
            continue

        add_to_tsv_list.append([read_name, repeat_length, length_to_trim_end_of_read])


    with open(f'{output_AC_read_f_name}', "w") as f_out:
        f_out.writelines(trimmed_reads_AC)

    with open(f'{output_TG_read_f_name}', "w") as f_out:
        f_out.writelines(trimmed_reads_TG)

    with open(f'{output_no_telo_repeat_AC_read_f_name}', "w") as f_out:
        f_out.writelines(no_telo_repeat_reads_AC)

    with open(f'{output_no_telo_repeat_TG_read_f_name}', "w") as f_out:
        f_out.writelines(no_telo_repeat_reads_TG)


    both_no_telomere_repeat_telomere_reads_list = no_telo_repeat_reads_AC + no_telo_repeat_reads_TG
    running_both_no_telo_repeat_reads_list = both_no_telomere_repeat_telomere_reads_list + running_both_no_telo_repeat_reads_list

    both_to_telomere_reads_list = trimmed_reads_AC + trimmed_reads_TG
    running_both_to_telomere_reads_list = both_to_telomere_reads_list + running_both_to_telomere_reads_list

    running_outside_trimmed_reads_list = outside_trimmed_reads + running_outside_trimmed_reads_list


with open(f'{output_both_read_f_name}', "w") as f_out:
    f_out.writelines(running_both_no_telo_repeat_reads_list)

with open(f'{output_both_read_f_name}', "w") as f_out:
    f_out.writelines(running_both_to_telomere_reads_list)

with open(f'{output_both_trim_f_name}', "w") as f_out:
    f_out.writelines(running_outside_trimmed_reads_list)


df_trimmed_repeat = pd.DataFrame(add_to_tsv_list, columns =['read_id', 'repeat_length', 'length_to_trim_end_of_read'])

df_trimmed_info = pd.merge(df_trimmed_repeat, df_info_filtered, on='read_id', how='inner')


df_best_anchor_overhang = pd.read_csv(best_anchor_telo_f_name, sep='\t')

columns_to_include = ['read_id', 'read_length_past_anchor']

df_best_anchor_overhang_selected = df_best_anchor_overhang[columns_to_include].copy()


df_trimmed_and_overhang_info = pd.merge(df_trimmed_info, df_best_anchor_overhang_selected, on='read_id', how='inner')

# Trims overhang length if possible, but if trim_length = None then trimmed_length = read_length_past_anchor
def trim_read_lengths_past_anchor(read_delta, trim_length) -> int:
    if trim_length == None:
        trimmed_length = trim_length
    else:
        trimmed_length = read_delta - trim_length
    return trimmed_length

df_trimmed_and_overhang_info['trimmed_read_length_past_anchor'] = df_trimmed_and_overhang_info.apply(lambda row: trim_read_lengths_past_anchor(row["read_length_past_anchor"], row["length_to_trim_end_of_read"]), axis=1)


columns_in_order = ['read_id', 'anchor_name', 'Repeat_Type', 'repeat_length', 'Adapter_After_Telomere', 'Tag_After_Telomere',
                    'both_adapter_and_tag', 'either_adapter_or_tag', 'trimmed_read_length_past_anchor', 'length_to_trim_end_of_read']


df_final_trimmed_and_overhang_info = df_trimmed_and_overhang_info[columns_in_order].copy()


print(df_final_trimmed_and_overhang_info)

df_final_trimmed_and_overhang_info.to_csv(output_tsv_f_name, sep = '\t')

