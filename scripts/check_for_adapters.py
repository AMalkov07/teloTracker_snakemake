import pandas as pd
import os
import sys
import re

# Usage: check_for_adapters.py <base_name> <porechop_log> <output_dir>
base_name  = sys.argv[1]
log_file   = sys.argv[2]
output_dir = sys.argv[3]

out_f_name = os.path.join(output_dir, 'porechop', f'{base_name}_porechopped.tsv')

print("Starting check_for_adapters.py")

read_id_pattern             = r"\w+-\w+-\w+-\w+ (AC|TG)"
start_of_read_pattern       = "start:"
start_of_alignment_pattern  = "start alignments:"
end_of_read_pattern         = "end:"
end_of_alignment_pattern    = "end alignments:"
poly_A_alignment_pattern    = "AAAAAAAAAA"
poly_T_alignment_pattern    = "TTTTTTTTTT"

read_adapter_check_list = []
with open(log_file, 'r') as f_in:
    first_read = True
    read_strand = "none"
    for line in f_in:
        if re.search(read_id_pattern, line) is not None:
            if not first_read:
                adapter_details = [read_id, adapter_after_telomere_repeat, poly_tag_after_telomere_repeat, read_strand, telomere_sequence]
                read_adapter_check_list.append(adapter_details)
            first_read = False
            line = line.strip('\n')
            read_id     = line.split(' ')[0]
            read_strand = line.split(' ')[1]
            adapter_after_telomere_repeat = False
            continue

        elif start_of_read_pattern in line:
            if read_strand == 'AC':
                telomere_sequence = line.split(': ')[1].split('...')[0].split('\t')[0]
                try:
                    if telomere_sequence[0] in 'ATCG':
                        telomere_sequence = re.findall('[ACTG]+', telomere_sequence)[0]
                    else:
                        telomere_sequence = '-'.join(re.findall('[ACTG]+', telomere_sequence))
                    poly_tag_after_telomere_repeat = poly_T_alignment_pattern in telomere_sequence[:70]
                except IndexError:
                    telomere_sequence = 'N/A'
                    poly_tag_after_telomere_repeat = 'N/A'
            continue

        elif start_of_alignment_pattern in line:
            if read_strand == 'AC':
                adapter_after_telomere_repeat = True
            continue

        elif end_of_read_pattern in line:
            if read_strand == 'TG':
                telomere_sequence = line.split(': ')[1].split('...')[1].split('\t')[0]
                try:
                    if telomere_sequence[-1] in 'ATCG':
                        telomere_sequence = re.findall('[ACTG]+', telomere_sequence)[0]
                    else:
                        telomere_sequence = '-'.join(re.findall('[ACTG]+', telomere_sequence))
                    poly_tag_after_telomere_repeat = poly_A_alignment_pattern in telomere_sequence[-70:]
                except IndexError:
                    telomere_sequence = 'N/A'
                    poly_tag_after_telomere_repeat = 'N/A'
            continue

        elif end_of_alignment_pattern in line:
            if read_strand == 'TG':
                adapter_after_telomere_repeat = True
            continue

# Append the last read
adapter_details = [read_id, adapter_after_telomere_repeat, poly_tag_after_telomere_repeat, read_strand, telomere_sequence]
read_adapter_check_list.append(adapter_details)

df = pd.DataFrame(read_adapter_check_list,
                  columns=['read_id', 'Adapter_After_Telomere', 'Tag_After_Telomere', 'Repeat_Type', 'Telomere_Sequence'])
print(df)

df = df[df['Telomere_Sequence'] != 'N/A']
print(df)

df.to_csv(out_f_name, sep='\t')
