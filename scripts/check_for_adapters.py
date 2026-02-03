import pandas as pd
import os
import sys
import re


base_name = f'{sys.argv[1]}' # base='guppy-6991_day0_PromethION_no_tag_no_rejection-trim-and-detect-pass'
log_file = f'{sys.argv[2]}' # log_file=f'{all_chrs_split_to_telomere}_porechopped.log'
out_f_name = f'results/{base_name}/{base_name}_porechopped.tsv'

print("Starting check_for_adapters.py")

read_id_pattern = r"\w+-\w+-\w+-\w+ (AC|TG)"
start_of_read_pattern = "start:"
start_of_alignment_pattern = "start alignments:"
end_of_read_pattern = "end:"
end_of_alignment_pattern = "end alignments:"

poly_A_alignment_pattern = "AAAAAAAAAA"
poly_T_alignment_pattern = "TTTTTTTTTT"


read_adapter_check_list = []
with open(log_file, 'r') as f_in:
    first_read = True
    read_strand = "none"
    for line in f_in:
        # If line matches the read id pattern -- start of read 
        if re.search(read_id_pattern, line) != None:
            if first_read == False: # Only flips to False once (when read seciton of logs starts)
                adapter_details = [read_id, adapter_after_telomere_repeat, poly_tag_after_telomere_repeat, read_strand, telomere_sequence]
                read_adapter_check_list.append(adapter_details)            

            first_read = False # Only flips to Flase once (when read seciton of logs starts)
            line = line.strip("\n")
            read_id = line.split(" ")[0]
            read_strand = line.split(" ")[1]     
            adapter_after_telomere_repeat = False
            continue
        elif start_of_read_pattern in line:
            if read_strand == "AC":
                telomere_sequence = line.split(": ")[1]
                telomere_sequence = telomere_sequence.split("...")[0]
                telomere_sequence = telomere_sequence.split("\t")[0]
                try:
                    if telomere_sequence[0] in 'ATCG':
                        telomere_sequence = re.findall("[ACTG]+",telomere_sequence)[0]
                    else:
                        telomere_sequence = list(re.findall("[ACTG]+",telomere_sequence))
                        telomere_sequence = "-".join(telomere_sequence)

                    if poly_T_alignment_pattern in telomere_sequence[:70]:
                        poly_tag_after_telomere_repeat = True
                    else:
                        poly_tag_after_telomere_repeat = False
                # Handles index error where read is split to 0 bp's and shows up as 0 in porechop logs - removed later
                except IndexError:
                    telomere_sequence = 'N/A'
                    poly_tag_after_telomere_repeat = 'N/A'    
            continue
        elif start_of_alignment_pattern in line:
            if read_strand == "AC":
                adapter_after_telomere_repeat = True
            continue
        elif end_of_read_pattern in line:
            if read_strand == "TG":            
                telomere_sequence = line.split(": ")[1]
                telomere_sequence = telomere_sequence.split("...")[1]
                telomere_sequence = telomere_sequence.split("\t")[0]
                try:
                    if telomere_sequence[-1] in 'ATCG':
                        telomere_sequence = re.findall("[ACTG]+",telomere_sequence)[0]
                    else:
                        telomere_sequence = list(re.findall("[ACTG]+",telomere_sequence))
                        telomere_sequence = "-".join(telomere_sequence)


                    if poly_A_alignment_pattern in telomere_sequence[-70:]:
                        poly_tag_after_telomere_repeat = True
                    else:
                        poly_tag_after_telomere_repeat = False
                # Handles index error where read is split to 0 bp's and shows up as 0 in porechop logs - removed later
                except IndexError:
                    telomere_sequence = 'N/A'
                    poly_tag_after_telomere_repeat = 'N/A'       
            continue
        elif end_of_alignment_pattern in line:
            if read_strand == "TG":
                adapter_after_telomere_repeat = True
            continue
        else:
            continue

# Append the last read
adapter_details = [read_id, adapter_after_telomere_repeat, poly_tag_after_telomere_repeat, read_strand, telomere_sequence]
read_adapter_check_list.append(adapter_details)




df = pd.DataFrame(read_adapter_check_list, columns =['read_id', 'Adapter_After_Telomere', 'Tag_After_Telomere', 'Repeat_Type', 'Telomere_Sequence'])
print(df)


# Removing the telomeres with 0 bp's from earlier
df = df[df["Telomere_Sequence"] != 'N/A']

print(df)

df.to_csv(out_f_name, sep = '\t')
