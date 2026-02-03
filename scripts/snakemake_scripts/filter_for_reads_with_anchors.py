import pandas as pd
import os
import matplotlib.pyplot as plt
import sys

#base_name=f'{sys.argv[1]}'
#anchor_set=f'{sys.argv[2]}'    # "new_all_best_anchors"  # f'{sys.argv[2]}' # Default = new_all_best_anchors
anchor_blast_f_name=f'{sys.argv[1]}'
output_f_name_all=f'{sys.argv[2]}'
output_f_name_top=f'{sys.argv[3]}'


#anchor_blast_input=f'{base_name}_blasted_{anchor_set}'
#anchor_blast_f_name=f'{base_name}/{anchor_blast_input}.tsv'
#all_reads_f_name=f'{base_name}/{anchor_blast_input.split("_blasted")[0]}.fasta'

#output_f_name_all=f'{base_name}/all_matches_{anchor_blast_input}.tsv'
#output_f_name_top=f'{base_name}/top_matches_{anchor_blast_input}.tsv'

print("Starting filter_for_reads_with_anchors.py")
print(f"Input: {anchor_blast_f_name}")

# Get the read direction in relation to the reference
def read_direction(match_start_on_anchor, match_end_on_anchor):
    if match_start_on_anchor < match_end_on_anchor:
        return 'forward'
    else:
        return 'reverse'

# Get the length past the end of the anchor of the read
def read_length_past_anchor_calc(total_read_length, match_start_on_read, match_end_on_read, match_start_on_anchor, match_end_on_anchor, l_end_chr, alignment_direction):
    # Length of read length past anchor = (Length of read past anchor) - (any missing distance missing from end of anchor in the blast match)
    if l_end_chr == True:
        if alignment_direction == 'forward':
            missing_distance_from_blast_match = match_start_on_anchor - 1   # Will = 0 if blast match span the end of the anchor towards the telomere
            read_length_past_anchor = match_start_on_read - missing_distance_from_blast_match
        else:
            missing_distance_from_blast_match = match_end_on_anchor - 1   # Will = 0 if blast match span the end of the anchor towards the telomere
            read_length_past_anchor = (total_read_length - match_end_on_read) - missing_distance_from_blast_match
    else:
        if alignment_direction == 'forward':
            missing_distance_from_blast_match = 5040 - match_end_on_anchor   # Will = 0 if blast match span the end of the anchor towards the telomere
            read_length_past_anchor = (total_read_length - match_end_on_read) - missing_distance_from_blast_match
        else:
            missing_distance_from_blast_match = 5040 - match_start_on_anchor   # Will = 0 if blast match span the end of the anchor towards the telomere
            read_length_past_anchor = match_start_on_read - missing_distance_from_blast_match
    return read_length_past_anchor

# Label where in the read we will want to look (for overhang)
def wanted_sequence(l_end_chr, alignment_direction):
    if l_end_chr == True:
        if alignment_direction == 'forward':
            wanted = 'before_match_start_on_read'
        else:
            wanted = 'after_match_end_on_read'
    else:
        if alignment_direction == 'forward':
            wanted = 'after_match_end_on_read'
        else:
            wanted = 'before_match_start_on_read'
    return wanted

# Label the repeat type of the read (AC/TG)
def repeat_type_of_read(l_end_chr, alignment_direction):
    if l_end_chr == True:
        if alignment_direction == 'forward':
            repeat_type = 'AC'
        else:
            repeat_type = 'TG'
    else:
        if alignment_direction == 'forward':
            repeat_type = 'TG'
        else:
            repeat_type = 'AC'
    return repeat_type


# Get each read's length of the anchor it matches
def match_length_calc(match_start_on_anchor, match_end_on_anchor):
    match_length = abs(match_start_on_anchor - match_end_on_anchor)
    return match_length


chr_r_end_list = ['chr10R_anchor',
'chr11R_anchor',
'chr12R_anchor',
'chr13R_anchor',
'chr14R_anchor',
'chr15R_anchor',
'chr16R_anchor',
'chr17R_anchor',
'chr1R_anchor',
'chr2R_anchor',
'chr3R_anchor',
'chr4R_anchor',
'chr5R_anchor',
'chr6R_anchor',
'chr7R_anchor',
'chr8R_anchor',
'chr9R_anchor'
]

chr_l_end_list = ['chr10L_anchor',
'chr11L_anchor',
'chr12L_anchor',
'chr13L_anchor',
'chr14L_anchor',
'chr15L_anchor',
'chr16L_anchor',
'chr17L_anchor',
'chr1L_anchor',
'chr2L_anchor',
'chr3L_anchor',
'chr4L_anchor',
'chr5L_anchor',
'chr6L_anchor',
'chr7L_anchor',
'chr8L_anchor',
'chr9L_anchor'
]

# Load TSV file into a Pandas DataFrame with custom column headers
df = pd.read_csv(anchor_blast_f_name, sep='\t')
# Names from blast = ['qseqid', 'qlen', 'length', 'qstart', 'qend', 'sseqid', 'slen', 'sstart', 'send', 'pident', 'bitscore', 'evalue'])
# ['read_id', 'total_read_length', 'read_bp_used_for_match', 'match_start_on_read', 'match_end_on_read', 'anchor_name', 'total_anchor_length', 'match_start_on_anchor', 'match_end_on_anchor', 'pident', 'bitscore', 'evalue'])

# Processing the table

# Filter for only reads with more than 2500 bp matching the anchor
df_filtered = df[df['read_bp_used_for_match'] > df['total_anchor_length']/2]

# Sort the values by best matches (match_length is consistent with bitscore)
df_filtered = df_filtered.sort_values(['read_id', 'bitscore', 'read_bp_used_for_match'], ascending=False)

# Keep only one row for each unique entry in anchor_name (chromosome) per read_id (read) e.g. One mapping for a given chr_anchor per read
df_unique = df_filtered.drop_duplicates(subset=['read_id', 'anchor_name'], keep="first").copy()

# Get each read's length of the anchor it matches
df_unique["match_length"] = df_unique.apply(lambda row: match_length_calc(row["match_start_on_anchor"], row["match_end_on_anchor"]), axis=1)

# Get the read direction in relation to the reference
df_unique['l_end_chr'] = df_unique['anchor_name'].apply(lambda x : x in chr_l_end_list)

# Get each read's alignment direction (direction in relation to ref)
df_unique["alignment_direction"] = df_unique.apply(lambda row: read_direction(row["match_start_on_anchor"], row["match_end_on_anchor"]), axis=1)

# Get each read's Repeat Type AC/TG
df_unique["repeat_type"] = df_unique.apply(lambda row: repeat_type_of_read(row["l_end_chr"], row["alignment_direction"]), axis=1)

# Get each read's overhang length
df_unique["read_length_past_anchor"] = df_unique.apply(lambda row: read_length_past_anchor_calc(row["total_read_length"], row["match_start_on_read"], row["match_end_on_read"], row["match_start_on_anchor"], row["match_end_on_anchor"], row["l_end_chr"], row["alignment_direction"]), axis=1)

# Get each read's wanted section
df_unique["wanted_section_of_read"] = df_unique.apply(lambda row: wanted_sequence(row["l_end_chr"], row["alignment_direction"]), axis=1)

# Sort the values by best matches (match_length is consistent with bitscore)
df_unique = df_unique.sort_values(['read_id', 'bitscore', 'match_length', 'read_bp_used_for_match'], ascending=False)

# Gets best match (top remaining) for the anchor
df_top = df_unique.drop_duplicates(subset=['read_id'], keep=False).copy()

#print(df_top)


# Make files
df_unique.to_csv(output_f_name_all, sep="\t",index=False)
df_top.to_csv(output_f_name_top, sep="\t",index=False)

print('df')
print(df['anchor_name'].value_counts())
print('')

print('df_filtered')
print(df_filtered['anchor_name'].value_counts())
print('')

print('df_unique')
print(df_unique['anchor_name'].value_counts())
print('')

print('df_top')
print(df_top)
print(df_top['anchor_name'].value_counts())
print('')

# Load the FAIDX file
#print(f'Making file and Indexing Reads...')
#os.system(f'samtools faidx {all_reads_f_name}')


