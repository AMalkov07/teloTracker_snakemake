import pandas as pd
import os
import sys

# Usage: filter_for_reads_with_anchors.py <base_name> <anchor_set> <output_dir> [foldback_ids.txt]
base_name  = sys.argv[1]
anchor_set = sys.argv[2]
output_dir = sys.argv[3]
# Optional: read IDs flagged by detect_foldback_reads.py (hairpin artifacts --
# the same anchor twice in opposite orientations). Excluding them here removes
# them from every downstream step, including scaffold selection in create_ref.
foldback_ids_file = sys.argv[4] if len(sys.argv) > 4 else None

anchor_blast_input  = f'{base_name}_blasted_{anchor_set}'
anchor_blast_f_name = os.path.join(output_dir, 'blast', f'{anchor_blast_input}.tsv')
all_reads_f_name    = os.path.join(output_dir, f'{base_name}.fasta')

output_f_name_all = os.path.join(output_dir, 'blast', f'all_matches_{anchor_blast_input}.tsv')
output_f_name_top = os.path.join(output_dir, 'blast', f'top_matches_{anchor_blast_input}.tsv')

print("Starting filter_for_reads_with_anchors.py")

def read_direction(match_start_on_anchor, match_end_on_anchor):
    return 'forward' if match_start_on_anchor < match_end_on_anchor else 'reverse'

def read_length_past_anchor_calc(total_read_length, match_start_on_read, match_end_on_read,
                                  match_start_on_anchor, match_end_on_anchor, l_end_chr, alignment_direction):
    if l_end_chr:
        if alignment_direction == 'forward':
            missing = match_start_on_anchor - 1
            return match_start_on_read - missing
        else:
            missing = match_end_on_anchor - 1
            return (total_read_length - match_end_on_read) - missing
    else:
        if alignment_direction == 'forward':
            missing = 5040 - match_end_on_anchor
            return (total_read_length - match_end_on_read) - missing
        else:
            missing = 5040 - match_start_on_anchor
            return match_start_on_read - missing

def wanted_sequence(l_end_chr, alignment_direction):
    if l_end_chr:
        return 'before_match_start_on_read' if alignment_direction == 'forward' else 'after_match_end_on_read'
    else:
        return 'after_match_end_on_read' if alignment_direction == 'forward' else 'before_match_start_on_read'

def repeat_type_of_read(l_end_chr, alignment_direction):
    if l_end_chr:
        return 'AC' if alignment_direction == 'forward' else 'TG'
    else:
        return 'TG' if alignment_direction == 'forward' else 'AC'

def match_length_calc(match_start_on_anchor, match_end_on_anchor):
    return abs(match_start_on_anchor - match_end_on_anchor)


chr_r_end_list = [f'chr{n}R_anchor' for n in range(1, 18)]
chr_l_end_list = [f'chr{n}L_anchor' for n in range(1, 18)]

df = pd.read_csv(anchor_blast_f_name, sep='\t')

# --- Exclude fold-back (hairpin) reads ---------------------------------------
# Detected by detect_foldback_reads.py, which re-BLASTs each anchored read
# against ONLY its assigned anchor without -min_raw_gapped_score (that flag
# filters on the preliminary gapped score and hides the degraded second copy).
# This can only remove reads -- it never adds or reassigns any.
if foldback_ids_file and os.path.exists(foldback_ids_file):
    with open(foldback_ids_file) as fh:
        foldback_ids = {ln.strip() for ln in fh if ln.strip()}
    if foldback_ids:
        before = df['read_id'].nunique()
        df = df[~df['read_id'].isin(foldback_ids)]
        after = df['read_id'].nunique()
        print(f'Fold-back filter: removed {before - after} of {before} reads '
              f'({len(foldback_ids)} flagged in {os.path.basename(foldback_ids_file)})')
    else:
        print('Fold-back filter: no reads flagged.')
elif foldback_ids_file:
    print(f'WARNING: fold-back ID file not found: {foldback_ids_file}')
else:
    print('Fold-back filter: not supplied (skipping).')
# ------------------------------------------------------------------------------

df_filtered = df[df['read_bp_used_for_match'] > df['total_anchor_length'] / 2]
df_filtered = df_filtered.sort_values(['read_id', 'bitscore', 'read_bp_used_for_match'], ascending=False)
df_unique   = df_filtered.drop_duplicates(subset=['read_id', 'anchor_name'], keep='first').copy()

df_unique['match_length']          = df_unique.apply(lambda r: match_length_calc(r['match_start_on_anchor'], r['match_end_on_anchor']), axis=1)
df_unique['l_end_chr']             = df_unique['anchor_name'].apply(lambda x: x in chr_l_end_list)
df_unique['alignment_direction']   = df_unique.apply(lambda r: read_direction(r['match_start_on_anchor'], r['match_end_on_anchor']), axis=1)
df_unique['repeat_type']           = df_unique.apply(lambda r: repeat_type_of_read(r['l_end_chr'], r['alignment_direction']), axis=1)
df_unique['read_length_past_anchor'] = df_unique.apply(
    lambda r: read_length_past_anchor_calc(r['total_read_length'], r['match_start_on_read'], r['match_end_on_read'],
                                           r['match_start_on_anchor'], r['match_end_on_anchor'], r['l_end_chr'], r['alignment_direction']), axis=1)
df_unique['wanted_section_of_read'] = df_unique.apply(lambda r: wanted_sequence(r['l_end_chr'], r['alignment_direction']), axis=1)

df_unique = df_unique.sort_values(['read_id', 'bitscore', 'match_length', 'read_bp_used_for_match'], ascending=False)
df_top    = df_unique.drop_duplicates(subset=['read_id'], keep=False).copy()

df_unique.to_csv(output_f_name_all, sep='\t', index=False)
df_top.to_csv(output_f_name_top, sep='\t', index=False)

print('df anchor_name counts:');           print(df['anchor_name'].value_counts())
print('\ndf_filtered anchor_name counts:'); print(df_filtered['anchor_name'].value_counts())
print('\ndf_unique anchor_name counts:');   print(df_unique['anchor_name'].value_counts())
print('\ndf_top:');                         print(df_top); print(df_top['anchor_name'].value_counts())

print(f'Making file and Indexing Reads...')
os.system(f'samtools faidx {all_reads_f_name}')
