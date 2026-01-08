import pandas as pd
import os
import sys

print("Starting get_summary_stats_for_y_prime_repeatmasker.py")

file_names = sys.argv[1:]

df_combined = pd.DataFrame()
for file_name in file_names:
    if not os.path.exists(file_name):
        print(f"File {file_name} does not exist.")
        continue

    print(f"Processing file: {file_name}")

    strain_id = file_name.split("_")[1]
    day = file_name.split("_")[2]

    repeat = file_name.split("_")[3]
    if repeat.startswith("Pro"):
        repeat = 'original'
        
    print(f"strain_id: {strain_id}, day: {day}, repeat: {repeat}")
    
    # Read the CSV file
    df = pd.read_csv(file_name, sep="\t")

    df['strain_id'] = strain_id
    df['day'] = day
    df['repeat'] = repeat

    df_combined = pd.concat([df_combined, df], ignore_index=True)
    
    print(len(df_combined))

df_combined = df_combined[df_combined['sub_match'] == False]
df_combined = df_combined[df_combined['SW_score'] >= 10000]

# Reorder columns and drop unnecessary ones
df_combined = df_combined[['strain_id', 'day', 'repeat', 'read_id', 'chr_end', 'y_prime_id', 'y_prime_group',
                        'SW_score', 'match_start_on_read', 'match_end_on_read', 'leftover_on_read',
                        'Repeat_Type', 'telomere_side',
                        'match_start_on_y_prime', 'match_end_on_y_prime']]

print("Combined DataFrame:")
print(df_combined.head())

print(len(df_combined))



# Save the combined DataFrame to a new CSV file
df_combined.to_csv(f"wt_repeatmasker_results_to_check_for_recombinations_with_telomere_info.tsv",
                   sep="\t", index=False)

print(f"Saved combined results to wt_repeatmasker_results_to_check_for_recombinations_with_telomere_info.tsv")


output_file = "wt_recombinations_with_breaks_with_telomere_info.tsv"

with open(output_file, "w") as f:
    previous_read_id = None

    for idx, row in df_combined.iterrows():
        # Add the header only once
        if idx == 0:
            header = "\t".join(df_combined.columns)
            f.write(header + "\n")
        current_read_id = row['read_id']
        if current_read_id != previous_read_id and previous_read_id is not None:
            f.write("-----\t-----\t-----\t-----\t-----\t-----\t-----\t-----\t-----\t-----\t-----\t-----\t-----\t-----\t\n")  # Write separator line
        line = "\t".join(str(value) for value in row)
        f.write(line + "\n")
        previous_read_id = current_read_id

#   python create_summary_file_for_y_prime_repeatmasker_to_check_for_recombinations.py dorado_7172_day3*Pro*_good_end_y_repeatmasker.tsv dorado_7172_day4*Pro*_good_end_y_repeatmasker.tsv dorado_7302_day3*Pro*_good_end_y_repeatmasker.tsv dorado_7302_day4*Pro*_good_end_y_repeatmasker.tsv


