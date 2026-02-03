import pandas as pd
import glob

# Define a pattern that matches all your target files
file_pattern = "dorado*y_prime_recombination.tsv"

# Get all matching files in the current directory
file_list = glob.glob(file_pattern)

# Read and concatenate all files
df_all = pd.concat(
    [pd.read_csv(f, sep='\t').assign(source_file=f) for f in file_list],
    ignore_index=True
)

# Optional: Save the combined dataframe
df_all.to_csv("output_combined_y_prime_recombination.tsv", sep='\t', index=False)

print(f"Combined {len(file_list)} files into one DataFrame with {len(df_all)} rows.")

