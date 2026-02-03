import pandas as pd

df1 = pd.read_csv('output_combined_y_prime_recombination.tsv', sep='\t')
df2 = pd.read_csv('all_spacer_switches-all_changes.tsv', sep='\t')

print((len(df1), len(df2)))

pd_combined = pd.merge(df1, df2, how='left', on=['read_id'])

print((len(pd_combined), pd_combined.columns))

pd_combined.to_csv('output_combined_y_prime_recombination_with_spacer_switches.tsv', sep='\t', index=False)
print("Output saved to 'output_combined_y_prime_recombination_with_spacer_switches.tsv'")

# Drop specific columns that are not needed
columns_to_drop = ['strain', 'genotype', 'day', 'repeat']
pd_combined.drop(columns=columns_to_drop, inplace=True)


pd_combined['strain'] = pd_combined['source_file'].str.split('_').str[1]
pd_combined['day'] = pd_combined['source_file'].str.split('_').str[2]
pd_combined['repeat'] = pd_combined['source_file'].str.split('_').str[3]
pd_combined['tagging'] = pd_combined['source_file'].str.split('_').str[4]



# Filter out columns that are not needed
columns_to_keep = [
    'read_id', 'strain', 'day', 'repeat', 'tagging',
    'chr_end', 'Repeat_Type', 'repeat_length', 'reference_y_primes',
    'y_prime_probe_count', 'better_repeatmasker_y_prime_count',
    'y_prime_recombination_status', 'y_prime_recombination_events',
    'anchor_and_spacer_switch_to_chr_end_pair', 'has_x_end', 'has_its_in_donor',
]

df = pd_combined[columns_to_keep]

# Change 1st Y' Change to Pre-Y' Recombination
df['y_prime_recombination_status'] = df['y_prime_recombination_status'].replace({
    "1st Y' Change": "Pre-Y' Recombination"
})


df.to_csv('output_combined_y_prime_recombination_with_spacer_switches_filtered.tsv', sep='\t', index=False)

