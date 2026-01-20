"""
Shared utilities for parsing and processing RepeatMasker output files.

This module provides common functions used by:
- make_y_prime_repeatmasker_tsv.py
- make_x_element_ends_pairs_repeatmasker_tsv.py
- make_spacer_pairs_repeatmasker_tsv.py
"""

import pandas as pd
import os


def correct_repeatmasker_spacing(input_file):
    """
    Correct spacing issues in RepeatMasker output files.

    RepeatMasker adds leading spaces when SW_score < 10,000,
    which causes parsing issues. This function removes those spaces
    and adds a '-' for missing sub_match values.

    Args:
        input_file (str): Path to the original RepeatMasker .ssv file

    Returns:
        str: Path to the corrected .ssv file
    """
    corrected_file = input_file.replace('results.ssv', 'results_corrected.ssv')

    # Skip if already corrected
    if os.path.isfile(corrected_file):
        return corrected_file

    with open(input_file, "r") as original:
        corrected_data = []
        for line_number, line in enumerate(original.readlines()):
            if line_number == 0:
                # Keep header as is
                corrected_data.append(line)
            else:
                line = line.strip()
                # Add '-' if sub_match column is missing
                if line and line[-1] != "*":
                    line = f'{line} -'
                corrected_data.append(line)

    with open(corrected_file, "w") as corrected:
        for fixed_line in corrected_data:
            corrected.write(f'{fixed_line}\n')

    return corrected_file


def load_and_aggregate_repeatmasker_results(results_dir, file_pattern_end, chr_end_extractor=None):
    """
    Load and aggregate multiple RepeatMasker result files.

    Args:
        results_dir (str): Directory containing RepeatMasker results
        file_pattern_end (str): Ending pattern for files to include (e.g., 'repeatmasker_results.ssv')
        chr_end_extractor (callable, optional): Function to extract chr_end from filename.
                                                If None, no chr_end column is added.

    Returns:
        pd.DataFrame: Aggregated RepeatMasker results with corrected spacing
    """
    results_files = os.listdir(results_dir)
    matching_files = [
        os.path.join(results_dir, f)
        for f in results_files
        if f.endswith(file_pattern_end)
    ]

    print(f'Found {len(matching_files)} RepeatMasker result files')

    df_all_results = pd.DataFrame()

    for results_file in matching_files:
        try:
            # Correct spacing issues
            corrected_file = correct_repeatmasker_spacing(results_file)

            # Read the corrected file
            df_single = pd.read_csv(corrected_file, sep=r"\s+")

            # Extract chr_end if extractor function provided
            if chr_end_extractor:
                chr_end = chr_end_extractor(corrected_file)
                if chr_end:
                    df_single["chr_end"] = chr_end

            # Remove empty columns
            df_single = df_single.dropna(axis=1, how='all')

            # Concatenate to main dataframe
            df_all_results = pd.concat([df_all_results, df_single])

        except Exception as e:
            print(f'Error processing {results_file}: {e}')
            continue

    # Convert sub_match column: "*" -> True, "-" -> False
    if 'sub_match' in df_all_results.columns:
        df_all_results.loc[df_all_results['sub_match'] == '*', 'sub_match'] = True
        df_all_results.loc[df_all_results['sub_match'] == '-', 'sub_match'] = False

    return df_all_results


def filter_good_reads(df_all_results, y_prime_probe_file, min_sw_score=500, max_divergence=2.0):
    """
    Filter RepeatMasker results for high-quality reads.

    Args:
        df_all_results (pd.DataFrame): All RepeatMasker results
        y_prime_probe_file (str): Path to Y prime probe TSV file
        min_sw_score (int): Minimum Smith-Waterman score (default: 500)
        max_divergence (float): Maximum divergence percent (default: 2.0)

    Returns:
        tuple: (df_good, df_filter) - Filtered results and filter dataframe
    """
    # Filter by quality metrics
    df_good = df_all_results[df_all_results['sub_match'] == False].copy()
    df_good = df_good[df_good['SW_score'] >= min_sw_score]
    df_good = df_good[df_good['divergence_percent'] <= max_divergence]

    # Filter by good telomere ends from y_prime_probe file
    df_filter = pd.read_csv(y_prime_probe_file, sep='\t')
    df_filter = df_filter.dropna(subset=["repeat_length"])
    df_filter = df_filter[df_filter['Adapter_After_Telomere'] == True]

    good_reads = df_filter['read_id'].to_list()
    df_good['good_read'] = df_good['read_id'].apply(lambda x: x in good_reads)

    print(f"Good read filtering: {df_good['good_read'].value_counts()}")

    df_good = df_good[df_good['good_read'] == True]

    return df_good, df_filter


def filter_gained_y_prime_reads(df_good, df_filter):
    """
    Filter for reads that gained Y' elements (delta_y_prime_sign == '+').

    Args:
        df_good (pd.DataFrame): Good quality RepeatMasker results
        df_filter (pd.DataFrame): Filtered Y prime probe data

    Returns:
        pd.DataFrame: Results for reads with gained Y' elements
    """
    # Filter for gained Y' elements
    df_filter_gained = df_filter[df_filter['delta_y_prime_sign'] == '+']
    reads_with_gain = df_filter_gained['read_id'].to_list()

    print(f'Reads with gain of Y prime: {len(reads_with_gain)}')

    df_good['gained_y'] = df_good['read_id'].apply(lambda x: x in reads_with_gain)

    print(f"Gained Y filtering: {df_good['gained_y'].value_counts()}")

    df_gained = df_good[df_good['gained_y'] == True].copy()

    return df_gained
