"""Y' probe count per anchored read, compared with the Y' copies expected at the read's chromosome
end, plus the Y' plots and statistics.

Usage: y_prime_analysis.py <base_name> <output_dir> [--expected FILE] [--expected-label TEXT]
                           [--no-adapter-filter]

--expected is a label BED or a counts table giving the expected Y' copies per end (see
yprime_reference.py; run_pipeline.py picks it from yprime_reference). With it, each read gets its
delta (probe count minus expected) and the delta plots are drawn. Without it the reference and
delta columns are left empty, and the plots show absolute Y' counts per end and counts relative to
the sample's own per-end median instead.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from telomere_filters import qualifying_mask
from yprime_reference import load_expected_counts


def read_direction(match_start_on_anchor, match_end_on_anchor):
    return 'forward' if match_start_on_anchor < match_end_on_anchor else 'reverse'

def repeat_type_of_read(l_end_chr, alignment_direction):
    if l_end_chr:
        return 'AC' if alignment_direction == 'forward' else 'TG'
    else:
        return 'TG' if alignment_direction == 'forward' else 'AC'

def match_length_calc(match_start_on_anchor, match_end_on_anchor):
    return abs(match_start_on_anchor - match_end_on_anchor)


def assign_delta_y_prime_sign(y_primes_relative_to_ref):
    if y_primes_relative_to_ref > 0:   return '+'
    elif y_primes_relative_to_ref < 0: return '-'
    else:                               return 'same'

def assign_reference_y_prime_end_status(reference_y_primes):
    if reference_y_primes == 0:    return 'Reference Y Primes = 0'
    if reference_y_primes == 1:    return 'Reference Y Primes = 1'
    elif reference_y_primes >= 2:  return 'Reference Y Primes = 2+'

def assign_read_y_prime_end_status(y_prime_probe_count):
    if y_prime_probe_count == 0:    return 'Read Y Primes = 0'
    if y_prime_probe_count == 1:    return 'Read Y Primes = 1'
    elif y_prime_probe_count >= 2:  return 'Read Y Primes = 2+'


def expected_y_prime_column(anchor_names, expected):
    """Expected Y' copies for each read's end ('chr4R_anchor' -> expected['chr4R']). Plain ints
    when every end in the data is known; a nullable Int64 column (blank in the TSV) otherwise."""
    ends = anchor_names.str.replace('_anchor', '', regex=False)
    if ends.isin(list(expected)).all():
        return ends.map(expected).astype('int64')
    return ends.map(expected).astype('Int64')


# ---------------------------------------------------------------------------
# Y prime plots
# ---------------------------------------------------------------------------

def chr_end_violin_strip_plot(dataframe, y, title, ylabel, save_file_name, plot_scale=(22, 9)):
    """Violin + strip plot of one per-read value for each chromosome end, one colour per end."""

    sns.set(style="whitegrid")

    chr_order = [f'{n}{arm}' for n in range(1, 17) for arm in ('L', 'R')]
    chr_list  = [c for c in chr_order if c in dataframe['chr_end'].unique()]

    if not chr_list:
        print(f'Warning: no chr_end values found in data, skipping {os.path.basename(save_file_name)}')
        return

    full_color_palette = ['#641E16', '#512E5F', '#4A235A', '#154360', '#2E86C1', '#0E6251', '#0B5345', '#145A32',
                          '#7D6608', '#7E5109', '#D35400', '#626567', '#4D5656', '#212F3C', '#17202A', '#F93409',
                          '#FFED24', '#CD6155', '#D7BDE2', '#BB8FCE', '#A9CCE3', '#AED6F1', '#A3E4D7', '#45B39D',
                          '#229954', '#F7DC6F', '#F5B041', '#E59866', '#D7DBDD', '#95A5A6', '#5D6D7E', '#ABB2B9',
                          '#FFC2B4', '#FFFEB4']
    chr_color_dict = dict(zip(chr_order, full_color_palette))
    chr_color_code = [chr_color_dict[c] for c in chr_list]
    sns.set_palette(chr_color_code)

    fig, ax = plt.subplots(figsize=plot_scale)
    sns.violinplot(x='chr_end', hue='chr_end', order=chr_list, y=y, data=dataframe, gridsize=1000, cut=0, palette=chr_color_code,legend=False, ax=ax)
    sns.stripplot(x='chr_end', order=chr_list, y=y, data=dataframe, linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="#FF2400", ax=ax)

    ax.set_title(title, fontweight="bold", fontsize=20, color='k', pad=15)
    plt.xlabel("Chromosome End", fontweight="bold", fontsize=20)
    plt.ylabel(ylabel, fontweight="bold", fontsize=30)
    plt.xticks(fontweight="bold", fontsize=13)
    plt.yticks(fontweight="bold", fontsize=20)

    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=300, format="png")
    plt.close()


def y_prime_count_violin_strip_plot(dataframe, base_name, figures_dir):
    total_reads = len(dataframe)
    chr_end_violin_strip_plot(
        dataframe, 'y_primes_relative_to_ref',
        f'{base_name} Delta Y Primes (N = {total_reads} Reads, Ave. Y Prime to ref = {dataframe["y_primes_relative_to_ref"].mean():.3f})',
        "Delta Y Primes", os.path.join(figures_dir, f'{base_name}_delta_y_primes.png'))


def delta_sign_violin_strip_plot(dataframe, base_name, figures_dir, plot_scale=(16, 9)):

    sns.set(style="whitegrid")

    fig, ax = plt.subplots(figsize=plot_scale)
    sns.violinplot(x='delta_y_prime_sign', hue='delta_y_prime_sign', hue_order=["-", "same", "+"], order=["-", "same", "+"], y='y_primes_relative_to_ref',
                   data=dataframe, gridsize=1000, cut=0, palette=['tab:red', 'tab:blue', 'tab:green'], legend=False, ax=ax)
    ax.set(ylim=(-10, 20))
    sns.stripplot(x='delta_y_prime_sign', order=["-", "same", "+"], y='y_primes_relative_to_ref', data=dataframe,
                  linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="tab:gray", ax=ax)

    total_reads    = len(dataframe)
    total_same     = len(dataframe[dataframe['delta_y_prime_sign'] == "same"])
    total_positive = len(dataframe[dataframe['delta_y_prime_sign'] == "+"])
    total_negative = len(dataframe[dataframe['delta_y_prime_sign'] == "-"])

    ax.set_title(f'{base_name} Delta Y Primes (N = {total_reads} Reads, -:{total_negative}, =:{total_same}, & +:{total_positive})',
                 fontweight="bold", fontsize=20, color='k', pad=15)
    plt.xlabel("Sign of Delta Y Primes", fontweight="bold", fontsize=20)
    plt.ylabel("Number of Delta Y Primes", fontweight="bold", fontsize=30)
    plt.xticks(fontweight="bold", fontsize=13)
    plt.yticks(fontweight="bold", fontsize=20)

    save_file_name = os.path.join(figures_dir, f'{base_name}_sign_delta_y_primes.png')
    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=300, format="png")
    plt.close()


def absolute_count_plots(dataframe, base_name, figures_dir, with_median_view):
    """Y' copies per read for each end and, without a reference, the same relative to the sample's
    own per-end median, which shows spread and outliers the way the delta plot does."""
    df = dataframe.dropna(subset=['y_prime_probe_count', 'chr_end'])
    n = len(df)
    chr_end_violin_strip_plot(
        df, 'y_prime_probe_count',
        f'{base_name} Y Primes per Read (N = {n} Reads, Ave. = {df["y_prime_probe_count"].mean():.3f})',
        "Y Primes per Read", os.path.join(figures_dir, f'{base_name}_y_prime_counts.png'))
    if with_median_view:
        df = df.assign(y_primes_vs_median=df['y_prime_probe_count']
                       - df.groupby('chr_end')['y_prime_probe_count'].transform('median'))
        chr_end_violin_strip_plot(
            df, 'y_primes_vs_median',
            f'{base_name} Y Primes vs the Sample Median of Each End (N = {n} Reads)',
            "Y Primes vs Median", os.path.join(figures_dir, f'{base_name}_y_primes_vs_sample_median.png'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('base_name')
    ap.add_argument('output_dir')
    ap.add_argument('--expected', default='',
                    help="label BED or counts table of expected Y' copies per end (empty: none)")
    ap.add_argument('--expected-label', default='', help='where the expected counts came from (recorded)')
    ap.add_argument('--no-adapter-filter', action='store_true',
                    help='qualifying reads need no adapter after the telomere (the adapter check was off)')
    args = ap.parse_args()
    base_name, output_dir = args.base_name, args.output_dir

    print("Starting y_prime_analysis.py")

    y_prime_dir    = os.path.join(output_dir, 'blast', 'y_prime_probe')
    figures_dir    = os.path.join(output_dir, 'graphs', 'y_prime_figures')
    stats_dir      = os.path.join(output_dir, 'graphs', 'stats_for_y_primes')

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(stats_dir,   exist_ok=True)

    telomere_repeat_results_input = os.path.join(output_dir, f'{base_name}_post_telo_trimming.tsv')
    output_tsv_f_name             = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')

    expected = load_expected_counts(args.expected) if args.expected else None

    # ---- Load all per-chromosome Y' probe blast results ----
    probe_files = [f for f in os.listdir(y_prime_dir) if f.endswith('blasted_probe.tsv')]
    df = pd.DataFrame()
    for f in probe_files:
        df_temp = pd.read_csv(os.path.join(y_prime_dir, f), sep='\t')
        df = pd.concat([df, df_temp])

    # Save combined probe results
    combined_probe_out = os.path.join(y_prime_dir, f'all_{base_name}_probe_matches.tsv')
    df.to_csv(combined_probe_out, sep='\t')

    column_type_dict = {
        'read_id': str, 'total_read_length': int, 'read_bp_used_for_match': int,
        'match_start_on_read': int, 'match_end_on_read': int, 'anchor_name': str,
        'total_anchor_length': int, 'match_start_on_anchor': int, 'match_end_on_anchor': int,
        'pident': float, 'bitscore': float, 'evalue': float
    }
    df = df.astype(column_type_dict)

    df_filtered = df[df['read_bp_used_for_match'] > df['total_anchor_length'] / 2]
    df_filtered = df_filtered.sort_values(['read_id', 'match_start_on_read', 'bitscore'], ascending=False)
    df_filtered['match_length'] = df_filtered.apply(
        lambda r: match_length_calc(r['match_start_on_anchor'], r['match_end_on_anchor']), axis=1)

    chr_l_end_list = [f'chr{n}L_anchor' for n in range(1, 17)]
    df_filtered['l_end_chr']          = df_filtered['anchor_name'].apply(lambda x: x in chr_l_end_list)
    df_filtered['alignment_direction'] = df_filtered.apply(
        lambda r: read_direction(r['match_start_on_anchor'], r['match_end_on_anchor']), axis=1)
    df_filtered['repeat_type']        = df_filtered.apply(
        lambda r: repeat_type_of_read(r['l_end_chr'], r['alignment_direction']), axis=1)

    read_y_prime_counts      = df_filtered['read_id'].value_counts().reset_index()
    read_y_prime_counts.columns = ['read_id', 'y_prime_probe_count']

    df_repeat_results = pd.read_csv(telomere_repeat_results_input, sep='\t', index_col=0)
    df_repeat_results['read_id'] = df_repeat_results['read_id'].astype(str)

    df_combined = pd.merge(df_repeat_results, read_y_prime_counts, on='read_id', how='outer')

    df_combined['y_prime_probe_count'] = df_combined['y_prime_probe_count'].fillna(0)
    df_combined = df_combined.dropna(subset=['anchor_name'])

    if expected is not None:
        df_combined['reference_y_primes']       = expected_y_prime_column(df_combined['anchor_name'], expected)
        df_combined['y_primes_relative_to_ref'] = (df_combined['y_prime_probe_count']
                                                   - df_combined['reference_y_primes']).astype('float64')
        df_combined['delta_y_prime_sign']           = df_combined['y_primes_relative_to_ref'].apply(
            lambda d: assign_delta_y_prime_sign(d) if pd.notna(d) else '')
        df_combined['reference_y_prime_end_status'] = df_combined['reference_y_primes'].apply(
            lambda r: assign_reference_y_prime_end_status(r) if pd.notna(r) else '')
    else:
        # No reference: nothing to compare with, so the reference and delta columns stay empty.
        df_combined['reference_y_primes']           = float('nan')
        df_combined['y_primes_relative_to_ref']     = float('nan')
        df_combined['delta_y_prime_sign']           = ''
        df_combined['reference_y_prime_end_status'] = ''
    df_combined['read_y_prime_end_status']      = df_combined['y_prime_probe_count'].apply(assign_read_y_prime_end_status)
    df_combined['y_prime_type_delta_y_prime_sign'] = (
        (df_combined['delta_y_prime_sign'] + ' ' + df_combined['reference_y_prime_end_status'])
        .where(df_combined['delta_y_prime_sign'] != '', ''))

    print(df_combined)

    chr_end_list = ['1L','1R','2L','2R','3L','3R','4L','4R','5L','5R','6L','6R','7L','7R','8L','8R',
                    '9L','9R','10L','10R','11L','11R','12L','12R','13L','13R','14L','14R','15L','15R','16L','16R']
    chr_end_sorter = {c: i for i, c in enumerate(chr_end_list)}

    df_combined['chr_end']      = df_combined['anchor_name'].apply(lambda a: a.strip('chr _anchor'))
    df_combined['chr_end_rank'] = df_combined['chr_end'].map(chr_end_sorter)
    df_combined.sort_values('chr_end_rank', ascending=True, inplace=True)

    df_combined.to_csv(output_tsv_f_name, sep='\t', index=False)
    print(f'Written: {output_tsv_f_name}')

    # Where the expected counts came from, and the counts themselves (shown in the sample report).
    with open(os.path.join(stats_dir, f'{base_name}_yprime_reference.tsv'), 'w') as f:
        f.write(f'# source: {args.expected_label or (args.expected or "none")}\n')
        f.write(f'# file: {args.expected}\n')
        f.write('chr_end\texpected_y_primes\n')
        for end, n in (expected or {}).items():
            f.write(f'{end}\t{n}\n')

    # ---- Good telomere reads filter (the qualifying reads; see telomere_filters.py) ----
    df_good_reads = df_combined.loc[qualifying_mask(df_combined, require_adapter=not args.no_adapter_filter)]
    print('Reads in df_good_reads:')
    print(df_good_reads)
    print(df_good_reads['y_prime_probe_count'].value_counts())
    print(df_good_reads['y_primes_relative_to_ref'].value_counts())
    print(df_good_reads['anchor_name'].value_counts())

    # ---- Stats file output ----
    output_stats_f_name = os.path.join(stats_dir, f'{base_name}_stats_y_prime.txt')
    print(f'\nWriting {output_stats_f_name}')
    with open(output_stats_f_name, 'w') as f:
        if expected is not None:
            stats_neutral  = df_good_reads[df_good_reads['delta_y_prime_sign'] == 'same']['y_primes_relative_to_ref'].describe()
            stats_positive = df_good_reads[df_good_reads['delta_y_prime_sign'] == '+']['y_primes_relative_to_ref'].describe()
            stats_negative = df_good_reads[df_good_reads['delta_y_prime_sign'] == '-']['y_primes_relative_to_ref'].describe()
            f.write(f'Neutral Y prime ends:\n{stats_neutral}\n'
                    f'Positive Y prime ends:\n{stats_positive}\n'
                    f'Negative Y prime ends:\n{stats_negative}\n')
        else:
            f.write(f'No Y prime reference (expected copies per end unknown).\n'
                    f'Y primes per read:\n{df_good_reads["y_prime_probe_count"].describe()}\n')

    # ---- Call plots on good telomere reads ----
    # The delta plots come first: their strip-plot jitter draws from the global random state, so a
    # plot drawn before them would move their points.
    if expected is not None:
        df_y_prime = df_good_reads.dropna(subset=['y_primes_relative_to_ref', 'chr_end'])

        y_prime_count_violin_strip_plot(df_y_prime, base_name, figures_dir)

        delta_sign_violin_strip_plot(df_y_prime, base_name, figures_dir)

    absolute_count_plots(df_good_reads, base_name, figures_dir, with_median_view=expected is None)


if __name__ == '__main__':
    main()
