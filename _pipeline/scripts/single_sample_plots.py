import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

# Usage: single_sample_plots.py <base_name> <output_dir>
base_name  = sys.argv[1]
output_dir = sys.argv[2]

print("Starting single_sample_plots.py")
print(f'Opening {base_name}...')

figure_directory = os.path.join(output_dir, 'graphs')
telomere_figures = os.path.join(figure_directory, 'telomere_figures')
ind_chr_figures  = os.path.join(telomere_figures, 'individual_chromosomes')
y_prime_figures  = os.path.join(figure_directory, 'y_prime_figures')

os.makedirs(figure_directory, exist_ok=True)
os.makedirs(telomere_figures, exist_ok=True)
os.makedirs(ind_chr_figures, exist_ok=True)
os.makedirs(y_prime_figures, exist_ok=True)

input_tsv = os.path.join(output_dir, f'{base_name}_post_y_prime_probe.tsv')


# ---------------------------------------------------------------------------
# Single-sample Plotting functions
# ---------------------------------------------------------------------------

def histogram_plot_500bp(dataframe, size=(16, 11)):
    title_name = f'{base_name.split("_no_tag")[0]}'
    num_reads_in_graph = len(dataframe)
    sns.set(style="ticks")
    fig, ax = plt.subplots(figsize=size)
    sns.histplot(data=dataframe, x='repeat_length', binwidth=10, binrange=(0, 500), kde=True, stat="percent", color='royalblue')
    ax.set_title(f'{title_name} (N = {num_reads_in_graph})', fontweight="bold", fontsize=20, color='k', pad=15)
    ax.set_xlabel('Telomere Repeat Length (bp)', fontweight="bold", fontsize=30, color='k', labelpad=6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize=30, color='k', labelpad=6)
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)
    ax.set_yticks([2, 4, 6, 8, 10, 12])
    ax.set_xticks([0, 50, 100, 200, 300, 400, 500])
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=12)
    ax.set_facecolor('w')
    save_file_name = os.path.join(telomere_figures, f'{title_name}_500bp.png')
    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=300, format="png")
    plt.close()


def histogram_plot_500bp_chr_facet(dataframe, size=(16, 11)):
    title_name = f'{base_name.split("_no_tag")[0]}'
    num_reads_in_graph = len(dataframe)
    sns.set(style="ticks")

    chr_list = ['chr1R_anchor', 'chr2R_anchor', 'chr3R_anchor', 'chr4R_anchor', 'chr5R_anchor',
                'chr6R_anchor', 'chr7R_anchor', 'chr8R_anchor', 'chr9R_anchor', 'chr10R_anchor',
                'chr11R_anchor', 'chr12R_anchor', 'chr13R_anchor', 'chr14R_anchor', 'chr15R_anchor',
                'chr16R_anchor', 'chr17R_anchor', 'chr1L_anchor', 'chr2L_anchor', 'chr3L_anchor',
                'chr4L_anchor', 'chr5L_anchor', 'chr6L_anchor', 'chr7L_anchor', 'chr8L_anchor',
                'chr9L_anchor', 'chr10L_anchor', 'chr11L_anchor', 'chr12L_anchor', 'chr13L_anchor',
                'chr14L_anchor', 'chr15L_anchor', 'chr16L_anchor', 'chr17L_anchor']

    # Only include chromosomes that actually appear in the data
    chr_list = [c for c in chr_list if c in dataframe['anchor_name'].unique()]

    chr_color_code = sns.color_palette("husl", len(chr_list))

    g = sns.FacetGrid(dataframe, col="anchor_name", col_wrap=4, height=5, aspect=1.5,
                      hue="anchor_name", palette=chr_color_code, col_order=chr_list)

    g.map(sns.histplot, 'repeat_length', binwidth=25, binrange=(0, 500), stat="percent", alpha=0.5)

    for ax, (anchor_name, subset) in zip(g.axes.flat, dataframe.groupby("anchor_name", sort=False)):
        median_value = np.round(subset['repeat_length'].median())
        ax.axvline(median_value, color='black', linestyle='--', linewidth=2)
        ax.text(ax.get_xlim()[1] * 0.9, ax.get_ylim()[1] * 0.9, f'Median: {int(median_value)}',
                horizontalalignment='right', fontsize=12, fontweight='bold', color='black')
        max_y = ax.get_ylim()[1]
        ax.set_ylim(0, max_y * 1.1)

    g.set_titles("{col_name}")
    g.set_axis_labels('Telomere Repeat Length (bp)', "Frequency (%)")
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(f'{title_name} (N = {num_reads_in_graph})', fontweight="bold", fontsize=20, color='k')

    save_file_name = os.path.join(telomere_figures, f'{title_name}_facet_by_chr.png')
    print(f'Graphing {save_file_name}...')
    g.savefig(save_file_name, dpi=300, format="png")
    plt.close()


def histogram_plot_individual_chr_repeat(dataframe, size=(16, 11), chr_anchor=None):
    chr_color_dict = {
        'chr1R_anchor': '#641E16', 'chr2R_anchor': '#512E5F', 'chr3R_anchor': '#4A235A',
        'chr4R_anchor': '#154360', 'chr5R_anchor': '#2E86C1', 'chr6R_anchor': '#0E6251',
        'chr7R_anchor': '#0B5345', 'chr8R_anchor': '#145A32', 'chr9R_anchor': '#7D6608',
        'chr10R_anchor': '#7E5109', 'chr11R_anchor': '#D35400', 'chr12R_anchor': '#626567',
        'chr13R_anchor': '#4D5656', 'chr14R_anchor': '#212F3C', 'chr15R_anchor': '#17202A',
        'chr16R_anchor': '#F93409',
        'chr1L_anchor': '#CD6155', 'chr2L_anchor': '#D7BDE2', 'chr3L_anchor': '#BB8FCE',
        'chr4L_anchor': '#A9CCE3', 'chr5L_anchor': '#AED6F1', 'chr6L_anchor': '#A3E4D7',
        'chr7L_anchor': '#45B39D', 'chr8L_anchor': '#229954', 'chr9L_anchor': '#F7DC6F',
        'chr10L_anchor': '#F5B041', 'chr11L_anchor': '#E59866', 'chr12L_anchor': '#D7DBDD',
        'chr13L_anchor': '#95A5A6', 'chr14L_anchor': '#5D6D7E', 'chr15L_anchor': '#ABB2B9',
        'chr16L_anchor': '#FFC2B4'
    }
    title_name = base_name
    num_reads_in_graph = len(dataframe)
    sns.set(style="ticks")
    fig, ax = plt.subplots(figsize=size)
    sns.histplot(data=dataframe, x='repeat_length', binwidth=10, binrange=(0, 500), kde=True,
                 stat="percent", color=chr_color_dict.get(chr_anchor, 'royalblue'))
    ax.set_title(f'{title_name} {chr_anchor} (N = {num_reads_in_graph})', fontweight="bold", fontsize=20, color='k', pad=15)
    ax.set_xlabel('Telomere Repeat Length (bp)', fontweight="bold", fontsize=30, color='k', labelpad=6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize=30, color='k', labelpad=6)
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)
    ax.set_yticks([2, 4, 6, 8, 10, 12])
    ax.set_xticks([0, 50, 100, 200, 300, 400, 500])
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=12)
    ax.set_facecolor('w')
    chr_fig_dir = os.path.join(ind_chr_figures, chr_anchor.split('_')[0])
    os.makedirs(chr_fig_dir, exist_ok=True)
    save_file_name = os.path.join(chr_fig_dir, f'{title_name}_500bp_individual_chr.png')
    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=100, format="png")
    plt.close()


def histogram_plot_10kb(dataframe, size=(16, 11)):
    title_name = f'{base_name.split("_no_tag")[0]}'
    num_reads_in_graph = len(dataframe)
    sns.set(style="ticks")
    fig, ax = plt.subplots(figsize=size)
    sns.histplot(data=dataframe, x='repeat_length', binwidth=50, binrange=(0, 10000), stat="percent", color='darkorchid')
    ax.set_title(f'{title_name} (N = {num_reads_in_graph})', fontweight="bold", fontsize=20, color='k', pad=15)
    ax.set_xlabel('Telomere Repeat Length (bp)', fontweight="bold", fontsize=30, color='k', labelpad=6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize=30, color='k', labelpad=6)
    ax.tick_params(axis='both', which='major', direction='out', length=8, width=2, labelsize=25, pad=2)
    ax.set_yticks([1, 2, 3, 4, 5, 6, 7, 8])
    ax.set_xticks([0, 2000, 4000, 6000, 8000, 10000])
    ax.set_xlim(left=0, right=10000)
    ax.set_ylim(bottom=0, top=9)
    ax.set_facecolor('w')
    save_file_name = os.path.join(telomere_figures, f'{title_name}_10kb.png')
    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=300, format="png")
    plt.close()


def y_prime_count_violin_strip_plot(dataframe, x_plot='chr_end', y_plot='y_primes_relative_to_ref', plot_scale=(22, 9)):

    sns.set(style="whitegrid")

    chr_list = ([f'chr{n}R' for n in range(1, 17)] +
                [f'chr{n}L' for n in range(1, 17)])
    chr_list = [c for c in chr_list if c in dataframe[x_plot].unique()]

    chr_color_code = ['#641E16', '#512E5F', '#4A235A', '#154360', '#2E86C1', '#0E6251', '#0B5345', '#145A32', '#7D6608', '#7E5109',
                      '#D35400', '#626567', '#4D5656', '#212F3C', '#17202A', '#F93409', '#FFED24', '#CD6155', '#D7BDE2', '#BB8FCE',
                      '#A9CCE3', '#AED6F1', '#A3E4D7', '#45B39D', '#229954', '#F7DC6F', '#F5B041', '#E59866', '#D7DBDD', '#95A5A6',
                      '#5D6D7E', '#ABB2B9', '#FFC2B4', '#FFFEB4']
    chr_color_code = chr_color_code[:len(chr_list)]
    sns.set_palette(chr_color_code)

    fig, ax = plt.subplots(figsize=plot_scale)

    sns.violinplot(x=x_plot, order=chr_list, y=y_plot, data=dataframe, gridsize=1000, cut=0, palette=chr_color_code)
    sns.stripplot(x=x_plot, order=chr_list, y=y_plot, data=dataframe, linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="#FF2400", ax=ax)

    plt.legend().remove()
    total_reads_in_plot = len(dataframe)
    ax.set_title(f'{base_name} Delta Y Primes (N = {total_reads_in_plot} Reads, Ave. Y Prime to ref = {dataframe[y_plot].mean():.3f})',
                 fontweight="bold", fontsize=20, color='k', pad=15)
    plt.xlabel("Chromosome End", fontweight="bold", fontsize=20)
    plt.ylabel("Delta Y Primes", fontweight="bold", fontsize=30)
    plt.xticks(fontweight="bold", fontsize=13)
    plt.yticks(fontweight="bold", fontsize=20)

    for i, anchor in enumerate(chr_list):
        df_anchor = dataframe[dataframe[x_plot] == anchor]
        num_reads = df_anchor[y_plot].count()
        plt.text(i, ax.get_ylim()[1] - 0.4, f'{num_reads * 100 / total_reads_in_plot:.0f}%',
                 ha='center', fontsize=10, fontweight="bold")
        ref_y = df_anchor['reference_y_primes'].iloc[0] if len(df_anchor) > 0 else 'N/A'
        plt.text(i, ax.get_ylim()[0] + 0.2, f'Ref={ref_y}',
                 ha='center', fontsize=10, fontweight="bold")

    save_file_name = os.path.join(y_prime_figures, f'{base_name}_delta_y_primes.png')
    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=300, format="png")
    plt.close()


def delta_sign_violin_strip_plot(dataframe, x_plot='delta_y_prime_sign', y_plot='y_primes_relative_to_ref', plot_scale=(16, 9)):

    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=plot_scale)

    sns.violinplot(x=x_plot, hue=x_plot, hue_order=["-", "+"], order=["-", "+"], y=y_plot,
                   data=dataframe, gridsize=1000, cut=0, palette=['tab:red', 'tab:green'])
    ax.set(ylim=(-10, 20))
    sns.stripplot(x=x_plot, order=["-", "+"], y=y_plot, data=dataframe,
                  linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="tab:gray", ax=ax)

    plt.legend().remove()
    total_reads_in_plot = len(dataframe)
    total_positive = len(dataframe[dataframe['delta_y_prime_sign'] == "+"])
    total_negative = len(dataframe[dataframe['delta_y_prime_sign'] == "-"])

    ax.set_title(f'{base_name} Delta Y Primes (N = {total_reads_in_plot} Reads, - = {total_negative} & + = {total_positive})',
                 fontweight="bold", fontsize=20, color='k', pad=15)
    plt.xlabel("Sign of Delta Y Primes", fontweight="bold", fontsize=20)
    plt.ylabel("Number of Delta Y Primes", fontweight="bold", fontsize=30)
    plt.xticks(fontweight="bold", fontsize=13)
    plt.yticks(fontweight="bold", fontsize=20)

    save_file_name = os.path.join(y_prime_figures, f'{base_name}_sign_delta_y_primes.png')
    print(f'Graphing {save_file_name}...')
    plt.savefig(save_file_name, dpi=300, format="png")
    plt.close()


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

df_all   = pd.read_csv(input_tsv, sep='\t')
df_graph = df_all[(df_all['repeat_length'] >= 30) & (df_all['Adapter_After_Telomere'] == True)]

# Telomere length plots
histogram_plot_500bp(df_graph)
histogram_plot_500bp_chr_facet(df_graph)

for chr_value in df_graph['anchor_name'].unique():
    df_chr = df_graph[df_graph['anchor_name'] == chr_value]
    histogram_plot_individual_chr_repeat(df_chr, chr_anchor=chr_value)

histogram_plot_10kb(df_graph)

# Y prime plots
#df_y_prime = df_graph.dropna(subset=['y_primes_relative_to_ref', 'chr_end'])

#y_prime_count_violin_strip_plot(df_y_prime, section='all')

#df_delta = df_y_prime[df_y_prime['delta_y_prime_sign'].isin(['-', '+'])]
#if len(df_delta) > 0:
#    delta_sign_violin_strip_plot(df_delta)
