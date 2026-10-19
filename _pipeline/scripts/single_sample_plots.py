import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

from telomere_filters import qualifying_mask

# Usage: single_sample_plots.py <base_name> <output_dir> [--no-adapter-filter]
base_name  = sys.argv[1]
output_dir = sys.argv[2]
require_adapter = '--no-adapter-filter' not in sys.argv[3:]   # the adapter check was off

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

    # FacetGrid shares one y axis, so set_ylim inside the loop rescaled every
    # panel rather than the current one, compounding 1.1x over all 32 facets.
    y_max = 50

    for ax, anchor_name in zip(g.axes.flat, chr_list):
        subset = dataframe[dataframe['anchor_name'] == anchor_name]
        median_value = np.round(subset['repeat_length'].median())
        ax.axvline(median_value, color='black', linestyle='--', linewidth=2)
        ax.text(ax.get_xlim()[1] * 0.9, y_max * 0.9, f'Median: {int(median_value)}',
                horizontalalignment='right', fontsize=12, fontweight='bold', color='black')

    for ax in g.axes.flat:
        ax.set_ylim(0, y_max)

    g.set_titles("{col_name}")
    g.set_axis_labels('Telomere Repeat Length (bp)', "Frequency (%)")
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(f'{title_name} (N = {num_reads_in_graph})', fontweight="bold", fontsize=20, color='k')

    save_file_name = os.path.join(telomere_figures, f'{title_name}_facet_by_chr.png')
    print(f'Graphing {save_file_name}...')
    # 32 panels of 7.5 x 5 in make a ~30 x 40 in figure: at 300 dpi that was ~9,000 x 12,000 px
    # (~105 MP), too big for browsers / VS Code to show. 100 dpi keeps the layout (~12 MP).
    g.savefig(save_file_name, dpi=100, format="png")
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


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

df_all   = pd.read_csv(input_tsv, sep='\t')
df_graph = df_all[qualifying_mask(df_all, require_adapter)]

# Telomere length plots
histogram_plot_500bp(df_graph)
histogram_plot_500bp_chr_facet(df_graph)

for chr_value in df_graph['anchor_name'].unique():
    df_chr = df_graph[df_graph['anchor_name'] == chr_value]
    histogram_plot_individual_chr_repeat(df_chr, chr_anchor=chr_value)

histogram_plot_10kb(df_graph)
