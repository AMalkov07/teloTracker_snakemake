import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.ticker import FuncFormatter

# Get the read direction in relation to the reference
def read_direction(match_start_on_anchor, match_end_on_anchor):
    if match_start_on_anchor < match_end_on_anchor:
        return 'forward'
    else:
        return 'reverse'

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

# Label the expected Y primes for each end based on Day 0 oberservations, for a given genotype
def y_prime_change_calc(anchor_name, base_name):
    if '6991' in base_name: # WildType

        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0, 'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0, 'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 1, 'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 7,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1, 'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }

    elif 'KRLT3' in base_name or '7637' in base_name: # RNAse H Mutant

        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 1, 'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 4,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 2, 'chr6R_anchor': 0, 'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 0, 'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 7,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1, 'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }

    else: # No strain, use WildType '6991'

        y_prime_numbers = {
            'chr1L_anchor': 0, 'chr1R_anchor': 0, 'chr2L_anchor': 1, 'chr2R_anchor': 0, 'chr3L_anchor': 0, 'chr3R_anchor': 0, 'chr4L_anchor': 0, 'chr4R_anchor': 7,
            'chr5L_anchor': 1, 'chr5R_anchor': 1, 'chr6L_anchor': 1, 'chr6R_anchor': 0, 'chr7L_anchor': 0, 'chr7R_anchor': 1, 'chr8L_anchor': 1, 'chr8R_anchor': 1,
            'chr9L_anchor': 1, 'chr9R_anchor': 0, 'chr10L_anchor': 1, 'chr10R_anchor': 1, 'chr11L_anchor': 0, 'chr11R_anchor': 0, 'chr12L_anchor': 1, 'chr12R_anchor': 7,
            'chr13L_anchor': 1, 'chr13R_anchor': 0, 'chr14L_anchor': 5, 'chr14R_anchor': 1, 'chr15L_anchor': 0, 'chr15R_anchor': 1, 'chr16L_anchor': 1, 'chr16R_anchor': 1
        }

    ref_y_prime_num = y_prime_numbers[anchor_name]

    return ref_y_prime_num

### Creating new categorical columns for the dataframe
def assign_delta_y_prime_sign(y_primes_relative_to_ref):
    if y_primes_relative_to_ref > 0:
        delta_y_prime_sign = "+"
    elif y_primes_relative_to_ref < 0:
        delta_y_prime_sign = "-"
    else: # y_primes_relative_to_ref == 0:
        delta_y_prime_sign = "same"
    return delta_y_prime_sign

def assign_reference_y_prime_end_status(reference_y_primes):
    if reference_y_primes == 0:
        reference_end_status = "Reference Y Primes = 0"
    if reference_y_primes == 1:
        reference_end_status = "Reference Y Primes = 1"
    elif  reference_y_primes >= 2:
        reference_end_status = "Reference Y Primes = 2+"
    return reference_end_status

def assign_read_y_prime_end_status(reference_y_primes):
    if reference_y_primes == 0:
        reference_end_status = "Read Y Primes = 0"
    if reference_y_primes == 1:
        reference_end_status = "Read Y Primes = 1"
    elif  reference_y_primes >= 2:
        reference_end_status = "Read Y Primes = 2+"
    return reference_end_status

###

# Plots of Y Prime info
def y_prime_count_violin_strip_plot(dataframe, section='all', x_plot='chr_end', y_plot='y_primes_relative_to_ref', plot_scale=(22,9)):

    title_name = f'{base_name}_{section}' #title_name = f'{base_name}_{end_protection}_{repeat_measure.split("_")[1]}'
    sns.set(style="whitegrid")

    # Set the color palette using the custom colors
    chr_color_code = ['#641E16', '#512E5F', '#4A235A', '#154360', '#2E86C1', '#0E6251', '#0B5345', '#145A32', '#7D6608', '#7E5109', '#D35400', '#626567', '#4D5656', '#212F3C', '#17202A', '#F93409', '#FFED24',
                  '#CD6155', '#D7BDE2', '#BB8FCE', '#A9CCE3', '#AED6F1', '#A3E4D7', '#45B39D', '#229954', '#F7DC6F', '#F5B041', '#E59866', '#D7DBDD', '#95A5A6', '#5D6D7E', '#ABB2B9', '#FFC2B4', '#FFFEB4']
    sns.set_palette(chr_color_code)

    #Set the plotting size (make larger to fit swarm plot points)
    fig, ax = plt.subplots(figsize=(plot_scale))

    g1 = sns.violinplot(x=x_plot, order=chr_list, y=y_plot, data=dataframe, gridsize=1000, cut=0, palette=chr_color_code) #inner='quartile' , hue=x_plot, hue_order=chr_list

    g2 = sns.stripplot(x=x_plot, order=chr_list, y=y_plot, data=dataframe, linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="#FF2400", ax=g1) # ax=g1, hue=x_plot, dodge=False, density_norm='count'


    plt.legend().remove()
    total_reads_in_plot=len(dataframe)
    fig_title_name = title_name.strip('guppy-')
    fig_title_name = fig_title_name.split('_PromethION')[0]
    fig_title_name = fig_title_name.split('_MinION')[0]
    ax.set_title(f'{fig_title_name} Delta Y Primes (N = {total_reads_in_plot} Reads, Ave. Y Prime to ref = {(dataframe[y_plot].mean()):.3f})', fontweight="bold", fontsize = 20, color='k', pad=15)
    plt.xlabel("Chromosome End", fontweight="bold", fontsize = 20)   #fontsize = 0
    plt.ylabel("Delta Y Primes", fontweight="bold", fontsize = 30)
    plt.xticks(fontweight="bold", fontsize=13)   #fontsize = 0  fontsize=25
    plt.yticks(fontweight="bold", fontsize=20)

    # Write the read % and expected number of Y primes for each chr_end
    for i, anchor in enumerate(dataframe[x_plot].unique()):
        num_reads = dataframe[dataframe[x_plot] == anchor][y_plot].count()
        plt.text(i, (ax.get_ylim()[1]-0.4), f'{((num_reads)*100/(total_reads_in_plot)):.0f}%', ha='center', fontsize=10, fontweight="bold")

        ref_y_prime_num_in_plot = (dataframe[dataframe[x_plot] == anchor]['reference_y_primes'].unique())[0] # ['reference_y_primes'][0]

        plt.text(i, (ax.get_ylim()[0]+0.2), f'Ref={ref_y_prime_num_in_plot}', ha='center', fontsize=10, fontweight="bold")

    if '6991' in base_name: # WildType
        figure_directory = 'figures_for_y_primes/wild_type/'
    elif '7021' in base_name: # rad52
        figure_directory = 'figures_for_y_primes/rad52/'
    elif '7093' in base_name: # rad52-annealing
        figure_directory = 'figures_for_y_primes/rad52-annealing/'
    elif '7154' in base_name: # rad51
        figure_directory = 'figures_for_y_primes/rad51/'
    elif '7172' in base_name: # mph1
        figure_directory = 'figures_for_y_primes/mph1/'
    elif '7174' in base_name: # rad59
        figure_directory = 'figures_for_y_primes/rad59/'
    elif '7250' in base_name: # rad54
        figure_directory = 'figures_for_y_primes/rad54/'
    elif '7321' in base_name: #pif1-m2
        figure_directory = 'figures_for_y_primes/pif1_m2/'
    elif '7323' in base_name: #pif1-NLS
        figure_directory = 'figures_for_y_primes/pif1_nls/'
    elif '7324' in base_name: # rif1rif2
        figure_directory = 'figures_for_y_primes/rif1rif2/'
    elif '7302' in base_name: # mph1
        figure_directory = 'figures_for_y_primes/mph1/'
    else:
        figure_directory = 'figures_for_y_primes/'

    save_file_name = f'{figure_directory}/{title_name}_delta_y_primes.png' #{figure_directory}/
    print(f'Graphing {save_file_name}...')
    os.makedirs(figure_directory, exist_ok=True)
    plt.savefig(f'{save_file_name}', dpi=300, format="png")

    plt.close()

def delta_sign_violin_strip_plot(dataframe, x_plot='delta_y_prime_sign', y_plot='y_primes_relative_to_ref', plot_y_max=None, plot_scale=(16,9)):

    title_name = f'{base_name}' #title_name = f'{base_name}_{end_protection}_{repeat_measure.split("_")[1]}'
    sns.set(style="whitegrid")

    #Set the plotting size (make larger to fit swarm plot points)
    fig, ax = plt.subplots(figsize=(plot_scale))

    g1 = sns.violinplot(x=x_plot, hue=x_plot, hue_order=["-", "+"], order=["-", "+"], y=y_plot, data=dataframe, gridsize=1000, cut=0, palette=['tab:red', 'tab:green']) #inner='quartile' , hue=x_plot, hue_order=chr_list

    #Final figure y-axis set below
    g1.set(ylim=(-20, 20))


    g2 = sns.stripplot(x=x_plot, order=["-", "+"], y=y_plot, data=dataframe, linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="tab:gray", ax=g1) # ax=g1, hue=x_plot, dodge=False, density_norm='count'

    g2.set(ylim=(-10, 20))

    plt.legend().remove()
    total_reads_in_plot=len(dataframe)
    total_positive_reads_in_plot = len(dataframe.loc[(dataframe['delta_y_prime_sign'] == "+")])
    total_negative_reads_in_plot = len(dataframe.loc[(dataframe['delta_y_prime_sign'] == "-")])

    fig_title_name = title_name.strip('guppy-')
    fig_title_name = fig_title_name.split('_PromethION')[0]
    fig_title_name = fig_title_name.split('_MinION')[0]
    ax.set_title(f'{fig_title_name} Delta Y Primes (N = {total_reads_in_plot} Reads, - = {total_negative_reads_in_plot} & + = {total_positive_reads_in_plot})',
                    fontweight="bold", fontsize = 20, color='k', pad=15)
    plt.xlabel("Sign of Delta Y Primes", fontweight="bold", fontsize = 20)   #fontsize = 0
    plt.ylabel("Number of Delta Y Primes", fontweight="bold", fontsize = 30)
    plt.xticks(fontweight="bold", fontsize=13)   #fontsize = 0  fontsize=25
    plt.yticks(fontweight="bold", fontsize=20)

    if '6991' in base_name: # WildType
        figure_directory = 'figures_for_y_primes/wild_type/'
    elif '7021' in base_name: # rad52
        figure_directory = 'figures_for_y_primes/rad52/'
    elif '7093' in base_name: # rad52-annealing
        figure_directory = 'figures_for_y_primes/rad52-annealing/'
    elif '7154' in base_name: # rad51
        figure_directory = 'figures_for_y_primes/rad51/'
    elif '7172' in base_name: # mph1
        figure_directory = 'figures_for_y_primes/mph1/'
    elif '7174' in base_name: # rad59
        figure_directory = 'figures_for_y_primes/rad59/'
    elif '7250' in base_name: # rad54
        figure_directory = 'figures_for_y_primes/rad54/'
    elif '7321' in base_name: #pif1-m2
        figure_directory = 'figures_for_y_primes/pif1_m2/'
    elif '7323' in base_name: #pif1-NLS
        figure_directory = 'figures_for_y_primes/pif1_nls/'
    elif '7324' in base_name: # rif1rif2
        figure_directory = 'figures_for_y_primes/rif1rif2/'
    elif '7302' in base_name: # mph1
        figure_directory = 'figures_for_y_primes/mph1/'
    else:
        figure_directory = 'figures_for_y_primes/'

    save_file_name = f'{figure_directory}/{title_name}_sign_delta_y_primes.png' #{figure_directory}/
    print(f'Graphing {save_file_name}...')
    os.makedirs(figure_directory, exist_ok=True)
    plt.savefig(f'{save_file_name}', dpi=300, format="png")

    plt.close()

    #plt.show()

def all_single_violin_strip_plot(dataframe, x_plot='delta_y_prime_sign', y_plot='y_primes_relative_to_ref', plot_y_max=None, plot_scale=(16,9)):

    title_name = f'{base_name}' #title_name = f'{base_name}_{end_protection}_{repeat_measure.split("_")[1]}'
    sns.set(style="whitegrid")

    #Set the plotting size (make larger to fit swarm plot points)
    fig, ax = plt.subplots(figsize=(plot_scale))

    g1 = sns.violinplot(x=x_plot, hue=x_plot, hue_order=["-", "same", "+"], order=["-", "same", "+"], y=y_plot, data=dataframe, gridsize=1000, cut=0, palette=['tab:red', 'tab:blue', 'tab:green'])

    #Final figure y-axis set below
    g1.set(ylim=(-10, 20))


    g2 = sns.stripplot(x=x_plot, order=["-", "same", "+"], y=y_plot, data=dataframe, linewidth=0.5, alpha=0.6, edgecolor="k", s=4, color="tab:gray", ax=g1)

    g2.set(ylim=(-10, 20))

    plt.legend().remove()
    total_reads_in_plot=len(dataframe)
    total_positive_reads_in_plot = len(dataframe.loc[(dataframe['delta_y_prime_sign'] == "+")])
    total_negative_reads_in_plot = len(dataframe.loc[(dataframe['delta_y_prime_sign'] == "-")])
    total_same_reads_in_plot = len(dataframe.loc[(dataframe['delta_y_prime_sign'] == "same")])

    percent_positive_reads_in_plot = f'{(total_positive_reads_in_plot / total_reads_in_plot)*100:.2f}'
    percent_negative_reads_in_plot = f'{(total_negative_reads_in_plot / total_reads_in_plot)*100:.2f}'
    percent_same_reads_in_plot = f'{(total_same_reads_in_plot / total_reads_in_plot)*100:.2f}'

    fig_title_name = title_name.strip('guppy-')
    fig_title_name = fig_title_name.split('_PromethION')[0]
    fig_title_name = fig_title_name.split('_MinION')[0]
    ax.set_title(f'{fig_title_name} Delta Y Primes (N = {total_reads_in_plot} Reads, - = {percent_negative_reads_in_plot}%, same = {percent_same_reads_in_plot}%, + = {percent_positive_reads_in_plot}%)',
                    fontweight="bold", fontsize = 20, color='k', pad=15)
    plt.xlabel("Change in Y Primes", fontweight="bold", fontsize = 20)   #fontsize = 0
    plt.ylabel("Number of Delta Y Primes", fontweight="bold", fontsize = 30)
    plt.xticks(fontweight="bold", fontsize=13)   #fontsize = 0  fontsize=25
    plt.yticks(fontweight="bold", fontsize=20)

    if '6991' in base_name: # WildType
        figure_directory = 'figures_for_y_primes/wild_type/'
    elif '7021' in base_name: # rad52
        figure_directory = 'figures_for_y_primes/rad52/'
    elif '7093' in base_name: # rad52-annealing
        figure_directory = 'figures_for_y_primes/rad52-annealing/'
    elif '7154' in base_name: # rad51
        figure_directory = 'figures_for_y_primes/rad51/'
    elif '7172' in base_name: # mph1
        figure_directory = 'figures_for_y_primes/mph1/'
    elif '7174' in base_name: # rad59
        figure_directory = 'figures_for_y_primes/rad59/'
    elif '7250' in base_name: # rad54
        figure_directory = 'figures_for_y_primes/rad54/'
    elif '7321' in base_name: #pif1-m2
        figure_directory = 'figures_for_y_primes/pif1_m2/'
    elif '7323' in base_name: #pif1-NLS
        figure_directory = 'figures_for_y_primes/pif1_nls/'
    elif '7324' in base_name: # rif1rif2
        figure_directory = 'figures_for_y_primes/rif1rif2/'
    elif '7302' in base_name: # mph1
        figure_directory = 'figures_for_y_primes/mph1/'
    else:
        figure_directory = 'figures_for_y_primes/'

    save_file_name = f'{figure_directory}/{title_name}_all_y_prime_counts.png' #{figure_directory}/
    print(f'Graphing {save_file_name}...')
    os.makedirs(figure_directory, exist_ok=True)
    plt.savefig(f'{save_file_name}', dpi=300, format="png")

    plt.close()

    #plt.show()

def plot_ref_y_prime_type_delta_barplot(dataframe, size=(16,11)):

    title_name = f'{base_name}' #f'{base_name}_{end_protection}_{repeat_measure.split("_")[1]}'

    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Calculate normalized counts for each day and outcome_status
    normalized_counts = dataframe.groupby(['reference_y_prime_end_status', 'delta_y_prime_sign']).size() / dataframe.groupby('reference_y_prime_end_status').size() * 100
    normalized_counts = normalized_counts.reset_index(name='percentage')

    print(normalized_counts)

    # Plot histogram
    g1 = sns.barplot(data=normalized_counts, x='reference_y_prime_end_status', y='percentage', hue='delta_y_prime_sign',
                        palette=['#e80008', '#00b6ff', '#00ff16'],
                        hue_order=['-', 'same', '+'])

    plt.ylabel('Percentage of each Y Prime Type')
    plt.title(f'{title_name} Delta by Y Prime Type')
    plt.legend(title='Delta Y Prime')

    #ax.set_title(f'{title_name} Delta to Ref End (N = {num_reads_delta})', fontweight="bold", fontsize = 20, color='k', pad=15)
    #ax.set_xlabel("Delta to Reference End", fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    #ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    #ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)

    # Set the bounds of the x and y axes
    #ax.set_xlim(left=0, right=500)
    #ax.set_ylim(bottom=0, top=12)

    # Set background color
    ax.set_facecolor('w')

    if '6991' in base_name: # WildType
        figure_directory = 'figures_for_y_primes/wild_type/'
    elif '7021' in base_name: # rad52
        figure_directory = 'figures_for_y_primes/rad52/'
    elif '7093' in base_name: # rad52-annealing
        figure_directory = 'figures_for_y_primes/rad52-annealing/'
    elif '7154' in base_name: # rad51
        figure_directory = 'figures_for_y_primes/rad51/'
    elif '7172' in base_name: # mph1
        figure_directory = 'figures_for_y_primes/mph1/'
    elif '7174' in base_name: # rad59
        figure_directory = 'figures_for_y_primes/rad59/'
    elif '7250' in base_name: # rad54
        figure_directory = 'figures_for_y_primes/rad54/'
    elif '7321' in base_name: #pif1-m2
        figure_directory = 'figures_for_y_primes/pif1_m2/'
    elif '7323' in base_name: #pif1-NLS
        figure_directory = 'figures_for_y_primes/pif1_nls/'
    elif '7324' in base_name: # rif1rif2
        figure_directory = 'figures_for_y_primes/rif1rif2/'
    elif '7302' in base_name: # mph1
        figure_directory = 'figures_for_y_primes/mph1/'
    else:
        figure_directory = 'figures_for_y_primes/'

    os.makedirs(figure_directory, exist_ok=True)
    save_file_name = f'{figure_directory}/{title_name}_ref_y_prime_type_delta.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    #plt.show()
    plt.close()

if __name__ == "__main__":


    print("Starting y_prime_analysis.py")

    base_name =  f'{sys.argv[1]}'  #file_name         # f'{sys.argv[1]}'
    anchor_set = f'{sys.argv[2]}'  #'1003_anchors'        #f'{sys.argv[2]}'

    y_prime_blast_dir=f'results/{base_name}/y_prime_blast_results_from_chr_anchored_reads'

    y_prime_blast_dir_y_primes_mega=f'{y_prime_blast_dir}/y_primes_mega'
    y_prime_blast_dir_probe=f'{y_prime_blast_dir}/probes'

    top_anchor_blast_input=f'results/{base_name}/top_matches_{base_name}_blasted_{anchor_set}.tsv'
    telomere_repeat_results_input=f'results/{base_name}/{base_name}_post_telo_trimming.tsv'

    output_tsv_f_name=f'{y_prime_blast_dir_probe}/all_{base_name}_probe_matches.tsv'
    new_output_file=f'results/outputs/{base_name}_post_y_prime_probe.tsv'

    print(f"Starting y_prime_analysis.py for {base_name}")

    probe_blasting_files = os.listdir(f"{y_prime_blast_dir_probe}/")
    files_to_analyze_list = [f'{y_prime_blast_dir_probe}/{f}' for f in probe_blasting_files if f.endswith('blasted_probe.tsv')]   # y_primes_for_blasting.tsv


    df = pd.DataFrame()
    for file in files_to_analyze_list:
        df_temp = pd.read_csv(file, sep='\t')
        df = pd.concat([df, df_temp])

    df.to_csv(output_tsv_f_name, sep = '\t')

    column_type_dict = {'read_id': str, 'total_read_length': int, 'read_bp_used_for_match': int, 'match_start_on_read': int, 'match_end_on_read': int,
                        'anchor_name': str, 'total_anchor_length': int, 'match_start_on_anchor': int, 'match_end_on_anchor': int, 'pident': float, 'bitscore': float, 'evalue' : float}

    df = df.astype(column_type_dict)

    # Filter for only reads with more than 2500 bp matching the anchor
    df_filtered = df[df['read_bp_used_for_match'] > df['total_anchor_length']/2]

    # Sort the values by best matches (match_length is consistent with bitscore)
    df_filtered = df_filtered.sort_values(['read_id', 'match_start_on_read', 'bitscore'], ascending=False)

    # Get each read's length of the anchor it matches
    df_filtered["match_length"] = df_filtered.apply(lambda row: match_length_calc(row["match_start_on_anchor"], row["match_end_on_anchor"]), axis=1)

    # Get the read direction in relation to the reference
    chr_l_end_list = ['chr10L_anchor', 'chr11L_anchor', 'chr12L_anchor', 'chr13L_anchor', 'chr14L_anchor', 'chr15L_anchor', 'chr16L_anchor', 'chr1L_anchor',
                    'chr2L_anchor', 'chr3L_anchor', 'chr4L_anchor', 'chr5L_anchor', 'chr6L_anchor', 'chr7L_anchor', 'chr8L_anchor', 'chr9L_anchor' ]
    df_filtered['l_end_chr'] = df_filtered['anchor_name'].apply(lambda x : x in chr_l_end_list)

    # Get each read's alignment direction (direction in relation to ref)
    df_filtered["alignment_direction"] = df_filtered.apply(lambda row: read_direction(row["match_start_on_anchor"], row["match_end_on_anchor"]), axis=1)

    # Get each read's Repeat Type AC/TG
    df_filtered["repeat_type"] = df_filtered.apply(lambda row: repeat_type_of_read(row["l_end_chr"], row["alignment_direction"]), axis=1)

    read_y_prime_counts = df_filtered['read_id'].value_counts()

    # Converting the Series to a DataFrame
    df_read_y_prime_counts = pd.DataFrame(read_y_prime_counts.reset_index())

    df_read_y_prime_counts.columns = ['read_id', 'y_prime_probe_count']


    df_repeat_results = pd.read_csv(telomere_repeat_results_input, sep='\t', index_col=0)


    df_repeat_results['read_id'] = df_repeat_results['read_id'].astype(str)


    df_combined_results = pd.merge(df_repeat_results, df_read_y_prime_counts, on='read_id', how='outer')

    df_combined_results['y_prime_probe_count'].fillna(0, inplace=True)

    df_combined_results['reference_y_primes'] = df_combined_results["anchor_name"].apply(lambda anchor_name: y_prime_change_calc(anchor_name, base_name))
    df_combined_results['y_primes_relative_to_ref'] = df_combined_results['y_prime_probe_count'] - df_combined_results['reference_y_primes']
    df_combined_results['delta_y_prime_sign'] = df_combined_results.apply(lambda row: assign_delta_y_prime_sign(row["y_primes_relative_to_ref"]), axis=1)
    df_combined_results['reference_y_prime_end_status'] = df_combined_results.apply(lambda row: assign_reference_y_prime_end_status(row["reference_y_primes"]), axis=1)
    df_combined_results['read_y_prime_end_status'] = df_combined_results.apply(lambda row: assign_read_y_prime_end_status(row["y_prime_probe_count"]), axis=1)
    df_combined_results['y_prime_type_delta_y_prime_sign'] = df_combined_results['delta_y_prime_sign'] + " " + df_combined_results['reference_y_prime_end_status']

    print(df_combined_results)

    chr_list = ['1L', '1R', '2L', '2R', '3L', '3R', '4L', '4R', '5L', '5R', '6L', '6R', '7L', '7R', '8L', '8R', '9L', '9R',
                '10L', '10R', '11L', '11R', '12L', '12R', '13L', '13R', '14L', '14R', '15L', '15R', '16L', '16R']

    chr_sorter = dict(zip(chr_list, range(len(chr_list))))

    # Order the chromosomes from 1L, 1R, ..., 16L, 16R
    df_combined_results['chr_end'] = df_combined_results['anchor_name'].apply(lambda anchor_name: anchor_name.strip('chr _anchor'))
    df_combined_results['chr_end_rank'] = df_combined_results['chr_end'].map(chr_sorter)
    df_combined_results.sort_values(by='chr_end_rank', ascending=True, inplace=True)

    # Make files
    df_combined_results.to_csv(new_output_file, sep="\t", index=False)

    print(f'Lenth of old repeat results: {len(df_repeat_results)}')
    print(f'Lenth of y prime results: {len(df_read_y_prime_counts)}')
    print(f'Lenth of combined: {len(df_combined_results)}')

    print(df_combined_results)

    print('All df_combined_results')
    print(df_combined_results['y_prime_probe_count'].value_counts())
    print(df_combined_results['y_primes_relative_to_ref'].value_counts())
    print('\n\n')

    df_good_telomeres = df_combined_results.loc[(df_combined_results['repeat_length'] > 40) & (df_combined_results['Adapter_After_Telomere'] == True )]
    print('Reads to Telomere df_good_telomeres')
    print(df_good_telomeres)
    print(df_good_telomeres['y_prime_probe_count'].value_counts())
    print(df_good_telomeres['y_primes_relative_to_ref'].value_counts())
    print(df_good_telomeres['anchor_name'].value_counts())
    print('\n\n')




    df_good_telomeres.sort_values(by='chr_end_rank', ascending=True, inplace=True)

    # Make plots for the Y prime analysis
    y_prime_count_violin_strip_plot(df_good_telomeres)

    all_single_violin_strip_plot(dataframe=df_good_telomeres)

    df_good_telomeres_positive = df_good_telomeres[df_good_telomeres['y_primes_relative_to_ref'] > 0]

    y_prime_count_violin_strip_plot(dataframe=df_good_telomeres_positive, section = "+")

    plot_ref_y_prime_type_delta_barplot(dataframe=df_good_telomeres)

    df_delta_sign = df_good_telomeres[df_good_telomeres['delta_y_prime_sign'] != 'same']

    delta_sign_violin_strip_plot(dataframe=df_delta_sign)

    # Write the stats down
    stats_of_neutral_sign = df_good_telomeres[df_good_telomeres['delta_y_prime_sign'] == "same"]['y_primes_relative_to_ref'].describe()
    stats_of_postive_sign = df_good_telomeres[df_good_telomeres['delta_y_prime_sign'] == "+"]['y_primes_relative_to_ref'].describe()
    stats_of_negative_sign = df_good_telomeres[df_good_telomeres['delta_y_prime_sign'] == "-"]['y_primes_relative_to_ref'].describe()

    if "guppy" in base_name:
        output_stats_f_name = f'stats_for_y_primes/{base_name.split("-")[1]}_new_stats_y_prime.txt'
    elif "dorado" in base_name:
        output_stats_f_name = f'stats_for_y_primes/{base_name.split("dorado_")[1]}_new_stats_y_prime.txt'
    else:
        print("Error")

    print(f'Writting {output_stats_f_name}')
    with open (output_stats_f_name, 'w') as f:
        f.write(f'Neutral Y prime ends:\n{stats_of_neutral_sign}\nPositive Y prime ends:\n{stats_of_postive_sign}\nNegative Y prime ends:\n{stats_of_negative_sign}\n')

    plt.close('all')
    sys.exit(0)

