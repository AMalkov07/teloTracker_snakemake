import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from scipy.stats import mannwhitneyu
import pingouin as pg

base_name = sys.argv[1]
input_f_name = sys.argv[2]
plot_output_dir = sys.argv[3]

#input_base_name = sys.argv[1:]

# Ensure the output directory exists
os.makedirs(plot_output_dir, exist_ok=True)

print("Starting make_telomere_graphs.py")

def check_for_tag(base_name, dataframe=None):
    if "no_tag" in base_name:
        has_tag = False
    elif "TeloTag" in base_name:
        has_tag = True
    else:
        # Automated checking
        if type(dataframe) != type(None):
            total_good = len(dataframe)
            good_with_tag = len(dataframe[dataframe['Tag_After_Telomere'] == True])
            if (good_with_tag/total_good) >= 0.02:
                has_tag = True
            else:
                has_tag = False
        else:
            has_tag = None
            
    return has_tag


def histogram_plot_repeat(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=10, binrange=(0,500), kde=True, stat="percent", color='royalblue')

        
    ax.set_title(f'{title_name} (N = {num_reads_repeat})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)

    ax.set_yticks([2,4,6,8,10,12])
    ax.set_xticks([0,50,100,200,300,400,500])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=12)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_500bp.png'
    print(f'Graphing {save_file_name}...')
    #plt.savefig(f"{save_file_name}", dpi=300, format="png")
    #plt.close()
    #plt.show()
    save_path = os.path.join(plot_output_dir, f'{base_name}_histogram.png')
    plt.savefig(save_path)
    plt.close() # Always close plots to save memory in pipelines

    
def histogram_plot_delta(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=10, binrange=(0,500), kde=True, stat="percent", color='royalblue')

        
    ax.set_title(f'{title_name} (N = {num_reads_repeat})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)

    ax.set_yticks([2,4,6,8,10,12])
    ax.set_xticks([0,50,100,200,300,400,500])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=12)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_500bp.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_10kb(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    sns.set(style="ticks")

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=50, binrange=(0,10000), stat="percent", color='darkorchid', ax=ax)
        
    ax.set_title(f'{title_name} (N = {num_reads_10kb})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)

    '''
    ax.set_yticks([1,2,3,4,5,6,7,8,9,10])
    ax.set_xticks([0,100,200,300,400,500])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=10)
    '''
    
    ax.set_yticks([5,10,15,20,25,30,35,40,45])
    ax.set_xticks([0,2000,4000,6000,8000,10000])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=10000)
    ax.set_ylim(bottom=0, top=45)

    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_10kb.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_10kb_zoom(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    sns.set(style="ticks")

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=50, binrange=(0,10000), stat="percent", color='darkorchid')
    
    
    
    ax.set_title(f'{title_name} (N = {num_reads_10kb})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='major', direction='out', length=8, width=2, labelsize=25, pad=2)


    ax.set_yticks([1,2,3,4,5,6,7,8])
    ax.set_xticks([0,2000,4000,6000,8000,10000])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=10000)
    ax.set_ylim(bottom=0, top=9)



    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_10kb_zoom.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_10kb_close_zoom(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    sns.set(style="ticks")

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=50, binrange=(0,10000), stat="percent", color='darkorchid')
    

    
    ax.set_title(f'{title_name} (N = {num_reads_10kb})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='major', direction='out', length=8, width=2, labelsize=25, pad=2)


    ax.set_yticks([0.25,0.5,0.75,1,1.25,1.5,1.75,2])
    ax.set_xticks([0,2000,4000,6000,8000,10000])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=10000)
    ax.set_ylim(bottom=0, top=2)



    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_10kb_close_zoom.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_all_close_zoom(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    sns.set(style="ticks")

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=500, stat="percent", color='deeppink') #  binrange=(0,10000)
    

    
    ax.set_title(f'{title_name} (N = {num_reads_all})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 10, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='major', direction='out', length=8, width=2, labelsize=18, pad=2)


    ax.set_yticks([0.25,0.5,0.75,1,1.25,1.5,1.75,2])

    #ax.set_xticks([0,10000,20000,30000,40000,50000,6000,70000,80000,90000,100000])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=100000)
    ax.set_ylim(bottom=0, top=2)
    


    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_all_close_zoom.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_all_close_zoom_25(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    sns.set(style="ticks")

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=500, stat="percent", color='deeppink') #  binrange=(0,10000)
    

    
    ax.set_title(f'{title_name} (N = {num_reads_all})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 10, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='major', direction='out', length=8, width=2, labelsize=18, pad=2)


    ax.set_yticks([0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5])

    ax.set_xticks([0,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0) #, right=10000)
    ax.set_ylim(bottom=0, top=2.5)
    


    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_all_close_zoom_25.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()


def histogram_plot_repeat_chr(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name}_{end_protection}_{repeat_measure.split("_")[1]}_chr_anchored'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))



    chr_list = ['chr1R_anchor',
                'chr2R_anchor',
                'chr3R_anchor',
                'chr4R_anchor',
                'chr5R_anchor',
                'chr6R_anchor',
                'chr7R_anchor',
                'chr8R_anchor',
                'chr9R_anchor',
                'chr10R_anchor',
                'chr11R_anchor',
                'chr12R_anchor',
                'chr13R_anchor',
                'chr14R_anchor',
                'chr15R_anchor',
                'chr16R_anchor',
                'chr1L_anchor',
                'chr2L_anchor',
                'chr3L_anchor',
                'chr4L_anchor',
                'chr5L_anchor',
                'chr6L_anchor',
                'chr7L_anchor',
                'chr8L_anchor',
                'chr9L_anchor',
                'chr10L_anchor',
                'chr11L_anchor',
                'chr12L_anchor',
                'chr13L_anchor',
                'chr14L_anchor',
                'chr15L_anchor',
                'chr16L_anchor',
		]



    # Set the color palette using the custom colors
    chr_color_code = ['#641E16', '#512E5F', '#4A235A', '#154360', '#2E86C1', '#0E6251', '#0B5345', '#145A32', '#7D6608', '#7E5109', '#D35400', '#626567', '#4D5656', '#212F3C', '#17202A', '#F93409',
                  '#CD6155', '#D7BDE2', '#BB8FCE', '#A9CCE3', '#AED6F1', '#A3E4D7', '#45B39D', '#229954', '#F7DC6F', '#F5B041', '#E59866', '#D7DBDD', '#95A5A6', '#5D6D7E', '#ABB2B9', '#FFC2B4']
    sns.set_palette(chr_color_code)             
                    
    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', hue="anchor_name", hue_order=chr_list, multiple="stack", binwidth=10, binrange=(0,500), stat="percent", ) #  palette='bright'

        
    ax.set_title(f'{title_name} (N = {num_reads_repeat})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)


    ax.set_yticks([2,4,6,8,10,12])
    ax.set_xticks([0,50,100,200,300,400,500])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=12)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_500bp_chr.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_10kb_close_zoom_chr(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name}_{end_protection}_{repeat_measure.split("_")[1]}_chr_anchored'
    # Set the graph background
    sns.set(style="ticks")

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))


    chr_list = ['chr1R_anchor',
                'chr2R_anchor',
                'chr3R_anchor',
                'chr4R_anchor',
                'chr5R_anchor',
                'chr6R_anchor',
                'chr7R_anchor',
                'chr8R_anchor',
                'chr9R_anchor',
                'chr10R_anchor',
                'chr11R_anchor',
                'chr12R_anchor',
                'chr13R_anchor',
                'chr14R_anchor',
                'chr15R_anchor',
                'chr16R_anchor',
                'chr1L_anchor',
                'chr2L_anchor',
                'chr3L_anchor',
                'chr4L_anchor',
                'chr5L_anchor',
                'chr6L_anchor',
                'chr7L_anchor',
                'chr8L_anchor',
                'chr9L_anchor',
                'chr10L_anchor',
                'chr11L_anchor',
                'chr12L_anchor',
                'chr13L_anchor',
                'chr14L_anchor',
                'chr15L_anchor',
                'chr16L_anchor']
            
    # Set the color palette using the custom colors
    chr_color_code = ['#641E16', '#512E5F', '#4A235A', '#154360', '#2E86C1', '#0E6251', '#0B5345', '#145A32', '#7D6608', '#7E5109', '#D35400', '#626567', '#4D5656', '#212F3C', '#17202A', '#F93409',
                  '#CD6155', '#D7BDE2', '#BB8FCE', '#A9CCE3', '#AED6F1', '#A3E4D7', '#45B39D', '#229954', '#F7DC6F', '#F5B041', '#E59866', '#D7DBDD', '#95A5A6', '#5D6D7E', '#ABB2B9', '#FFC2B4']
    sns.set_palette(chr_color_code)       
    
    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', hue="anchor_name", hue_order=chr_list, multiple="stack", binwidth=50, binrange=(0,10000), stat="percent", palette='bright')
    

    
    ax.set_title(f'{title_name} (N = {num_reads_10kb})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='major', direction='out', length=8, width=2, labelsize=25, pad=2)


    ax.set_yticks([0.25,0.5,0.75,1,1.25,1.5,1.75,2])
    ax.set_xticks([0,2000,4000,6000,8000,10000])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=10000)
    ax.set_ylim(bottom=0, top=2)



    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_10kb_close_zoom_chr.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_repeat_overlay(dataframe, size=(16,11)) -> None:

    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x='trimmed_overhang_past_ref_telo_start', binwidth=10, binrange=(0,500), kde=True, stat="percent", hue= "end_type")

        
    ax.set_title(f'Overlay Adapter vs Tag)', fontweight="bold", fontsize = 10, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)

    '''
    ax.set_yticks([1,2,3,4,5,6,7,8,9,10])
    ax.set_xticks([0,100,200,300,400,500])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=10)
    '''

    #ax.set_yticks([1,2,3,4,5,6,7,8,9,10,11,12])
    ax.set_yticks([2,4,6,8,10,12])
    ax.set_xticks([0,100,200,300,400,500])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=8)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_vs_repeat-select.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()


def histogram_plot_individual_chr_repeat(dataframe, size=(16,11), chr_anchor=None) -> None:

    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    chr_color_dict = {'chr1R_anchor': '#641E16', 'chr2R_anchor': '#512E5F', 'chr3R_anchor': '#4A235A',
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
                    'chr16L_anchor': '#FFC2B4'}
    
    # 'chr17R_anchor': '#D4C202', 'chr17L_anchor': '#FFED24'

    title_name = base_name
    sns.set(style="ticks")

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=10, binrange=(0,500), kde=True,
                      stat="percent", color=chr_color_dict[chr_anchor])


    ax.set_title(f'{title_name} (N = {num_individul_chr_reads_repeat})', fontweight="bold", fontsize = 20, color='k', pad=15)
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)

    ax.set_yticks([2,4,6,8,10,12])
    ax.set_xticks([0,50,100,200,300,400,500])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=500)
    ax.set_ylim(bottom=0, top=12)

    # Set background color
    ax.set_facecolor('w')


    save_file_name = f'{genotype_chr_figure_directory}/{title_name}_500bp_individual_chr.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=100, format="png")
    plt.close()
    #plt.show()

def histogram_plot_repeat_chr_facet(dataframe, size=(16, 11)) -> None:

    title_name = base_name
    sns.set(style="ticks")


    chr_color_dict = {'chr1R_anchor': '#641E16', 'chr2R_anchor': '#512E5F', 'chr3R_anchor': '#4A235A',
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
                    'chr16L_anchor': '#FFC2B4'}





    # 'chr17R_anchor': '#D4C202', 'chr17L_anchor': '#FFED24'

    # Create the FacetGrid with the explicit mapping
    g = sns.FacetGrid(
        dataframe,
        col="anchor_name",
        col_wrap=4,
        height=5,
        aspect=1.5,
        col_order=chr_color_dict.keys(),
    )

    # Map the histograms and set colors based on the column name
    def plot_hist(x, color, **kwargs):
        # Extract just the chromosome name from the title
        title = plt.gca().get_title()
        # Find the chromosome name that appears in the title
        for chr_name in chr_color_dict.keys():
            if chr_name in title:
                chr_key = chr_name
                break

        sns.histplot(x=x, binwidth=10, binrange=(0, 500), stat="percent",
                    kde=True, line_kws={'linewidth': 3},
                    color=chr_color_dict[chr_key], **kwargs)

    # Apply the custom plotting function
    g.map(plot_hist, f'{repeat_measure}')

    # Add the median lines and annotations for each facet
    #for ax, (anchor_name, subset) in zip(g.axes.flat, dataframe.groupby("anchor_name")):
        # Calculate the median for the subset of the data
        #median_value = np.round(subset[f'{repeat_measure}'].median())

        # Add a dashed line at the rounded median
        #ax.axvline(median_value, color='black', linestyle='--', linewidth=2)

        # Annotate the median value in the middle-right corner of the facet
        #ax.text(ax.get_xlim()[1] * 0.9, ax.get_ylim()[1] * 0.9, f'Median: {int(median_value)}',
        #        horizontalalignment='right', fontsize=12, fontweight='bold', color='black')

    # Add title and labels to the overall plot
    g.set_titles("{col_name}")
    g.set_axis_labels(f'Telomere Repeat Length (bp)', "Frequency (%)")

    # Set the overall plot title
    plt.subplots_adjust(top=0.95)
    g.fig.suptitle(f'{title_name} (N = {num_reads_repeat})', fontweight="bold", fontsize=20, color='k')

    # Save the facet plot as an image file
    save_file_name = f'{figure_directory}/{title_name}_facet_by_chr.png'
    print(f'Graphing {save_file_name}...')
    g.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()


def histogram_plot_repeat_chr_facet_new(dataframe, repeat_measure, base_name, num_reads_repeat, figure_directory, size=(16, 11)) -> None:
    title_name = base_name
    sns.set(style="ticks")

    # Define the chromosome list (anchor names)
    chr_list = ['chr1R_anchor', 'chr2R_anchor', 'chr3R_anchor', 'chr4R_anchor', 'chr5R_anchor',
                'chr6R_anchor', 'chr7R_anchor', 'chr8R_anchor', 'chr9R_anchor', 'chr10R_anchor',
                'chr11R_anchor', 'chr12R_anchor', 'chr13R_anchor', 'chr14R_anchor', 'chr15R_anchor',
                'chr16R_anchor', 'chr17R_anchor', 'chr1L_anchor', 'chr2L_anchor', 'chr3L_anchor',
                'chr4L_anchor', 'chr5L_anchor', 'chr6L_anchor', 'chr7L_anchor', 'chr8L_anchor',
                'chr9L_anchor', 'chr10L_anchor', 'chr11L_anchor', 'chr12L_anchor', 'chr13L_anchor',
                'chr14L_anchor', 'chr15L_anchor', 'chr16L_anchor', 'chr17L_anchor']

    # Set the color palette
    chr_color_code = sns.color_palette("husl", len(chr_list))

    # Create the FacetGrid for facet plotting
    g = sns.FacetGrid(dataframe, col="anchor_name", col_wrap=4, height=5, aspect=1.5,
                      hue="anchor_name", palette=chr_color_code, col_order=chr_list)

    # Map the histograms onto the FacetGrid
    g.map(sns.histplot, repeat_measure, binwidth=25, binrange=(0, 2000), stat="percent", alpha=0.5)

    # Add the median lines and annotations for each facet
    for ax, (anchor_name, subset) in zip(g.axes.flat, dataframe.groupby("anchor_name")):
        # Calculate the median for the subset of the data
        median_value = np.round(subset[repeat_measure].median())

        # Add a dashed line at the rounded median
        ax.axvline(median_value, color='black', linestyle='--', linewidth=2)

        # Annotate the median value in the middle-right corner of the facet
        ax.text(ax.get_xlim()[1] * 0.9, ax.get_ylim()[1] * 0.9, f'Median: {int(median_value)}',
                horizontalalignment='right', fontsize=12, fontweight='bold', color='black')

        # Adjust y-axis limits based on the maximum value of the histogram
        max_y = ax.get_ylim()[1]  # Get the current maximum y limit
        ax.set_ylim(0, max_y * 1.1)  # Increase y limit slightly for better visibility

    # Add title and labels to the overall plot
    g.set_titles("{col_name}")
    g.set_axis_labels(f'Telomere Repeat Length (bp)', "Frequency (%)")

    # Set the overall plot title
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(f'{title_name} (N = {num_reads_repeat})', fontweight="bold", fontsize=20, color='k')

    # Save the facet plot as an image file
    save_file_name = f'{figure_directory}/{title_name}_facet_by_chr_new.png'
    print(f'Graphing {save_file_name}...')
    g.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()

def histogram_plot_repeat_ref_y_prime_layer(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', hue='reference_y_prime_end_status', multiple='layer', common_norm=True, binwidth=10, binrange=(0,750),
                      stat="percent", palette='tab10', alpha=0.4, kde=True, hue_order=['Reference Y Primes = 0', 'Reference Y Primes = 1', 'Reference Y Primes = 2+'])
        
    ax.set_title(f'{title_name} Ref Y Primes (N = {num_reads_repeat_good})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)


    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=750)
    ax.set_ylim(bottom=0, top=20)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_layer_reference_y_prime_telomere_repeat.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()
        
def histogram_plot_repeat_read_y_prime_layer(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', hue='read_y_prime_end_status', multiple='layer', common_norm=True, binwidth=10, binrange=(0,750),
                      stat="percent", palette='Set2', alpha=0.4, kde=True, hue_order=['Read Y Primes = 0', 'Read Y Primes = 1', 'Read Y Primes = 2+'])

        
    ax.set_title(f'{title_name} Read Y Primes (N = {num_reads_repeat_good})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)


    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=750)
    ax.set_ylim(bottom=0, top=20)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_layer_read_y_prime_telomere_repeat.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_repeat_delta_y_prime_layer(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', hue='delta_y_prime_sign', multiple='layer', binwidth=10, binrange=(0,750), stat="percent", palette=['tab:red','tab:green'],
                      hue_order=['-', '+'])
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', hue='delta_y_prime_sign', multiple='layer', common_norm=True, binwidth=10, binrange=(0,750),
                      stat="percent", palette=['tab:red','tab:green'], alpha=0.4, kde=True, hue_order=['-', '+'])
        
    ax.set_title(f'{title_name} Delta Y Primes (N = {num_reads_repeat_good_delta_y_prime})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)


    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=750)
    ax.set_ylim(bottom=0, top=20)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_layer_y_prime_delta_telomere_repeat.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def histogram_plot_repeat_y_prime_type_delta_y_prime_layer(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', hue='y_prime_type_delta_y_prime_sign', multiple='layer', common_norm=True, binwidth=10, binrange=(0,1000),
                      stat="percent", palette='tab10', alpha=0.4, kde=True,
                      hue_order=['+ Reference Y Primes = 0', '- Reference Y Primes = 1', '+ Reference Y Primes = 1', '- Reference Y Primes = 2+', '+ Reference Y Primes = 2+'])
        
    ax.set_title(f'{title_name} Y Prime Type Delta Y Primes (N = {num_reads_repeat_good_delta_y_prime})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)


    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=1000)
    ax.set_ylim(bottom=0, top=20)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_layer_y_prime_type_delta_y_primes_telomere_repeat.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()


def histogram_plot_delta(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x='trimmed_delta_length_distance_to_ref_end', kde=True, stat="percent", color='forestgreen') # binwidth=10, binrange=(0,500),

        
    ax.set_title(f'{title_name} Delta to Ref End (N = {num_reads_delta})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel("Delta to Reference End", fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    #ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)

    # Set the bounds of the x and y axes
    #ax.set_xlim(left=0, right=500)
    #ax.set_ylim(bottom=0, top=12)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_delta_to_ref_end.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    #plt.close()
    #plt.show()
    
def histogram_plot_telo_loss_y_prime(dataframe, size=(16,11)) -> None:

    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Plot histogram
    g1 = sns.histplot(data=dataframe, x=f'{repeat_measure}', binwidth=20, binrange=(0,1000), kde=True, stat="percent", color='darkblue')

        
    ax.set_title(f'{title_name} with Y Prime Loss (N = {num_reads_telo_loss_y_prime})', fontweight="bold", fontsize = 20, color='k', pad=15) 
    ax.set_xlabel(f'Telomere Repeat Length (bp)', fontweight="bold", fontsize = 30, color='k', labelpad = 6)
    ax.set_ylabel("Frequency (%)", fontweight="bold", fontsize = 30, color='k', labelpad = 6)

    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)

    ax.set_yticks([2,4,6,8,10,12,14,16,18,20])
    ax.set_xticks([0,50,100,200,300,400,500,600,700,800,900,1000])

    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=1000)
    ax.set_ylim(bottom=0, top=20)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_1kb_with_y_prime_loss.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    #plt.close()
    #plt.show()


def ridgeplot_repeat_y_prime_type_delta_y_prime(dataframe, size=(16,11)) -> None:
    
    title_name = f'{base_name.split("_no_tag")[0]}'
    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    plot_order = ['same Reference Y Primes = 0', '+ Reference Y Primes = 0','- Reference Y Primes = 1', 'same Reference Y Primes = 1', '+ Reference Y Primes = 1', 
                  '- Reference Y Primes = 2+', 'same Reference Y Primes = 2+', '+ Reference Y Primes = 2+']

    g1 = sns.FacetGrid(data=dataframe, row="y_prime_type_delta_y_prime_sign", row_order=plot_order, aspect=15, height=1.5,)
    g1.map_dataframe(sns.kdeplot, x=f'{repeat_measure}', fill=True, hue='y_prime_type_delta_y_prime_sign', common_norm=False, alpha=0.2,
                     palette=['tab:gray', 'tab:green','tab:red', 'tab:gray', 'tab:green','tab:red', 'tab:gray', 'tab:green'], hue_order=plot_order, bw_adjust=0.5, cut=0)

    # Modify set_titles to include the count of items in each facet
    for i, sign_type in enumerate(plot_order):
        df_temp_sign_type = dataframe[dataframe['y_prime_type_delta_y_prime_sign'] == sign_type]
        
        count = len(df_temp_sign_type)
        g1.axes[i, 0].set_title(f'{sign_type} (N = {count})')

        # Add red line for median in each facet
        lower_05_val = df_temp_sign_type[repeat_measure].quantile(0.05)  # 5 Percentile
        q1 = df_temp_sign_type[repeat_measure].quantile(0.25)  # First quartile (Q1)
        median_val = df_temp_sign_type[repeat_measure].quantile(0.50)  # Second quartile (median) (Q2)
        q3 = df_temp_sign_type[repeat_measure].quantile(0.75)  # Third quartile (Q3)
        upper_95_val = df_temp_sign_type[repeat_measure].quantile(0.95)  # 95 Percentile
        
        g1.axes[i, 0].axvline(x=lower_05_val, color='red', linestyle='--', linewidth=2)
        g1.axes[i, 0].axvline(x=q1, color='black', linestyle='--', linewidth=2)
        g1.axes[i, 0].axvline(x=median_val, color='blue', linestyle='--', linewidth=2)
        g1.axes[i, 0].axvline(x=q3, color='black', linestyle='--', linewidth=2)
        g1.axes[i, 0].axvline(x=upper_95_val, color='red', linestyle='--', linewidth=2)
        
        g1.axes[i, 0].set_xlim(left=0, right=600)

    
    # Adjust subplot layout to create space for the title
    plt.subplots_adjust(top=0.85)
    g1.fig.suptitle(f'{title_name} Y Prime Type & Delta Y (N = {num_reads_repeat_good_ridge})', fontweight="bold", fontsize = 20, color='k') 


    # Adjust tick mark size and length of the hash marks
    ax.tick_params(axis='both', which='both', direction='out', length=8, width=2, labelsize=25, pad=2)


    # Set the bounds of the x and y axes
    ax.set_xlim(left=0, right=1000)
    ax.set_ylim(bottom=0, top=20)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_ridgeplot_y_prime_type_delta_y_primes_telomere_repeat.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()

def plot_y_prime_type_delta_gain_y_prime_barplot(dataframe, size=(16,11)):
    
    title_name = f'{base_name.split("_no_tag")[0]}'

    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    #### fix here
    # Calculate normalized counts for each day and outcome_status
    normalized_counts = dataframe.groupby(['reference_y_prime_end_status', 'y_primes_relative_to_ref']).size() / dataframe.groupby('reference_y_prime_end_status').size() * 100
    normalized_counts = normalized_counts.reset_index(name='percentage') 

    # Plot histogram
    g1 = sns.barplot(data=normalized_counts, x='reference_y_prime_end_status', y='percentage', hue='y_primes_relative_to_ref', palette='bright')
    
    plt.ylabel('Percentage of each Y Prime Gain Amount in a Read')
    plt.title(f'{title_name} Y Prime Gain by Ref Y Prime Type ({percentage_y_prime_gain_reads:.2f}% of Total Reads)')
    plt.legend(title='Delta Y Prime')
    
    # Get counts for each reference_y_prime_end_status
    counts_per_group = dataframe['reference_y_prime_end_status'].value_counts().to_dict()
    #print(len(dataframe))
    #print(counts_per_group)
    
    if len(counts_per_group) == 3:
        annotation_list = [(0,counts_per_group['Reference Y Primes = 0']), (1,counts_per_group['Reference Y Primes = 1']), (2,counts_per_group['Reference Y Primes = 2+'])]
    elif len(counts_per_group) == 2:
        if 'Reference Y Primes = 0' not in counts_per_group:
            annotation_list = [(0,counts_per_group['Reference Y Primes = 1']), (1,counts_per_group['Reference Y Primes = 2+'])]
        if 'Reference Y Primes = 1' not in counts_per_group:
            annotation_list = [(0,counts_per_group['Reference Y Primes = 0']), (1,counts_per_group['Reference Y Primes = 2+'])]    
        if 'Reference Y Primes = 2+' not in counts_per_group:
            annotation_list = [(0,counts_per_group['Reference Y Primes = 0']), (1,counts_per_group['Reference Y Primes = 1'])]
    elif len(counts_per_group) == 1:
        if 'Reference Y Primes = 0' in counts_per_group:
            annotation_list = [(0,counts_per_group['Reference Y Primes = 0'])]
        if 'Reference Y Primes = 1' in counts_per_group:
            annotation_list = [(0,counts_per_group['Reference Y Primes = 1'])]    
        if 'Reference Y Primes = 2+' in counts_per_group:
            annotation_list = [(0,counts_per_group['Reference Y Primes = 2+'])]
    else:
        (0,0)

    max_y_height = int(normalized_counts['percentage'].max())

    # Annotate bars with counts
    for index, count in annotation_list:
        g1.text(index, max_y_height, f"N: {count}", color='black', ha="center")

    # Set the bounds of the x and y axes
    #ax.set_xlim(left=0, right=500)
    #ax.set_ylim(bottom=0, top=12)
    
    # Set background color
    ax.set_facecolor('w')

    save_file_name = f'{figure_directory}/{title_name}_y_prime_gain_by_ref_y_prime_type.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()  
    
    

def determine_outcome_status(repeat_length, y_primes_relative_to_ref):           
    # Is rif1rif2    
    if file_name.startswith("dorado_7324"):
        short_telo_cutoff_length = 100
        long_telo_cutoff_length = 1000
    # Default telomere cutoffs
    else:
        short_telo_cutoff_length = 100
        long_telo_cutoff_length = 500
    
    if y_primes_relative_to_ref > 0:
        if repeat_length > long_telo_cutoff_length:
            status = 'gain_Y, long_telo'
        elif repeat_length < short_telo_cutoff_length:
            status = 'gain_Y, short_telo'
        else:
            status = 'gain_Y, normal_telo'
    elif y_primes_relative_to_ref < 0:
        if repeat_length > long_telo_cutoff_length:
            status = 'loss_Y, long_telo'
        elif repeat_length < short_telo_cutoff_length:
            status = 'loss_Y, short_telo'
        else:
            status = 'loss_Y, normal_telo'
    else: # y_primes_relative_to_ref == 0:
        if repeat_length > long_telo_cutoff_length:
            status = 'same_Y, long_telo'
        elif repeat_length < short_telo_cutoff_length:
            status = 'same_Y, short_telo'
        else:
            status = 'same_Y, normal_telo'
    return status

def plot_read_outcome_status(dataframe, size=(16,11)):
    
    title_name = f'{base_name.split("_no_tag")[0]}'

    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Calculate normalized counts for each day and outcome_status
    normalized_counts = dataframe.groupby(['day', 'outcome_status']).size() / dataframe.groupby('day').size() * 100
    normalized_counts = normalized_counts.reset_index(name='percentage')
    
    print(normalized_counts)
    
    day_list = normalized_counts['day'].unique()
    day_value_list = [int(day.strip('day')) for day in day_list]
    day_value_list = sorted(day_value_list)
    day_plot_order = [f'day{day}' for day in day_value_list]

    # Plot histogram
    g1 = sns.barplot(data=normalized_counts, x='day', y='percentage', hue='outcome_status', order=day_plot_order,
                        palette=['#fa758b','#e80008','#8a0005', '#98defa','#00b6ff','#00628a', '#75ff81','#00ff16','#008c0c'],
                        hue_order=['loss_Y, short_telo', 'loss_Y, normal_telo', 'loss_Y, long_telo',
                                'same_Y, short_telo', 'same_Y, normal_telo', 'same_Y, long_telo',
                                'gain_Y, short_telo', 'gain_Y, normal_telo', 'gain_Y, long_telo'])
    
    ax.set_ylim(bottom=0, top=100)
    
    plt.xlabel('Days')
    plt.ylabel('Percentage of each Day')
    plt.title(f'{title_name} Per Day Status for Each Day')
    plt.legend(title='Outcome Status')
    
    # Set background color
    ax.set_facecolor('w')


    save_file_name = f'{timecourse_figure_directory}/{title_name}_status_outcome_barplot.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()            

def determine_survivor_status(repeat_length, y_prime_probe_count, short_telo_cutoff_length, long_telo_cutoff_length):           
    
    # Short cutoff set by 10%-tile, and long cutoff set by 99%-tile + 100bp
    
    if y_prime_probe_count >= 2: # Type 1 or Hybrid
        if repeat_length < long_telo_cutoff_length:
            survivor_status = 'Type_1'
        else: survivor_status = 'Hybrid'
    else: # Type 2, Healthy, or Short
        if repeat_length >= long_telo_cutoff_length:
            survivor_status = 'Type_2'
        else:
            if repeat_length > short_telo_cutoff_length:
                survivor_status = 'Healthy'
            else: survivor_status = 'Short'

    return survivor_status

def plot_read_survivor_status(dataframe, size=(16,11)):
    
    title_name = f'{base_name.split("_no_tag")[0]}'

    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    # Calculate normalized counts for each day and survivor_status
    normalized_counts = dataframe.groupby(['day', 'survivor_status']).size() / dataframe.groupby('day').size() * 100
    normalized_counts = normalized_counts.reset_index(name='percentage')
    
    print(normalized_counts)
    
    day_list = normalized_counts['day'].unique()
    day_value_list = [int(day.strip('day')) for day in day_list]
    day_value_list = sorted(day_value_list)
    day_plot_order = [f'day{day}' for day in day_value_list]
    
    print(day_plot_order)

    # Plot histogram
    g1 = sns.barplot(data=normalized_counts, x='day', y='percentage', hue='survivor_status', order=day_plot_order,
                        palette=['tab:green','tab:red', 'tab:blue','tab:olive','tab:orange'],
                        hue_order=['Healthy', 'Short', 'Type_1', 'Type_2', 'Hybrid'])

    ax.set_ylim(bottom=0, top=100)

    plt.xlabel('Days')
    plt.ylabel('Percentage of each Day')
    plt.title(f'{title_name} Per Day Survivor Status')
    plt.legend(title='Survivor Status')
    
    # Set background color
    ax.set_facecolor('w')


    save_file_name = f'{timecourse_figure_directory}/{title_name}_survivor_status_barplot.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=300, format="png")
    plt.close()
    #plt.show()            



def timecourse_plot_ridgeplot_repeat_y_prime_type_delta_y_prime(dataframe, size=(16,11)):
    
    plot_order = ['same Reference Y Primes = 0', '+ Reference Y Primes = 0','- Reference Y Primes = 1', 'same Reference Y Primes = 1', '+ Reference Y Primes = 1', 
                  '- Reference Y Primes = 2+', 'same Reference Y Primes = 2+', '+ Reference Y Primes = 2+']

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))


    for plot_type in plot_order:
        
        
        df_plot_type = dataframe[dataframe['y_prime_type_delta_y_prime_sign'] == plot_type]
        
        title_name = f'{base_name.split("_no_tag")[0]}_{plot_type.replace(" ","_")}'

        # Set the graph background
        #custom_params = {}
        sns.set(style="ticks")#, rc=custom_params)

        # Set the figure size
        #fig, ax = plt.subplots(figsize=(size))

        day_list = df_plot_type['day'].unique()
        day_value_list = [int(day.strip('day')) for day in day_list]
        day_value_list = sorted(day_value_list)
        day_plot_order = [f'day{day}' for day in day_value_list]

        print(day_plot_order)
        g1 = sns.FacetGrid(data=df_plot_type, row="day", row_order=day_plot_order, aspect=15, height=1.5,)
        g1.map_dataframe(sns.kdeplot, x=f'{"repeat_length"}', fill=True, hue='day', common_norm=False, alpha=0.2,
                        palette='tab10', hue_order=day_plot_order, bw_adjust=0.5, cut=0)

        # Modify set_titles to include the count of items in each facet
        for i, day_num in enumerate(day_plot_order):
            df_temp_day_num = df_plot_type[df_plot_type['day'] == day_num]
            
            count = len(df_temp_day_num)
            g1.axes[i, 0].set_title(f'{day_num} (N = {count})')

            # Add red line for median in each facet
            lower_05_val = df_temp_day_num["repeat_length"].quantile(0.05)  # 5 Percentile
            q1 = df_temp_day_num["repeat_length"].quantile(0.25)  # First quartile (Q1)
            median_val = df_temp_day_num["repeat_length"].quantile(0.50)  # Second quartile (median) (Q2)
            q3 = df_temp_day_num["repeat_length"].quantile(0.75)  # Third quartile (Q3)
            upper_95_val = df_temp_day_num["repeat_length"].quantile(0.95)  # 95 Percentile
            
            
            g1.axes[i, 0].axvline(x=lower_05_val, color='red', linestyle='--', linewidth=2)
            g1.axes[i, 0].axvline(x=q1, color='black', linestyle='--', linewidth=2)
            g1.axes[i, 0].axvline(x=median_val, color='blue', linestyle='--', linewidth=2)
            g1.axes[i, 0].axvline(x=q3, color='black', linestyle='--', linewidth=2)
            g1.axes[i, 0].axvline(x=upper_95_val, color='red', linestyle='--', linewidth=2)
            
            
            # Annotate percentile values above the lines
            for x_val, percentile_val in zip([lower_05_val, q1, median_val, q3, upper_95_val],
                                            [0.05, 0.25, 0.50, 0.75, 0.95]):
                if i == 0:
                    g1.axes[i, 0].text(x_val, g1.axes[i, 0].get_ylim()[1] * 1.2, f'{percentile_val*100:.0f}%',
                                    ha='center', va='bottom', color='black', fontsize=10, fontweight='bold')                    
                g1.axes[i, 0].text(x_val, g1.axes[i, 0].get_ylim()[1] * 1, f'{x_val:.0f}',
                                ha='center', va='bottom', color='black', fontsize=10)
            
            
            g1.axes[i, 0].set_xlim(left=0, right=600)

            #print('stats')
            if plot_type != 'same Reference Y Primes = 0':
                if plot_type == '+ Reference Y Primes = 0':
                    comparison_to = 'same Reference Y Primes = 0'
                    df_comparison = dataframe[(dataframe['day'] == day_num) &
                                            (dataframe['y_prime_type_delta_y_prime_sign'] == comparison_to)]
                    null_hypothesis = 'less'  
                if plot_type in ['- Reference Y Primes = 1','+ Reference Y Primes = 1']:
                    comparison_to = 'same Reference Y Primes = 1'
                    df_comparison = dataframe[(dataframe['day'] == day_num) &
                                            (dataframe['y_prime_type_delta_y_prime_sign'] == comparison_to)]
                    null_hypothesis = 'less' if plot_type == '+ Reference Y Primes = 1' else 'greater'
                if plot_type in ['- Reference Y Primes = 2+','+ Reference Y Primes = 2+']:
                    comparison_to = 'same Reference Y Primes = 2+'
                    df_comparison = dataframe[(dataframe['day'] == day_num) &
                                            (dataframe['y_prime_type_delta_y_prime_sign'] == comparison_to)]         
                    null_hypothesis = 'less' if plot_type == '+ Reference Y Primes = 2+' else 'greater'
                if plot_type in ['same Reference Y Primes = 1','same Reference Y Primes = 2+']:
                    comparison_to = 'same Reference Y Primes = 0'
                    df_comparison = dataframe[(dataframe['day'] == day_num) &
                                            (dataframe['y_prime_type_delta_y_prime_sign'] == comparison_to)]         
                    null_hypothesis = 'less'
                    
                                
                null_hypothesis = 'two-sided'
                
                test_results = pg.mwu(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], alternative=null_hypothesis)
                cohens_d = pg.compute_effsize(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], paired=False, eftype='cohen')
                #hedges_g = pg.compute_effsize(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], paired=False, eftype='hedges')
                
                u1, p = mannwhitneyu(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], alternative=null_hypothesis)
                u2 = (len(df_temp_day_num["repeat_length"]) * len(df_comparison["repeat_length"])) - u1
                
                if u1 <= u2:
                    plot_description = "Lesser average rank than comparison"
                else:
                    plot_description = "Greater average rank than comparison"
                
                
                #print(test_results)
                #print(cohens_d)
                #print(hedges_g)
                
                """
                print(f'{day_num}: {plot_type} mean = {df_temp_day_num["repeat_length"].mean()},\
                      same Reference Y Primes = 0 mean = {df_comparison["repeat_length"].mean()}')
                
                print(f'{day_num}: {plot_type} median = {df_temp_day_num["repeat_length"].quantile(0.50)},\
                      same Reference Y Primes = 0 median = {df_comparison["repeat_length"].quantile(0.50)}')
                
                print(f'{day_num}: {plot_type} 95% = {df_temp_day_num["repeat_length"].quantile(0.95)},\
                      same Reference Y Primes = 0 95% = {df_comparison["repeat_length"].quantile(0.95)}')
                
                print(f'{day_num}: {plot_type} max = {df_temp_day_num["repeat_length"].max()},\
                      same Reference Y Primes = 0 max = {df_comparison["repeat_length"].max()}')                                
                  
                print(f'Cohens D: {cohens_d}')
                """               
                
                # Prints the results for each ridge, at the 500bp mark, taking the ymax times either 1, 0.5, or 0.3
                g1.axes[i, 0].text(500, g1.axes[i, 0].get_ylim()[1] * 1, f'Comparison to {comparison_to}:\n{test_results}', ha='center', va='center', color='black', fontsize=10)
                g1.axes[i, 0].text(500, g1.axes[i, 0].get_ylim()[1] * 0.5, f'{plot_description}; Cohens D: {cohens_d:.2f}', ha='center', va='center', color='black', fontsize=10)
            

        # Adjust subplot layout to create space for the title
        plt.subplots_adjust(top=0.85)
        g1.fig.suptitle(f'{title_name} (N = {len(df_plot_type)})', fontweight="bold", fontsize = 20, color='k') 

        # Set background color
        ax.set_facecolor('w')
        save_file_name = f'{timecourse_figure_directory}/{title_name}_ridgeplot_y_prime_type_telomere_repeat.png'
        print(f'Graphing {save_file_name}...')
        plt.savefig(f"{save_file_name}", dpi=100, format="png")
        plt.close()
        #plt.show()

def timecourse_plot_ridgeplot_read_y_prime_type_delta_y_prime(dataframe, size=(16,11)):
    
    plot_order = ['Read Y Primes = 0', 'Read Y Primes = 1','Read Y Primes = 2+']
    
    df_baseline = dataframe[dataframe['read_y_prime_end_status'] == 'Read Y Primes = 0']

    for plot_type in plot_order:
        
        df_plot_type = dataframe[dataframe['read_y_prime_end_status'] == plot_type]
        
        title_name = f'{base_name.split("_no_tag")[0]}_{plot_type.replace(" ","_")}'

        # Set the graph background
        #custom_params = {}
        sns.set(style="ticks")#, rc=custom_params)

        # Set the figure size
        fig, ax = plt.subplots(figsize=(size))

        day_list = df_plot_type['day'].unique()
        day_value_list = [int(day.strip('day')) for day in day_list]
        day_value_list = sorted(day_value_list)
        day_plot_order = [f'day{day}' for day in day_value_list]

        g1 = sns.FacetGrid(data=df_plot_type, row="day", row_order=day_plot_order, aspect=15, height=1.5,)
        g1.map_dataframe(sns.kdeplot, x=f'{"repeat_length"}', fill=True, hue='day', common_norm=False, alpha=0.2,
                        palette='tab10', hue_order=day_plot_order, bw_adjust=0.5, cut=0)

        # Modify set_titles to include the count of items in each facet
        for i, day_num in enumerate(day_plot_order):
            df_temp_day_num = df_plot_type[df_plot_type['day'] == day_num]
            
            count = len(df_temp_day_num)
            g1.axes[i, 0].set_title(f'{day_num} (N = {count})')

            # Add red line for median in each facet
            lower_05_val = df_temp_day_num["repeat_length"].quantile(0.05)  # 5 Percentile
            q1 = df_temp_day_num["repeat_length"].quantile(0.25)  # First quartile (Q1)
            median_val = df_temp_day_num["repeat_length"].quantile(0.50)  # Second quartile (median) (Q2)
            q3 = df_temp_day_num["repeat_length"].quantile(0.75)  # Third quartile (Q3)
            upper_95_val = df_temp_day_num["repeat_length"].quantile(0.95)  # 95 Percentile
            
            g1.axes[i, 0].axvline(x=lower_05_val, color='red', linestyle='--', linewidth=2)
            g1.axes[i, 0].axvline(x=q1, color='black', linestyle='--', linewidth=2)
            g1.axes[i, 0].axvline(x=median_val, color='blue', linestyle='--', linewidth=2)
            g1.axes[i, 0].axvline(x=q3, color='black', linestyle='--', linewidth=2)
            g1.axes[i, 0].axvline(x=upper_95_val, color='red', linestyle='--', linewidth=2)
            
            
            # Annotate percentile values above the lines
            for x_val, percentile_val in zip([lower_05_val, q1, median_val, q3, upper_95_val],
                                            [0.05, 0.25, 0.50, 0.75, 0.95]):
                if i == 0:
                    g1.axes[i, 0].text(x_val, g1.axes[i, 0].get_ylim()[1] * 1.2, f'{percentile_val*100:.0f}%',
                                    ha='center', va='bottom', color='black', fontsize=10, fontweight='bold')                    
                g1.axes[i, 0].text(x_val, g1.axes[i, 0].get_ylim()[1] * 1, f'{x_val:.0f}',
                                ha='center', va='bottom', color='black', fontsize=10)
                
            g1.axes[i, 0].set_xlim(left=0, right=600)
            
            if plot_type != 'Read Y Primes = 0':
                df_comparison = df_baseline[df_baseline['day'] == day_num]
                
                test_results = pg.mwu(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], alternative='two-sided')
                cohens_d = pg.compute_effsize(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], paired=False, eftype='cohen')
                #hedges_g = pg.compute_effsize(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], paired=False, eftype='hedges')
                
                u1, p = mannwhitneyu(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], alternative='two-sided')
                u2 = (len(df_temp_day_num["repeat_length"]) * len(df_comparison["repeat_length"])) - u1
                
                if u1 <= u2:
                    plot_description = "Lesser average rank than comparison"
                else:
                    plot_description = "Greater average rank than comparison"                
                
                #print(test_results)
                
                # Prints the results for each ridge, at the 500bp mark, taking the ymax times either 1, 0.5, or 0.3
                g1.axes[i, 0].text(500, g1.axes[i, 0].get_ylim()[1] * 1, f'Comparison to Read Y Primes = 0:\n{test_results}', ha='center', va='center', color='black', fontsize=10)
                g1.axes[i, 0].text(500, g1.axes[i, 0].get_ylim()[1] * 0.5, f'{plot_description}; Cohens D: {cohens_d:.2f}', ha='center', va='center', color='black', fontsize=10)

        
        # Adjust subplot layout to create space for the title
        plt.subplots_adjust(top=0.85)
        g1.fig.suptitle(f'{title_name} (N = {len(df_plot_type)})', fontweight="bold", fontsize = 20, color='k') 

        # Set background color
        ax.set_facecolor('w')
        save_file_name = f'{timecourse_figure_directory}/{title_name}_ridgeplot_y_prime_type_telomere_repeat.png'
        print(f'Graphing {save_file_name}...')
        plt.savefig(f"{save_file_name}", dpi=100, format="png")
        plt.close()
        #plt.show()

def timecourse_plot_ridgeplot_read_y_prime_type_delta_y_prime_overlaying(dataframe, size=(16,11)):
    
    title_name = f'{base_name.split("_no_tag")[0]}_timecourse_read_y_prime_repeat_overlay'

    # Set the graph background
    #custom_params = {}
    sns.set(style="ticks")#, rc=custom_params)

    # Set the figure size
    fig, ax = plt.subplots(figsize=(size))

    day_list = dataframe['day'].unique()
    day_value_list = [int(day.strip('day')) for day in day_list]
    day_value_list = sorted(day_value_list)
    day_plot_order = [f'day{day}' for day in day_value_list]

    g1 = sns.FacetGrid(data=dataframe, row="day", row_order=day_plot_order, aspect=15, height=1.5,)
    g1.map_dataframe(sns.kdeplot, x=f'{"repeat_length"}', fill=True, hue='read_y_prime_end_status', hue_order=['Read Y Primes = 0', 'Read Y Primes = 1', 'Read Y Primes = 2+'],
                     common_norm=False, alpha=0.1, linewidth=1.5, palette='Set2', bw_adjust=0.5, cut=0, legend=True)
    g1.add_legend(title='Read Y Prime End Status', label_order=['Green: Read Y Primes = 0', 'Orange: Read Y Primes = 1', 'Blue: Read Y Primes = 2+'])

    # Modify set_titles to include the count of items in each facet
    for i, day_num in enumerate(day_plot_order):
        df_temp_day_num = dataframe[dataframe['day'] == day_num]
        
        count = len(df_temp_day_num)
        y_prime_read_0_percent = len(df_temp_day_num[df_temp_day_num['read_y_prime_end_status'] == 'Read Y Primes = 0'])/count
        y_prime_read_1_percent = len(df_temp_day_num[df_temp_day_num['read_y_prime_end_status'] == 'Read Y Primes = 1'])/count
        y_prime_read_min_2_percent = len(df_temp_day_num[df_temp_day_num['read_y_prime_end_status'] == 'Read Y Primes = 2+'])/count
        
        g1.axes[i, 0].set_title(f"{day_num} (N = {count})\n" \
                                f"0={y_prime_read_0_percent*100:.0f}%, " \
                                f"1={y_prime_read_1_percent*100:.0f}%, " \
                                f"2+={y_prime_read_min_2_percent*100:.0f}%",
                                y=0.7)

        g1.axes[i, 0].set_xlim(left=0, right=600)

    # Adjust subplot layout to create space for the title
    plt.subplots_adjust(top=0.85)
    g1.fig.suptitle(f'{title_name}', fontweight="bold", fontsize = 20, color='k') 
    
    # Set background color
    ax.set_facecolor('w')
    save_file_name = f'{timecourse_figure_directory}/{title_name}.png'
    print(f'Graphing {save_file_name}...')
    plt.savefig(f"{save_file_name}", dpi=100, format="png")
    plt.close()
    #plt.show()



def compare_timecourse_plot_ridgeplot_read_y_prime_type_delta_y_prime(dataframe1, dataframe2, genotype1='wild_type', genotype2='mph1', size=(16,11)):
    
    plot_order = ['Read Y Primes = 0', 'Read Y Primes = 1','Read Y Primes = 2+']
    
    dataframe1['genotype'] = genotype1
    dataframe2['genotype'] = genotype2
    
    df_combined = pd.concat([dataframe1, dataframe2], ignore_index=True) 

    for plot_type in plot_order:
        print(plot_type)
        
        df_plot_type = df_combined[df_combined['read_y_prime_end_status'] == plot_type]
        
        title_name = f'{genotype1}_vs_{genotype2}_{plot_type.replace(" ","_")}'

        # Set the graph background
        #custom_params = {}
        sns.set(style="ticks")#, rc=custom_params)

        # Set the figure size
        fig, ax = plt.subplots(figsize=(size))

        day_list = df_plot_type['day'].unique()
        day_value_list = [int(day.strip('day')) for day in day_list]
        day_value_list = sorted(day_value_list)
        day_plot_order = [f'day{day}' for day in day_value_list]

        print(day_plot_order)
        g1 = sns.FacetGrid(data=df_plot_type, row="day", row_order=day_plot_order, aspect=15, height=1.5,)
        g1.map_dataframe(sns.kdeplot, x=f'{"repeat_length"}', fill=True, hue='genotype', hue_order=[genotype1, genotype2], common_norm=False, alpha=0.4,
                        palette=['#64d4fa', '#fc3df3'], bw_adjust=0.5, cut=0)
        g1.add_legend(title='Genotype', label_order=[f'Blue: {genotype1}', f'Pink: {genotype2}'])

        # Modify set_titles to include the count of items in each facet
        for i, day_num in enumerate(day_plot_order):
            df_temp_day_num = df_plot_type[df_plot_type['day'] == day_num]       
        
            count = len(df_temp_day_num)
            g1.axes[i, 0].set_title(f'{day_num} (N = {count})')
            
            for genotype in [genotype1, genotype2]:
                df_temp_day_num_genotype = df_temp_day_num[df_temp_day_num['genotype'] == genotype]
                
                # Add red line for median in each facet
                #lower_05_val = df_temp_day_num_genotype["repeat_length"].quantile(0.05)  # 5 Percentile
                #q1 = df_temp_day_num_genotype["repeat_length"].quantile(0.25)  # First quartile (Q1)
                median_val = df_temp_day_num_genotype["repeat_length"].quantile(0.50)  # Second quartile (median) (Q2)
                #q3 = df_temp_day_num_genotype["repeat_length"].quantile(0.75)  # Third quartile (Q3)
                #upper_95_val = df_temp_day_num_genotype["repeat_length"].quantile(0.95)  # 95 Percentile

                if genotype == genotype1:
                    g1.axes[i, 0].axvline(x=median_val, color='blue', linestyle='--', linewidth=2)
                else: g1.axes[i, 0].axvline(x=median_val, color='pink', linestyle='--', linewidth=2)
                
                # Annotate percentile values above the lines
                for x_val, percentile_val in zip([median_val],
                                                [0.50]):
                    if i == 0:
                        g1.axes[i, 0].text(x_val, g1.axes[i, 0].get_ylim()[1] * 1.2, f'{percentile_val*100:.0f}%',
                                        ha='center', va='bottom', color='black', fontsize=10, fontweight='bold')                    
                    g1.axes[i, 0].text(x_val, g1.axes[i, 0].get_ylim()[1] * 1, f'{x_val:.0f}',
                                    ha='center', va='bottom', color='black', fontsize=10)
                    
                g1.axes[i, 0].set_xlim(left=0, right=600)
            
            df_temp_day_num_genotype1 = df_temp_day_num[df_temp_day_num['genotype'] == genotype1]
            df_temp_day_num_genotype2 = df_temp_day_num[df_temp_day_num['genotype'] == genotype2]
            
            test_results = pg.mwu(df_temp_day_num_genotype1["repeat_length"], df_temp_day_num_genotype2["repeat_length"], alternative='two-sided')
            cohens_d = pg.compute_effsize(df_temp_day_num_genotype1["repeat_length"], df_temp_day_num_genotype2["repeat_length"], paired=False, eftype='cohen')
            #hedges_g = pg.compute_effsize(df_temp_day_num["repeat_length"], df_comparison["repeat_length"], paired=False, eftype='hedges')
            
            u1, p = mannwhitneyu(df_temp_day_num_genotype1["repeat_length"], df_temp_day_num_genotype2["repeat_length"], alternative='two-sided')
            u2 = (len(df_temp_day_num_genotype1["repeat_length"]) * len(df_temp_day_num_genotype2["repeat_length"])) - u1
            
            if u1 <= u2:
                plot_description = "Lesser average rank than comparison"
            else:
                plot_description = "Greater average rank than comparison"                
            
            # Prints the results for each ridge, at the 500bp mark, taking the ymax times either 1, 0.5, or 0.3
            g1.axes[i, 0].text(500, g1.axes[i, 0].get_ylim()[1] * 1, f'Comparison of {genotype1} to {genotype2}:\n{test_results}', ha='center', va='center', color='black', fontsize=10)
            g1.axes[i, 0].text(500, g1.axes[i, 0].get_ylim()[1] * 0.5, f'{plot_description}; Cohens D: {cohens_d:.2f}', ha='center', va='center', color='black', fontsize=10)

        
        # Adjust subplot layout to create space for the title
        plt.subplots_adjust(top=0.85)
        g1.fig.suptitle(f'{title_name} (N = {len(df_plot_type)})', fontweight="bold", fontsize = 20, color='k') 

        # Set background color
        ax.set_facecolor('w')
        save_file_name = f'{timecourse_figure_directory}/{title_name}_genotype_comparison_ridgeplot_y_prime_type_telomere_repeat.png'
        print(f'Graphing {save_file_name}...')
        plt.savefig(f"{save_file_name}", dpi=100, format="png")
        plt.close()
        #plt.show()


def timecourse_plot_ridgeplot_telomere_repeat_lengths(dataframe, size=(16,11)):
    
    for f_name in plot_files_to_analyze_list:
        if 'day10' in f_name: # it is a survivor
            plot_surviors = True
            break
        else:
            plot_surviors = False
    
    if plot_surviors: 
        
        title_name = f'{base_name.split("_no_tag")[0]}_survivor_telomere_repeat_lengths'

        print(dataframe['day'].unique())

        # Make a suvivor only dataframe
        dataframe['is_survivor'] = dataframe['day'].apply(lambda x: int(x.strip('day')) >= 10)
        dataframe = dataframe[dataframe['is_survivor'] == True]
        # Make a good survivor numbering scheme
        dataframe['IT_num'] = dataframe['day'].apply(lambda x: f"IT{x.strip('day')}")
        IT_num_order = sorted(dataframe['IT_num'].unique(), key=lambda x: int(x.strip('IT')))
        survivor_id_dict = {it_num: f's{i+1}' for i, it_num in enumerate(IT_num_order)}
        dataframe['survivor_id'] = dataframe['IT_num'].apply(lambda x: survivor_id_dict[x])
        
        survivor_id_plot_order = sorted(dataframe['survivor_id'].unique(), key=lambda x: int(x.strip('s')))
        #median = dataframe[dataframe['day'] == 'day0']['repeat_length'].median()

        # Set the graph background
        sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

        num_divisions = len(survivor_id_plot_order)
        gradient_colors = sns.color_palette("tab10", n_colors=num_divisions)
        color_map = dict(zip(survivor_id_plot_order, gradient_colors))
        print(color_map)
        
        g1 = sns.FacetGrid(dataframe, row="survivor_id", row_order=survivor_id_plot_order, hue="survivor_id", aspect=15, height=0.5, palette=color_map)
        
        g1.map(sns.kdeplot, "repeat_length", clip_on=False, fill=True, alpha=0.5, lw=1, cut=0, bw_adjust=0.5)
        #g1.map(sns.kdeplot, "repeat_length", clip_on=False, color="white", lw=1.5, cut=0, bw_adjust=0.5)
        #g1.map(plt.axvline, x=median, linestyle='--', color='k', linewidth=1)
        
        for ax, survivor_id_num, grad_color in zip(g1.axes.flat, survivor_id_plot_order, gradient_colors):
            median_val = dataframe[dataframe['survivor_id'] == survivor_id_num]["repeat_length"].median()
            ax.axvline(x=median_val, color=grad_color, linestyle='--', linewidth=2)
            #ax.axvline(x=median_val, color='k', linestyle='--', linewidth=2)
            print(grad_color)
            
        # Set the subplots to overlap
        g1.fig.subplots_adjust(hspace=-.40)

        # Define and use a simple function to label the plot in axes coordinates
        def label(x, color, label):
            ax = plt.gca()
            ax.text(0, 0.4, label, fontweight="bold", fontsize=24, color=color,
                    ha="left", va="center", transform=ax.transAxes)

        g1.map(label, "repeat_length")

        # Remove axes details that don't play well with overlap
        g1.set_titles("")
        g1.set(yticks=[], ylabel="")
        g1.despine(bottom=True, left=True)
        g1.set(xlim=(0, 1000))
        
        #set the titles
        plt.suptitle(f'{title_name}', fontweight="bold", fontsize = 20, color='k')
        g1.set_xlabels(label="Repeat Length (bp)", fontweight="bold", fontsize=32)
        g1.set_ylabels(label='', fontsize=0)
        g1.tick_params(axis='x', which='major', labelsize=24)
        
        # set the dimensions of the plot
        g1.fig.set_figwidth(size[0])
        g1.fig.set_figheight(size[1])
        g1.fig.patch.set_facecolor('white')
        
        # Save the figure
        save_file_name = f'{timecourse_figure_directory}/{title_name}.png'
        print(f'Graphing {save_file_name}...')
        plt.savefig(f"{save_file_name}", dpi=300, format="png")
        plt.close()
        #plt.show()
        
    else:
        
        title_name = f'{base_name.split("_no_tag")[0]}_timecourse_telomere_repeat_lengths'

        day_plot_order = sorted(dataframe['day'].unique(), key=lambda x: int(x.strip('day')))
        #median = dataframe[dataframe['day'] == 'day0']['repeat_length'].median()

        # Set the graph background
        sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
        #pal = sns.color_palette("flare", n_colors=len(day_plot_order), as_cmap=True)
        
        cmap = plt.get_cmap('flare')
        num_divisions = len(day_plot_order)
        gradient_colors = [cmap(i / num_divisions) for i in range(num_divisions)]
        color_map = dict(zip(day_plot_order, gradient_colors))
        print(color_map)
        
        g1 = sns.FacetGrid(dataframe, row="day", row_order=day_plot_order, hue="day", aspect=15, height=0.5, palette=color_map)
        
        g1.map(sns.kdeplot, "repeat_length", clip_on=False, fill=True, alpha=0.5, lw=1, cut=0, bw_adjust=0.5)
        #g1.map(sns.kdeplot, "repeat_length", clip_on=False, color="white", lw=1.5, cut=0, bw_adjust=0.5)
        #g1.map(plt.axvline, x=median, linestyle='--', color='k', linewidth=1)
        
        for ax, day_num, grad_color in zip(g1.axes.flat, day_plot_order, gradient_colors):
            median_val = dataframe[dataframe['day'] == day_num]["repeat_length"].median()
            #ax.axvline(x=median_val, color=grad_color, linestyle='--', linewidth=2)
            ax.axvline(x=median_val, color='k', linestyle='--', linewidth=2)
            print(grad_color)
            
        # Set the subplots to overlap
        g1.fig.subplots_adjust(hspace=-.40)

        # Define and use a simple function to label the plot in axes coordinates
        def label(x, color, label):
            ax = plt.gca()
            ax.text(0, 0.3, label, fontweight="bold", fontsize=24, color=color,
                    ha="left", va="center", transform=ax.transAxes)

        g1.map(label, "repeat_length")

        # Remove axes details that don't play well with overlap
        g1.set_titles("")
        g1.set(yticks=[], ylabel="")
        g1.despine(bottom=True, left=True)
        g1.set(xlim=(0, 500))
        
        #set the titles
        plt.suptitle(f'{title_name}', fontweight="bold", fontsize = 20, color='k')
        g1.set_xlabels(label="Repeat Length (bp)", fontweight="bold", fontsize=32)
        g1.set_ylabels(label='', fontsize=0)
        g1.tick_params(axis='x', which='major', labelsize=24)
        
        # set the dimensions of the plot
        g1.fig.set_figwidth(size[0])
        g1.fig.set_figheight(size[1])
        g1.fig.patch.set_facecolor('white')
        
        # Save the figure
        save_file_name = f'{timecourse_figure_directory}/{title_name}.png'
        print(f'Graphing {save_file_name}...')
        plt.savefig(f"{save_file_name}", dpi=300, format="png")
        plt.close()
        #plt.show()



files_in_dir = os.listdir('./')

#if input_base_name[0] == 'all':
    #files_to_analyze_list = [f for f in files_in_dir if f.endswith('post_y_prime_probe.tsv')]  # repeat_results.tsv, combined_results.tsv
#else:

files_to_analyze_list = [input_f_name]

for file_name in files_to_analyze_list:
    
    print(f'Starting {file_name}...')
    
    # Use the base_name from command line arguments
    base_name_for_file = base_name
    figure_directory = plot_output_dir
    os.makedirs(figure_directory, exist_ok=True)
    
    # SKIP TIME-COURSE SECTION - jump directly to single file analysis
    # Load TSV file into a Pandas DataFrame
    df_all = pd.read_csv(file_name, sep='\t')
    #files_to_analyze_list = input_base_name
files_to_analyze_list = [input_f_name]


for file_name in files_to_analyze_list:
    
    print(f'Starting {file_name}...') 
    
    base_name = f'{file_name.split("_post_y_prime_probe.tsv")[0]}'
    base_name = f'{base_name.split("dorado_")[1]}'    
    
    figure_directory = f'telomere_figures/{base_name}'
    if os.path.exists(figure_directory) == True:
        pass
    else:
        os.system(f'mkdir -p {figure_directory}')

    
    if "day0" in file_name:

        if "PromethION" in file_name:
            group_run_id = file_name.split("_PromethION")[0]
        elif "MinION" in file_name:
            group_run_id = file_name.split("_MinION")[0]
        else:
            group_run_id = file_name.split("_post_y_prime_probe.tsv")[0]
            #break

        # Check if "repeat" is in the timecourse name or not and then get the correct group
        if "repeat" in group_run_id:
            group_run_prefix = group_run_id.split("_day0_repeat")[0]
            group_run_repeat_id = f'repeat{group_run_id.split("repeat")[1]}'
        else:
            group_run_prefix = group_run_id.split("_day0")[0]

        group_files_to_analyze_list = [f for f in files_in_dir if f.startswith(group_run_prefix) and f.endswith('post_y_prime_probe.tsv')]

        ## Check if "no_tag" or "TeloTag" is in the timecourse name or not and then get the correct group
        #if "repeat" in file_name:
            #if "no_tag" in file_name:
                #plot_files_to_analyze_list = [f for f in group_files_to_analyze_list if group_run_repeat_id in f and "no_tag" in f]
                #group_run_timecourse_id = f'{group_run_prefix}_{group_run_repeat_id}_no_tag'
            #elif "TeloTag" in file_name:
                #plot_files_to_analyze_list = [f for f in group_files_to_analyze_list if group_run_repeat_id in f and "TeloTag" in f]
                #group_run_timecourse_id = f'{group_run_prefix}_{group_run_repeat_id}_TeloTag'
            #else:
                #plot_files_to_analyze_list = [f for f in group_files_to_analyze_list if group_run_repeat_id in f]
                #group_run_timecourse_id = f'{group_run_prefix}_{group_run_repeat_id}'
        #else: # "repeat" not in group_run_id:
            #if "no_tag" in file_name:
                #plot_files_to_analyze_list = [f for f in group_files_to_analyze_list if "repeat" not in f and "no_tag" in f]
                #group_run_timecourse_id = f'{group_run_prefix}_no_tag'
            #elif "TeloTag" in file_name:
                #plot_files_to_analyze_list = [f for f in group_files_to_analyze_list if "repeat" not in f and "TeloTag" in f]
                #group_run_timecourse_id = f'{group_run_prefix}_TeloTag'
            #else:
                #plot_files_to_analyze_list = [f for f in group_files_to_analyze_list if "repeat" not in f]
                #group_run_timecourse_id = f'{group_run_prefix}'

        #print(plot_files_to_analyze_list)
    
        #df_group_plot = pd.DataFrame()
        #for f_name in plot_files_to_analyze_list:
            #df_adding = pd.read_csv(f_name, sep='\t')
            #if "IT" in f_name:
                #day_id = f_name.split("_IT")[1]
                #day_id = f'day{day_id.split("-")[0]}'
                #df_adding["day"] = day_id
            #else:
                #df_adding["day"] = f_name.split("_")[2]
            #df_group_plot = pd.concat([df_group_plot, df_adding])

        #df_group_plot['outcome_status'] = df_group_plot.apply(lambda row: determine_outcome_status(row["repeat_length"], row["y_primes_relative_to_ref"]), axis=1)
        
        #df_group_plot_adpt = df_group_plot[df_group_plot['Adapter_After_Telomere'] == True]
        
        #timecourse_figure_directory = f'telomere_figures/timecourses/{group_run_timecourse_id}'
        #if os.path.exists(timecourse_figure_directory) == True:
            #pass
        #else:
            #os.system(f'mkdir -p {timecourse_figure_directory}')
        
        #plot_read_outcome_status(df_group_plot_adpt, size=(16,11))
        
        #df_group_plot_good = df_group_plot_adpt[(df_group_plot_adpt["repeat_length"] >= 30)]
        
        #df_group_plot_good_day0 = df_group_plot_good[df_group_plot_good['day'] == 'day0']
        
        #short_cutoff = 75   #100    #df_group_plot_good_day0["repeat_length"].quantile(0.10)
        #long_cutoff = df_group_plot_good_day0["repeat_length"].quantile(0.99) + 100     #500  #df_group_plot_good_day0["repeat_length"].quantile(0.99) + 100
        
        #print(f'short= {short_cutoff} and long= {long_cutoff}')
        
        #df_group_plot_good['survivor_status'] = df_group_plot_good.apply(lambda row: determine_survivor_status(row["repeat_length"], row["y_prime_probe_count"], 
                                                                                                               #short_cutoff, long_cutoff), axis=1)
        #plot_read_survivor_status(df_group_plot_good, size=(16,11))
        
        #timecourse_plot_ridgeplot_repeat_y_prime_type_delta_y_prime(df_group_plot_good, size=(16,11))
        #timecourse_plot_ridgeplot_read_y_prime_type_delta_y_prime(df_group_plot_good, size=(16,11))
        
        #timecourse_plot_ridgeplot_read_y_prime_type_delta_y_prime_overlaying(df_group_plot_good, size=(16,11))
        
        #timecourse_plot_ridgeplot_telomere_repeat_lengths(df_group_plot_good, size=(16,11))

    #continue

    # Load TSV file into a Pandas DataFrame with custom column headers
    df_all = pd.read_csv(file_name, sep='\t') 
    # ['qseqid', 'qlen', 'length', 'qstart', 'qend', 'anchor_name', 'slen', 'sstart', 'send', 'pident', 'bitscore', 'evalue', 'match_length', 'l_ended', 'alignment_direction', 'trimmed_overhang_past_ref_telo_start', 'correct_overhang', 'wanted_section', 'telomere_at_end']
    # ['read_id', 'total_read_length', 'read_bp_used_for_match', 'match_start_on_read', 'match_end_on_read', 'anchor_name', 'total_anchor_length', 'match_start_on_anchor', 'match_end_on_anchor', 'pident', 'bitscore', 'evalue', 'wanted_section_of_read', 'l_end_chr', 'trimmed_overhang_past_ref_telo_start', 'expected_wt_overhang_distance'])
    
    # Selects for only "good" ends
    df_good = df_all.copy()
    ############################################## df_good = df_good[df_good['anchor_name'] != 'chr4R_anchor'] ##############################################
    ############################################## df_good = df_good[df_good['anchor_name'] != 'chr12R_anchor'] ##############################################
    ############################################## df_good = df_good[df_good['anchor_name'] != 'chr6L_anchor'] ##############################################
    
    # Selects for only reads with repeat at the end
    df_good = df_good.dropna(subset=["repeat_length"])
    
    df_adpt = df_good[df_good['Adapter_After_Telomere'] == True].copy()
    df_adpt['end_type'] = df_adpt.apply(lambda x: 'tag', axis = 1)
    end_protection = 'adpt'

    repeat_measure_list = ['repeat_length'] #'trimmed_overhang_past_ref_telo_start', 

    for repeat_measure in repeat_measure_list:
        print(repeat_measure)
        print(df_adpt)

        df_graph = df_adpt.dropna(subset=[repeat_measure])

        num_reads_all = len(df_adpt)
        
        histogram_plot_all_close_zoom(df_adpt, size=(16,11))
        histogram_plot_all_close_zoom_25(df_adpt, size=(16,11))


        ######## Below are the new functions
        df_adpt_repeat_good = df_adpt[(df_adpt[repeat_measure] >= 30)]
        num_reads_repeat_good = len(df_adpt_repeat_good)
        
        
        histogram_plot_repeat_ref_y_prime_layer(df_adpt_repeat_good, size=(16,11))
        
        histogram_plot_repeat_read_y_prime_layer(df_adpt_repeat_good, size=(16,11))


        df_adpt_repeat_good_ridge = df_adpt_repeat_good[(df_adpt_repeat_good[repeat_measure] <= 1000)]
        num_reads_repeat_good_ridge = len(df_adpt_repeat_good_ridge)
        
        ridgeplot_repeat_y_prime_type_delta_y_prime(df_adpt_repeat_good_ridge)


        df_adpt_repeat_good_delta_y_prime = df_adpt_repeat_good[df_adpt_repeat_good['delta_y_prime_sign'] != 'same']
        
        num_reads_repeat_good_delta_y_prime = len(df_adpt_repeat_good_delta_y_prime)
        
        histogram_plot_repeat_delta_y_prime_layer(df_adpt_repeat_good_delta_y_prime)
        
        histogram_plot_repeat_y_prime_type_delta_y_prime_layer(df_adpt_repeat_good_delta_y_prime)
        
        
        df_delta_y_prime_gain = df_adpt_repeat_good[df_adpt_repeat_good['delta_y_prime_sign'] == '+']
        total_y_prime_gain_good_repeat_reads = len(df_delta_y_prime_gain)      

        df_delta_y_prime_gain = df_delta_y_prime_gain.sort_values(by=['y_primes_relative_to_ref'])
        
        percentage_y_prime_gain_reads = (total_y_prime_gain_good_repeat_reads/num_reads_repeat_good) * 100
        
        plot_y_prime_type_delta_gain_y_prime_barplot(df_delta_y_prime_gain)            
        

        df_adpt_repeat_short = df_adpt[(df_adpt[repeat_measure] <= 500) & (df_adpt[repeat_measure] >= 0)]
        num_reads_repeat = len(df_adpt_repeat_short)
        
        histogram_plot_repeat(df_adpt_repeat_short, size=(16,11))
        histogram_plot_repeat_chr(df_adpt_repeat_short, size=(16,11))
        histogram_plot_repeat_chr_facet(df_adpt_repeat_short, size=(16,11))        


        for chr_value in df_adpt_repeat_short['anchor_name'].unique():
            df_individual_chr_repeat_short = df_adpt_repeat_short[df_adpt_repeat_short['anchor_name'] == chr_value]
            
            genotype_chr_figure_directory = f'telomere_figures/chromosomes/{chr_value.split("_")[0]}'
            if os.path.exists(genotype_chr_figure_directory) == True:
                pass
            else:
                os.system(f'mkdir -p {genotype_chr_figure_directory}')
            
            
            num_individul_chr_reads_repeat = len(df_individual_chr_repeat_short)
            
            histogram_plot_individual_chr_repeat(df_individual_chr_repeat_short, size=(16,11), chr_anchor=chr_value)
        

        df_adpt_10kb = df_adpt[df_adpt[repeat_measure] <= 10000]
        num_reads_10kb = len(df_adpt_10kb)
        histogram_plot_10kb(df_adpt_10kb, size=(16,11))
        histogram_plot_10kb_zoom(df_adpt_10kb, size=(16,11))
        histogram_plot_10kb_close_zoom(df_adpt_10kb, size=(16,11))
        histogram_plot_10kb_close_zoom_chr(df_adpt_10kb, size=(16,11))
        

    has_tag = check_for_tag(base_name, df_good)

    df_tag_and_adpt = df_good[df_good['both_adapter_and_tag'] == True].copy()
    df_tag_and_adpt['end_type'] = df_tag_and_adpt.apply(lambda x: 'df_tag_and_adpt', axis =1)
    end_protection = 'tag_and_adpt'


    has_tag = check_for_tag(base_name, df_good)

    df_tag_and_adpt = df_good[df_good['both_adapter_and_tag'] == True].copy()
    df_tag_and_adpt['end_type'] = df_tag_and_adpt.apply(lambda x: 'df_tag_and_adpt', axis =1)
    end_protection = 'tag_and_adpt'

    if has_tag == True:
        for repeat_measure in repeat_measure_list:

            num_reads_all = len(df_tag_and_adpt)
            histogram_plot_all_close_zoom(df_tag_and_adpt, size=(16,11))
            histogram_plot_all_close_zoom_25(df_tag_and_adpt, size=(16,11))

            df_tag_and_adpt_repeat = df_tag_and_adpt[df_tag_and_adpt[repeat_measure] <= 500]
            num_reads_repeat = len(df_tag_and_adpt_repeat)
            histogram_plot_repeat(df_tag_and_adpt_repeat, size=(16,11))
            histogram_plot_repeat_chr(df_tag_and_adpt_repeat, size=(16,11))

            df_tag_and_adpt_10kb = df_tag_and_adpt[df_tag_and_adpt[repeat_measure] <= 10000]
            num_reads_10kb = len(df_tag_and_adpt_10kb)
            histogram_plot_10kb(df_tag_and_adpt_10kb, size=(16,11))
            histogram_plot_10kb_zoom(df_tag_and_adpt_10kb, size=(16,11))
            histogram_plot_10kb_close_zoom(df_tag_and_adpt_10kb, size=(16,11))
            histogram_plot_10kb_close_zoom_chr(df_tag_and_adpt_10kb, size=(16,11))


###############################################################################################
"""
Below is the skeleton for comparing tag vs no tag distributions

# Load TSV file into a Pandas DataFrame with custom column headers
df_tag = pd.read_csv(input_f_name, sep='\t') # , names=['qseqid', 'qlen', 'length', 'qstart', 'qend', 'anchor_name', 'slen', 'sstart', 'send', 'pident', 'bitscore', 'evalue', 'match_length', 'l_ended', 'alignment_direction', 'trimmed_overhang_past_ref_telo_start', 'correct_overhang', 'wanted_section', 'telomere_at_end']
df_tag = df_tag.drop_duplicates()
df_tag = df_tag[df_tag['anchor_name'] != 'chr4R_anchor']
df_tag = df_tag[df_tag['anchor_name'] != 'chr12R_anchor']
df_tag = df_tag[df_tag['anchor_name'] != 'chr6L_anchor']

df_tag_and_adpt = df_good[df_good['either_adapter_or_tag'] == True]
df_tag_and_adpt['end_type'] = 'tag'




histogram_plot_repeat_overlay(df_plot_repeat)

"""

