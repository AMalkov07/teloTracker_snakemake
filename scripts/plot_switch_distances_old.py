import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# wt = 39417
wt = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 144, 213, 286, 400, 455, 455, 457, 490, 523,
      607, 624, 655, 694, 751, 775, 801, 908, 921, 937, 1018, 1023, 1024, 1024, 1039, 1045, 1046, 1245, 1266, 1607, 1935, 2266, 2274, 2519, 2521, 3171, 3487]
wt_r1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 142, 340, 429, 1025, 1025, 1398, 2705]
wt_r2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 278, 397, 1019, 1024, 1526, 2522]

mph_7172 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            204, 292, 459, 556, 601, 693, 696, 776, 1330, 1338, 1525, 1525, 1525, 1527, 1863, 1916, 2024, 2068, 2268, 2320]
mph_7302 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 221, 319, 754, 1370, 1559, 1690, 2518]

all_wt = wt + wt_r1 + wt_r2

all_mph = mph_7172 + mph_7302

df_wt = pd.DataFrame(data={'switch_distance': all_wt, 'genotype': 'WT'})
df_mph = pd.DataFrame(data={'switch_distance': all_mph, 'genotype': 'MPH'})



df_switches = pd.concat([df_wt, df_mph])

print(df_switches['genotype'].value_counts())

df_switches['switch_distance_rev'] = df_switches['switch_distance'] * -1

df_no_0 = df_switches.loc[df_switches['switch_distance'] > 0]
print(df_no_0['genotype'].value_counts())

figure_directory = 'spacer_outputs/switch_location_figures'

sns.histplot(data=df_switches, x='switch_distance_rev', hue='genotype', cumulative=True, stat="density", common_norm=False, kde=True, binwidth=50)

save_file_name = f'cumulative_density_curve_of_switch_location.png'
plt.title(f"{save_file_name.removesuffix('.png')}")
print(f'Graphing {figure_directory}/{save_file_name}...')
plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
plt.show()


sns.histplot(data=df_switches, x='switch_distance_rev', hue='genotype', cumulative=False, stat="percent", common_norm=False, kde=True, binwidth=250)

save_file_name = f'frequency_curve_of_switch_location.png'
plt.title(f"{save_file_name.removesuffix('.png')}")
print(f'Graphing {figure_directory}/{save_file_name}...')
plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
plt.show()


sns.histplot(data=df_no_0, x='switch_distance_rev', hue='genotype', cumulative=True, stat="density", common_norm=False, kde=True, binwidth=50)

save_file_name = f'cumulative_density_curve_no_0_of_switch_location.png'
plt.title(f"{save_file_name.removesuffix('.png')}")
print(f'Graphing {figure_directory}/{save_file_name}...')
plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
plt.show()


sns.histplot(data=df_no_0, x='switch_distance_rev', hue='genotype', cumulative=False, stat="percent", common_norm=False, kde=True, binwidth=250)

save_file_name = f'frequency_curve_no_0_of_switch_location.png'
plt.title(f"{save_file_name.removesuffix('.png')}")
print(f'Graphing {figure_directory}/{save_file_name}...')
plt.savefig(f"{figure_directory}/{save_file_name}", dpi=300, format="png")
plt.show()