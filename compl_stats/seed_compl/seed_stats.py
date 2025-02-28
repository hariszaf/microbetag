# %% load libraries
import os, json
import gzip
import pickle
import numpy as np
import pandas as pd
import patchworklib as pw
from plotnine import ggplot, aes, element_text, labs, theme, geom_bar, geom_histogram

# %% load files
with gzip.open("seeds_binary_per_patric.pckl.gz", "rb") as f:
    seeds = pickle.load(f)

with gzip.open("non_seeds_binary_per_patric.pckl.gz", "rb") as f:
    self_prod_mets = pickle.load(f)

"""
452 unique metabolites as seeds
5454 unique metabolites as non-seeds
33755 genomes in total
"""

# %% Arrays
seed_intervals     = seeds.sum()
nonseeds_intervals = self_prod_mets.sum()


# Create custom bins with a separate bin for zero
# Separate the zero values from the rest
zero_values     = seed_intervals[seed_intervals == 0]
non_zero_values = seed_intervals[seed_intervals != 0]

#%% Distribution of the number of seeds per genome

seeds_num_per_genome = seeds.sum(axis=1)
df           = pd.DataFrame(seeds_num_per_genome)
df.columns   = ["number_of_seeds"]
counts, bins = np.histogram(df, bins=25)
df           = pd.DataFrame({'bins': bins[:-1], 'counts': counts})

p1 = (
    ggplot(df, aes(x='bins', y='counts')) +
    geom_bar(stat='identity', fill='#28579E', alpha=0.7) +
    labs(
        title="A. Histogram of number of seeds per genome",
        x="Number of seeds",
        y="Number of genomes") +
    theme(
        plot_title=element_text(size=16, weight="bold"),
        axis_text=element_text(size=14),   # Adjust axis tick labels font size
        axis_title=element_text(size=14)   # Adjust axis title font size
    )
)

seeds_per_genomes_plt = p1

#%% Distribution of the number of non-seeds per genome

nonseeds_num_per_genome = self_prod_mets.sum(axis=1)

df           = pd.DataFrame(nonseeds_num_per_genome, columns = ["number_of_nonseeds"])
counts, bins = np.histogram(df, bins=25)
df           = pd.DataFrame({'bins': bins[:-1], 'counts': counts})

p2 = (
    ggplot(df, aes(x='bins', y='counts')) +
    geom_bar(stat='identity', fill='#28579E', alpha=0.7) +
    labs(
        title="B. Histogram of number of non-seeds per genome",
        x="Number of non-seeds",
        y="Number of genomes") +
    theme(
        plot_title=element_text(size=16, weight="bold"),
        axis_text=element_text(size=14),   # Adjust axis tick labels font size
        axis_title=element_text(size=14)   # Adjust axis title font size
    )
)

nonseeds_per_genomes_plt = p2
"""
124 min - 442 max
"""

# -----------------------
# %% Seeds ratios
# -----------------------

df           = pd.DataFrame( seed_intervals / seeds.shape[0], columns = ["ratio"])
counts, bins = np.histogram(df, bins=30)
df           = pd.DataFrame({'bins': bins[:-1], 'counts': counts})

# Create the histogram plot
p3 = (
    ggplot(df, aes(x='bins', y='counts')) +
    geom_bar(stat='identity', fill='#28579E', alpha=0.7) +
    labs(
        x='Percentage of genomes in which metabolite is a seed',
        y='Number of metabolites',
        title='C. Distribution of seed compounds coverage'
    ) +
    theme(
        plot_title=element_text(size=16, weight="bold"),
        axis_text=element_text(size=14),   # Adjust axis tick labels font size
        axis_title=element_text(size=14)   # Adjust axis title font size
    )
)

seeds_coverage_plt = p3


# -----------------------
# %% Non-seeds ratios
# -----------------------

df           = pd.DataFrame( nonseeds_intervals / self_prod_mets.shape[0], columns = ["ratio"])
counts, bins = np.histogram(df, bins=30)
df           = pd.DataFrame({'bins': bins[:-1], 'counts': counts})

# Create the histogram plot
p4 = (
    ggplot(df, aes(x='bins', y='counts')) +
    geom_bar(stat='identity', fill='#28579E', alpha=0.7) +
    labs(
        x='Percentage of genomes in which metabolite is a non-seed',
        y='Number of metabolites',
        title='D. Distribution of non-seed compounds coverage'
    ) +
    theme(
        plot_title=element_text(size=16, weight="bold"),
        axis_text=element_text(size=14),   # Adjust axis tick labels font size
        axis_title=element_text(size=14)   # Adjust axis title font size
    )
)

nonseeds_coverage_plt = p4


# %% Number of genomes (donors) potentially providing a seed

if not os.path.exists("counts_of_a_genomes_potential_nonseed_hits.json"):

    nonseed_hits_per_genome = {}
    for index, genome in enumerate(self_prod_mets.index):

        # potential donor
        # self_prod_mets.iloc[1,:]   gets a row
        # (self_prod_mets.iloc[1,:][self_prod_mets.iloc[1,:]==1]  gets colums where row is 1
        # gets its column names
        compounds_of_interest = list(self_prod_mets.iloc[index,:][self_prod_mets.iloc[index,:]==1].index)

        compounds_of_interest_presnt = [x for x in compounds_of_interest if x in seeds.columns]

        # times a metabolite of donor's nonseeds appears as a seed accross the 33K genomes as a seed
        count_ones_specific = seeds[compounds_of_interest_presnt]

        c = count_ones_specific.sum().sum()

        nonseed_hits_per_genome[genome] = c

    data = {k: int(v) if isinstance(v, np.integer) else v for k, v in nonseed_hits_per_genome.items()}

    with open("counts_of_a_genomes_potential_nonseed_hits.json", "w") as f:
        json.dump(data, f)

else:
    with open("counts_of_a_genomes_potential_nonseed_hits.json", "r") as f:
        data = json.load(f)


# %%

df = pd.DataFrame.from_dict([data])
df = df.T

total_potential_hits = seeds.shape[0] * seeds.shape[1]

df_perc = df / total_potential_hits * 100
df_perc.columns = ["percentage"]

counts, bins = np.histogram(df_perc, bins=30)

df = pd.DataFrame({'bins': bins[:-1], 'counts': counts})

q = (
    ggplot(df, aes(x='bins', y='counts')) +
    geom_bar(stat='identity', fill='#28579E', alpha=0.7) +
    labs(
        title="E. Per genome non seeds overlaps across genomes' seeds",
        x="Percentage of potentially total nonseed overlap",
        y="Number of genomes") +
    theme(
        plot_title=element_text(size=16, weight="bold"),
        axis_text=element_text(size=14),   # Adjust axis tick labels font size
        axis_title=element_text(size=14)   # Adjust axis title font size
    )
)

nonseed_across_seeds_percentage_plt = q

# %% Build merged image

ax1 = pw.load_ggplot(p1)  # + labs(title='A')  this would overwrite the initial plot title
ax2 = pw.load_ggplot(p2)
ax3 = pw.load_ggplot(p3)
ax4 = pw.load_ggplot(p4)
ax5 = pw.load_ggplot(q)


# %%
ax123 = ax1/ax2|ax3/ax4|ax5
ax123.savefig("seeds_stats.png")











# %%


# ratios = {}
# count_non_compl = 0
# for seed, number_of_genomes_with_the_seed in seed_intervals.items():

#     try:
#         number_of_genomes_with_seed_as_nonseed = nonseeds_intervals[seed]
#         ratios[seed] = number_of_genomes_with_the_seed / number_of_genomes_with_seed_as_nonseed
#     except:
#         count_non_compl += 1
#         pass

# print(f"Number of seeds for which no genome was found to be able to complete them: {count_non_compl}")

# ratios_s = pd.Series(ratios)
# ratios_df = pd.DataFrame(ratios_s)
# ratios_df.index.name = "cpd"
# ratios_df.columns = ["ratio"]

# # ratios_df['ratio'] = ratios_df['ratio'] / ratios_df['ratio'].max()

# ratios_df_up_to_1 = ratios_df[ratios_df['ratio'] <=1]
# ratios_between_1_and_100 = ratios_df[(ratios_df['ratio'] >1) & (ratios_df['ratio']  < 100)]
# ratios_df_higher_than_100 = ratios_df[ratios_df['ratio'] >= 100 ]


# # %%
# p2 = (
#     ggplot(ratios_df_up_to_1, aes(x="ratio"))
#     + geom_histogram(bins=50, alpha=0.7, fill="#28579E")

#     + labs(
#         title="Histogram of number of seeds per genome",
#         x="Number of seeds",
#         y="Number of genomes")
#     + theme(
#         plot_title=element_text(size=16, weight="bold"),
#         axis_text=element_text(size=14),   # Adjust axis tick labels font size
#         axis_title=element_text(size=14)   # Adjust axis title font size
#     )
# )

# p2

# # %% High seed

# from plotnine import geom_line, geom_point, element_blank

# ratios_df['ratio'] = ratios_df['ratio'] / ratios_df['ratio'].max()

# # Create the plot
# p = (
#     ggplot(ratios_df_higher_than_1, aes(x=ratios_df_higher_than_1.index, y='ratio')) +
#     geom_line(color='blue') +  # Line plot
#     geom_point(color='red') +  # Add points
#     labs(x='Index', y='Values', title='1-Col DataFrame Plot') +
#     theme()
# )

# p







# # %%  Genomes a met was found a seed in

# # Create histogram for non-zero values with 30 bins
# counts, bins = np.histogram(non_zero_values, bins=30)

# # Create a DataFrame with bins and frequencies
# df = pd.DataFrame({'bins': bins[:-1], 'counts': counts})

# # threshold_value = bins[threshold_bin]

# # Create the histogram plot
# p2 = (
#     ggplot(df, aes(x='bins', y='counts')) +
#     geom_bar(stat='identity', fill='#28579E', alpha=0.7) +
#     # geom_vline(xintercept=threshold_value, linetype='dashed', color='green')+
#     labs(
#         x='Number of genomes in which metabolite is a seed',
#         y='Number of metabolites',
#         title='Genomes a metabolite was found a seed in'
#     ) +
#     theme(
#         plot_title=element_text(size=16, weight="bold"),
#         axis_text=element_text(size=14),   # Adjust axis tick labels font size
#         axis_title=element_text(size=14)   # Adjust axis title font size
#     )
# )

# p2



# # %%
# ratios_df['ratio'] = ratios_df['ratio'] / ratios_df['ratio'].max()







# deprecated
# ---------

# # Count the number of zeros
# zero_count = len(zero_values)

# # Plot the zero bin manually
# plt.bar(0, zero_count, width=0.1, edgecolor='black', alpha=0.6, color='blue', label=f"Zero count: {zero_count}")

# # Calculate the cumulative sum of the non-zero histogram (cumulative frequency)
# cumulative_counts = np.cumsum(counts)

# # Set a threshold for the cumulative mass (e.g., 60% of the total sum)
# threshold = 0.7 * cumulative_counts[-1]  # 60% of the total mass

# # Ensure threshold is within the range of cumulative counts
# threshold = min(threshold, cumulative_counts[-1])

# # Find the bin that crosses the threshold
# threshold_bin = np.searchsorted(cumulative_counts, threshold)

# # Highlight the region up to the threshold
# # plt.fill_between(bins[1:threshold_bin + 1], 0, counts[:threshold_bin + 1], color='orange', alpha=0.4, label=f'Up to threshold ({threshold:.2f})')

# # Plot the cumulative histogram (second layer)
# plt.plot(bins[1:], cumulative_counts, color='red', label="Cumulative", linewidth=2)

# # Plot a vertical line for the threshold
# plt.axvline(x=bins[threshold_bin], color='green', linestyle='--', label=f"Threshold at {bins[threshold_bin]:.2f}")

# # Customize the plot
# plt.xlabel("Value")
# plt.ylabel("Frequency")
# plt.title("Histogram with Cumulative Mass Threshold")
# plt.legend()
# plt.show()
