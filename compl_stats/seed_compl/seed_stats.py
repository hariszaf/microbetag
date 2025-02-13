# %% load libraries
import os, sys
import gzip
import pickle
import numpy as np
import matplotlib.pyplot as plt

# %% load files
with gzip.open("seeds_binary_per_patric.pckl.gz", "rb") as f:
    seeds = pickle.load(f)

with gzip.open("non_seeds_binary_per_patric.pckl.gz", "rb") as f:
    self_prod_mets = pickle.load(f)


# %%
cpd_intervals = seeds.sum()


# Create custom bins with a separate bin for zero
# Separate the zero values from the rest
zero_values = cpd_intervals[cpd_intervals == 0]
non_zero_values = cpd_intervals[cpd_intervals != 0]

# Create histogram for non-zero values with 30 bins
counts, bins, patches = plt.hist(non_zero_values, bins=30, edgecolor='black', alpha=0.6, label="Histogram", orientation='vertical')

# Count the number of zeros
zero_count = len(zero_values)

# Plot the zero bin manually
plt.bar(0, zero_count, width=0.1, edgecolor='black', alpha=0.6, color='blue', label=f"Zero count: {zero_count}")

# Calculate the cumulative sum of the non-zero histogram (cumulative frequency)
cumulative_counts = np.cumsum(counts)

# Set a threshold for the cumulative mass (e.g., 60% of the total sum)
threshold = 0.6 * cumulative_counts[-1]  # 60% of the total mass

# Ensure threshold is within the range of cumulative counts
threshold = min(threshold, cumulative_counts[-1])

# Find the bin that crosses the threshold
threshold_bin = np.searchsorted(cumulative_counts, threshold)

# Highlight the region up to the threshold
# plt.fill_between(bins[1:threshold_bin + 1], 0, counts[:threshold_bin + 1], color='orange', alpha=0.4, label=f'Up to threshold ({threshold:.2f})')

# Plot the cumulative histogram (second layer)
plt.plot(bins[1:], cumulative_counts, color='red', label="Cumulative", linewidth=2)

# Plot a vertical line for the threshold
plt.axvline(x=bins[threshold_bin], color='green', linestyle='--', label=f"Threshold at {bins[threshold_bin]:.2f}")

# Customize the plot
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title("Histogram with Cumulative Mass Threshold")
plt.legend()
plt.show()


#%%
seeds.columns[(seeds == 1).all()]
#  cpd01772 (Succinylbenzoate, 4)
#  cpd00345 (5-Methyltetrahydrofolate, 4)
#  cpd00557 (Siroheme, 4)

# No metabolites were found as being not a seed globally
seeds.columns[(seeds == 0).all()]



