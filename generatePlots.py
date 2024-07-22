# python3 generatePlots.py recombigator_output_on_X

import os
import matplotlib.pyplot as plt
import sys

base_dir = sys.argv[1]

# Create directory for plots
plots_dir = base_dir + "/plots"
if not os.path.exists(plots_dir):
    os.mkdir(plots_dir)


### Plot 1: Read Lengths
len_dir = base_dir + "/misc/read_lengths"
readLengths = {}

# read in all of the read lengths
filenames = os.listdir(len_dir)
filenames = sorted(filenames)
filenames[7], filenames[8] = filenames[8], filenames[7]

for filename in filenames:
    lengths_file = open(len_dir + "/" + filename, "r")
    lengths = [float(line.strip()) for line in lengths_file.readlines()]
    category_name = filename.split('_')[0]
    readLengths[category_name] = lengths

# plot the read lengths
fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 8))

# colors = ['black','darkred', 'firebrick','indianred','forestgreen','olivedrab','fuchsia','orchid','mediumvioletred','lightseagreen']
colors = ['black','darkred', 'firebrick','indianred','forestgreen','olivedrab','fuchsia','orchid','mediumvioletred']

labels = list(readLengths.keys())
for i, label in enumerate(labels):
    labels[i] = label + " (n = " + str(len(readLengths[label])) + ")"

axes.hist(readLengths.values(), bins = 20, label = labels, histtype = 'bar', color = colors)
axes.legend(prop = {'size' : 10})
axes.set_xlabel("Read Length (bp)")
axes.set_ylabel("Number of Reads")
axes.set_title("Read Lengths of different categories")
axes.set_yscale('log')
plt.tight_layout()
plt.savefig(plots_dir + "/sizes_plot.png")


### Plot 2: Crossover positions along the reads
info_dir = base_dir + "/recombinant_info"

# read in all of the information
hq_info_file = open(info_dir + "/" + "HighQualityRecombinantInfo.tsv", "r")
hq_info_file.readline()
hq_breakpoints = [float((line.strip()).split()[3]) for line in hq_info_file.readlines()]

n = len(hq_breakpoints)

# plot the crossover positions along the reads
fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 8))
axes.hist(hq_breakpoints, bins = 20, histtype = 'bar', color='green', edgecolor='black')
axes.set_xlabel("Fraction along SNPs where crossover occurs")
axes.set_ylabel("Number of Reads")
axes.set_title("Distribution of Crossover Positions on High Quality Recombinant Reads (n = " + str(n) + ")")
plt.savefig(plots_dir + "/hqrecombinant_crossover.png")