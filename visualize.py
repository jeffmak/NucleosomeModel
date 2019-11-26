import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set()
import pandas as pd
from pybdm import BDM

# read input sequence
with open("input.fa") as in_file:
    for i, line in enumerate(in_file):
        print(line)
        if i == 1:  # second line contains the DNA sequence
            seq = line[:-1]  # remove last '\n' character
## GC content
window_size = 147
# numpy implementation (handles edge of string)
base2GC = {'A': 0, 'T': 0, 'C': 1, 'G': 1}
GCseq = [base2GC[nt] for nt in seq]
GC_counts = np.convolve(GCseq, [1.0 / window_size] * window_size, 'same')
norm_GC_counts =  (GC_counts - np.min(GC_counts)) / (np.max(GC_counts) - np.min(GC_counts))
plt.plot(GC_counts, label='GC Content')
plt.plot(norm_GC_counts, label='norm GC Content')

# TODO: Block decomposition method
base2num = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
num_seq = np.array([base2num[nt] for nt in seq])
x_pos = []
bdm_values = []
bdm = BDM(ndim=1, nsymbols=4)
for i in range(len(seq) - window_size):
    x_pos.append(i + window_size // 2)
    substr = num_seq[i : i + window_size]
    bdm_values.append(bdm.bdm(substr))
bdm_values = np.array(bdm_values)
norm_bdm_values = (bdm_values - np.min(bdm_values)) / (np.max(bdm_values) - np.min(bdm_values))
plt.plot(x_pos, norm_bdm_values, label='normBDM')
print(x_pos)
print(norm_bdm_values)

# Kaplan Model
# read output file
print("reading tab-delimited file...")
data = pd.read_csv('out.tab', sep='\t')

# plot curves
print('plotting curves...')
data['P occupied'].plot(label='Kaplan Model')
plt.xlabel('Base Position')
plt.ylabel('Nucleosome Occupancy')
plt.ylim(0,1)
plt.legend()
plt.savefig('plots/plot.png')
