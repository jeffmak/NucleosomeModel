import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set()
import pandas as pd
from pybdm import BDM
import math

# read input sequence
# input.fa - chr10 359283-360283 (length 1001)
print('Reading DNA sequence...')
with open("input.fa") as in_file:
    for i, line in enumerate(in_file):
        print(line)
        if i == 1:  # second line contains the DNA sequence
            seq = line[:-1]  # remove last '\n' character

# read raw signals from file
chr10_nuc_signals = pd.read_csv('yeast/chr10_processed', sep='\t', names=['Position', 'Signal Intensity'])
chr10_nuc_signals['Position'] -= 359283  # start position

interp_signals = chr10_nuc_signals.copy()
sig_inten = interp_signals['Signal Intensity']
interp_signals['Signal Intensity'] = (sig_inten - sig_inten.min()) / (sig_inten.max() - sig_inten.min())

norm_signals = chr10_nuc_signals.copy()
sig_inten = norm_signals['Signal Intensity']
norm_signals['Signal Intensity'] = (sig_inten - sig_inten.mean()) / sig_inten.std()

## GC content
print('Evaluating GC content...')
window_size = 147
# numpy implementation (handles edge of string)
base2GC = {'A': 0, 'T': 0, 'C': 1, 'G': 1}
GCseq = [base2GC[nt] for nt in seq]
GC_counts = np.convolve(GCseq, [1.0 / window_size] * window_size, 'same')
interp_GC_content = (GC_counts - np.min(GC_counts)) / (np.max(GC_counts) - np.min(GC_counts))
norm_GC_content = (GC_counts - np.mean(GC_counts)) / np.std(GC_counts)
# plt.plot(GC_counts, label='GC Content')


## Block decomposition method
print("Evaluating using Block Decomposition Method...")
base2num = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
raw_num_seq = [base2num[nt] for nt in seq]
# flank the input sequence with sequences of low algorithmic complexity
# flanked_num_seq = [2] * (window_size // 2)
# flanked_num_seq.extend(raw_num_seq)
# flanked_num_seq.extend([2] * (window_size // 2))
num_seq = np.array(raw_num_seq)

bdm_x_pos = []
bdm_values = []
bdm = BDM(ndim=1, nsymbols=4)
for i in range(len(seq) - window_size):
    bdm_x_pos.append(i + window_size // 2)
    substr = num_seq[i: i + window_size]
    bdm_values.append(bdm.bdm(substr))
bdm_values = np.array(bdm_values)
interp_bdm_values = (bdm_values - np.min(bdm_values)) / (np.max(bdm_values) - np.min(bdm_values))
norm_bdm_values = (bdm_values - np.mean(bdm_values)) / np.std(bdm_values)

# TODO: Van der Heijden Model
print("Evaluating using Van der Heijden Model...")
# parameters
B = 0.2  # amplitude
p = 10.1  # period
N = 147  # length bound to tetramer/octamer

heijden_x_pos = []
heijden_values = []
for i in range(len(seq) - N):
    heijden_x_pos.append(i + N // 2)
    substr = seq[i: i + N + 1]  # +1 for computing GC dinucleotides
    heijden_value = 4.0 ** N
    for j in range(N):
        s = j - N // 2
        delta_GC = 1 if substr[j:j + 2] == 'GC' else 0
        heijden_value *= 0.25 + B * math.cos(2 * math.pi * (s / p + delta_GC / 2))
    heijden_values.append(heijden_value)
interp_heijden_values = (heijden_values - np.min(heijden_values)) / (np.max(heijden_values) - np.min(heijden_values))
norm_heijden_values = (heijden_values - np.mean(heijden_values)) / np.std(heijden_values)

## Kaplan Model
print("Evaluating using Kaplan Model...")
# read output file
data = pd.read_csv('out.tab', sep='\t')
occup_col = data['P occupied']
data['P occupied (interp)'] = (occup_col - occup_col.min()) / (occup_col.max() - occup_col.min())
data['P occupied (norm)'] = (occup_col - occup_col.mean()) / (occup_col.std())

# plot curves
print('plotting curves...')
# interpolated plot
sns.lineplot(x='Position', y='Signal Intensity', data=interp_signals, label='Nuc Signal')
plt.plot(interp_GC_content, label='GC Content')
plt.plot(bdm_x_pos, interp_bdm_values, label='BDM')
plt.plot(heijden_x_pos, interp_heijden_values, label='Heijden')
data['P occupied (interp)'].plot(label='Kaplan Model')
plt.xlabel('Base Position')
plt.ylabel('Interpolated Nucleosome Occupancy Score')
plt.ylim(0, 1)
plt.legend()
plt.savefig('plots/interp_plot.png')
plt.close()

# norm plot
sns.lineplot(x='Position', y='Signal Intensity', data=norm_signals, label='Nuc Signal')
plt.plot(norm_GC_content, label='GC Content')
plt.plot(bdm_x_pos, norm_bdm_values, label='BDM')
plt.plot(heijden_x_pos, norm_heijden_values, label='Heijden')
data['P occupied (norm)'].plot(label='Kaplan Model')
plt.xlabel('Base Position')
plt.ylabel('Norm Nucleosome Occupancy Score')
plt.ylim(-3, 3)
plt.legend()
plt.savefig('plots/norm_plot.png')
plt.close()
