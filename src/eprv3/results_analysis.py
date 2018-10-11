# Python script to analyse all the results from the EPRV3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Results from the EPRV3 paper
eprv_broadZ = pd.read_csv('eprv3_broadZs.txt')
eprv_narrowZ = pd.read_csv('eprv3_narrowZ.txt')

# Path to my results
filepath = '/home/nunger/tesis/resultados/eprv3/'
datafiles = ['0001', '0002', '0003', '0004', '0005', '0006']

# Initialization of pandas DataFrames for medians
medianZ_broad = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])
medianZ_narrow = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])

# Broad Priors
for dfile in datafiles:
    medians = np.zeros(4)
    for nplanets in range(4):
        try:
            data = pd.read_csv(
                filepath + f'{dfile}/results{dfile}_{nplanets}a.txt', sep='\t')
            # medians[nplanets] = data.median()['log10Z']
            medians[nplanets] = data.mean()['log10Z']
        except:
            # That configuration hasn't been calculated yet
            # Just leave it at zero
            pass
    medianZ_broad[dfile] = medians

# Narrow Priors
for dfile in datafiles:
    medians = np.zeros(4)
    for nplanets in range(4):
        try:
            data = pd.read_csv(
                filepath + f'{dfile}/results{dfile}_{nplanets}b.txt', sep='\t')
            # medians[nplanets] = data.median()['log10Z']
            medians[nplanets] = data.mean()['log10Z']
        except:
            # That configuration hasn't been calculated yet
            pass
    medianZ_narrow[dfile] = medians

print('\nBroad Prior:')
print(medianZ_broad)

print('\nNarrow Prior:')
print(medianZ_narrow)

# Calculate difference
dif_broad = medianZ_broad - eprv_broadZ
dif_broad.rename_axis('Planets', axis=0, inplace=True)
dif_broad.rename_axis('Dataset', axis=1, inplace=True)
print('\nDifference of eprv res with my res: Broad Prior')
print(dif_broad)

plt.figure()
plt.hlines(0, 0.5, 6.5)
plt.vlines(np.arange(0.5, 6, 1), -4, 4, linestyles='dashed')
plt.plot(range(1, 7), dif_broad.loc[0], '.', label='0 Planet')
plt.plot(range(1, 7), dif_broad.loc[1], '.', label='1 Planets')
plt.plot(range(1, 7), dif_broad.loc[2], '.', label='2 Planets')
plt.plot(range(1, 7), dif_broad.loc[3], '.', label='3 Planets')
plt.xlim([0.5, 6.5])
plt.ylim([-5, 5])
plt.xlabel('Dataset')
plt.ylabel(r'log(Z) - <log(Z)>')
plt.legend()
plt.show()
