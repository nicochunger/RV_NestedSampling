# Python script to analyse all the results from the EPRV3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Results from the EPRV3 paper
eprv_broadZ = pd.read_csv('eprv3_broadZs.txt')
eprv_narrowZ = pd.read_csv('eprv3_narrowZ_norm.txt')

# Path to my results
filepath = '/home/nunger/tesis/resultados/eprv3/nlive2000_iter6/'
datafiles = ['0001', '0002', '0003', '0004', '0005', '0006']

# Initialization of pandas DataFrames for medians
medianZ_broad = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])
std_broad = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])
medianZ_narrow = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])
std_narrow = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])

# Broad Priors
for dfile in datafiles:
    medians = np.zeros(4)
    std = np.zeros(4)
    for nplanets in range(4):
        try:
            data = pd.read_csv(
                filepath + f'{dfile}/results{dfile}_{nplanets}a.txt', sep='\t')
            medians[nplanets] = data.median()['log10Z']
            std[nplanets] = data.std()['log10Z']
            # std[nplanets] = np.sqrt(data.median()[
            #                         'logZerr']**2 + np.median(np.abs(data['log10Z'].values - data.median()['log10Z']))**2)

        except:
            # That configuration hasn't been calculated yet
            # Just leave it at zero
            pass
    medianZ_broad[dfile] = medians
    std_broad[dfile] = std

# Narrow Priors
for dfile in datafiles:
    medians = np.zeros(4)
    std = np.zeros(4)
    for nplanets in range(4):
        try:
            data = pd.read_csv(
                filepath + f'{dfile}/results{dfile}_{nplanets}b.txt', sep='\t')
            medians[nplanets] = data.median()['log10Z']
            # medians[nplanets] = data.mean()['log10Z']
            std[nplanets] = data.std()['log10Z']
            # std[nplanets] = np.sqrt(data.median()[
            #                         'logZerr']**2 + np.median(np.abs(data['log10Z'].values - data.median()['log10Z']))**2)
        except:
            # That configuration hasn't been calculated yet
            pass
    medianZ_narrow[dfile] = medians
    std_narrow[dfile] = std

print('\nMy results of Broad Prior:')
print(medianZ_broad)

print('\nMy results of Narrow Prior:')
print(medianZ_narrow)

# Calculate difference
dif_broad = medianZ_broad - eprv_broadZ
dif_broad.rename_axis('Planets', axis=0, inplace=True)
dif_broad.rename_axis('Dataset', axis=1, inplace=True)
print('\nDifference of eprv res with my res: Broad Prior')
print(dif_broad)

# Calculate difference
dif_narrow = medianZ_narrow - eprv_narrowZ
dif_narrow.rename_axis('Planets', axis=0, inplace=True)
dif_narrow.rename_axis('Dataset', axis=1, inplace=True)
print('\nDifference of eprv res with my res: Broad Prior')
print(dif_narrow)

# Generar archivos con mismo formato que eprv3
columns = ['comptime', 'nplanets', 'mode_log10Z',
           'median_log10Z', '-2sigma', '-1sigma', '+1sigma', '+2sigma']
dframes = []
br_nr = ['a', 'b']
savepath = 'paper_results/PolyChord2/'
for typ in br_nr:
    if typ == 'a':
        Zdata = medianZ_broad
        stddata = std_broad
    else:
        Zdata = medianZ_narrow
        stddata = std_narrow
    for dset in datafiles:
        filetosave = pd.DataFrame()
        filetosave[columns[0]] = [1000]*4
        filetosave[columns[1]] = range(4)
        filetosave[columns[2]] = Zdata[dset]
        filetosave[columns[3]] = Zdata[dset]
        filetosave[columns[4]] = Zdata[dset] - 2*stddata[dset]
        filetosave[columns[5]] = Zdata[dset] - 1*stddata[dset]
        filetosave[columns[6]] = Zdata[dset] + 1*stddata[dset]
        filetosave[columns[7]] = Zdata[dset] + 2*stddata[dset]
        filetosave.to_csv(
            savepath+f'evidences_{dset}{typ}.txt', index=False, header=False)


# ---------------------- GRFICOS -----------------------------
# Broad
plt.figure()
plt.hlines(0, 0.5, 6.5, linestyles='--', linewidth=1)
plt.vlines(np.arange(0.5, 6, 1), -5, 5, linestyles='-', linewidth=1)
plt.errorbar(
    np.arange(1, 7, 1)-0.15, dif_broad.loc[0], yerr=std_broad.loc[0], fmt='.', label='0 Planets')
plt.errorbar(
    np.arange(1, 7, 1)-0.05, dif_broad.loc[1], yerr=std_broad.loc[1], fmt='.', label='1 Planet')
plt.errorbar(
    np.arange(1, 7, 1)+0.05, dif_broad.loc[2], yerr=std_broad.loc[2], fmt='.', label='2 Planets')
plt.errorbar(
    np.arange(1, 7, 1)+0.15, dif_broad.loc[3], yerr=std_broad.loc[3], fmt='.', label='3 Planets')
plt.xlim([0.5, 6.5])
plt.ylim([-5, 5])
plt.title('Broad Prior')
plt.xlabel('Dataset')
plt.ylabel(r'log(Z) - <log(Z)>')
plt.legend()

# Narrow
plt.figure()
plt.hlines(0, 0.5, 6.5, linestyles='--', linewidth=1)
plt.vlines(np.arange(0.5, 6, 1), -5, 5, linestyles='-', linewidth=1)
plt.errorbar(
    np.arange(1, 7, 1)-0.15, dif_narrow.loc[0], yerr=std_narrow.loc[0], fmt='.', label='0 Planets')
plt.errorbar(
    np.arange(1, 7, 1)-0.05, dif_narrow.loc[1], yerr=std_narrow.loc[1], fmt='.', label='1 Planet')
plt.errorbar(
    np.arange(1, 7, 1)+0.05, dif_narrow.loc[2], yerr=std_narrow.loc[2], fmt='.', label='2 Planets')
plt.errorbar(
    np.arange(1, 7, 1)+0.15, dif_narrow.loc[3], yerr=std_narrow.loc[3], fmt='.', label='3 Planets')
plt.xlim([0.5, 6.5])
plt.ylim([-5, 5])
plt.title('Narrow Prior')
plt.xlabel('Dataset')
plt.ylabel(r'log(Z) - <log(Z)>')
plt.legend()

plt.show()
# -----------------------------------------------------------
