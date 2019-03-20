import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Path to my results
nlive2000_filepath = '/home/nunger/tesis/resultados/eprv3/nlive2000_iter6/'
nlive25nd_filepath = '/home/nunger/tesis/resultados/eprv3/nlive25*d_iter50/'
datafiles = ['0001', '0002', '0003', '0004', '0005', '0006']

# Initialization of pandas DataFrames for medians
medianZ_2000 = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])
std_2000 = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])

medianZ_25nd = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])
std_25nd = pd.DataFrame(columns=datafiles, index=[0, 1, 2, 3])

# Case with nlive = 2000
for dfile in datafiles:
    medians = np.zeros(4)
    std = np.zeros(4)
    for nplanets in range(4):
        try:
            data = pd.read_csv(
                nlive2000_filepath + f'{dfile}/results{dfile}_{nplanets}a.txt', sep='\t')
            medians[nplanets] = data.median()['log10Z']
            std[nplanets] = data.std()['log10Z']
            # std[nplanets] = np.sqrt(data.median()[
            #                         'logZerr']**2 + np.median(np.abs(data['log10Z'].values - data.median()['log10Z']))**2)

        except:
            # That configuration hasn't been calculated yet
            # Just leave it at zero
            pass
    medianZ_2000[dfile] = medians
    std_2000[dfile] = std


# Case with nlive = 25*ndims
for dfile in datafiles:
    medians = np.zeros(4)
    std = np.zeros(4)
    for nplanets in range(4):
        try:
            data = pd.read_csv(
                nlive25nd_filepath + f'{dfile}/results{dfile}_{nplanets}a.txt', sep='\t')
            medians[nplanets] = data.median()['log10Z']
            std[nplanets] = data.std()['log10Z']
            # std[nplanets] = np.sqrt(data.median()[
            #                         'logZerr']**2 + np.median(np.abs(data['log10Z'].values - data.median()['log10Z']))**2)

        except:
            # That configuration hasn't been calculated yet
            # Just leave it at zero
            pass
    medianZ_25nd[dfile] = medians
    std_25nd[dfile] = std

# print("\nnlive = 2000")
# print(medianZ_2000)

# print("\nnlive = 25*ndims")
# print(medianZ_25nd)

oddsfactor_2000_nelson = pd.DataFrame(
    columns=datafiles, index=["1v0", "2v1", "3v2"])
oddsfactor_2000_nelson.name = "Odds Factors for nlive=2000 with nelsons prior."

oddsfactor_2000_uniform = pd.DataFrame(
    columns=datafiles, index=["1v0", "2v1", "3v2"])

oddsfactor_25nd_nelson = pd.DataFrame(
    columns=datafiles, index=["1v0", "2v1", "3v2"])

oddsfactor_25nd_uniform = pd.DataFrame(
    columns=datafiles, index=["1v0", "2v1", "3v2"])


# Calculate all the odds factor
nelson = [oddsfactor_2000_nelson, oddsfactor_25nd_nelson]
uniform = [oddsfactor_2000_uniform, oddsfactor_25nd_uniform]
medians = [medianZ_2000, medianZ_25nd]

# With uniform prior
for i, nlive in enumerate(uniform):
    for datafile in datafiles:
        for j in range(3):
            odds = np.exp(medians[i][datafile][j+1]) / \
                np.exp(medians[i][datafile][j])
            nlive[datafile][j] = round(np.log10(odds), 3)

# With nelson prior


def nelson_prior(n):
    beta = 1/3
    if n == 0:
        return (1 - sum([beta**i for i in range(1, 4)]))
    else:
        return beta**n


for i, nlive in enumerate(nelson):
    for datafile in datafiles:
        for j in range(3):
            prior_factor = nelson_prior(j+1) / nelson_prior(j)
            odds = prior_factor * (np.exp(medians[i][datafile][j+1]) /
                                   np.exp(medians[i][datafile][j]))
            nlive[datafile][j] = round(np.log10(odds), 3)

print("\n Uniform prior")
print(uniform[0])
print("")
print(uniform[1])

print("\n Nelson prior")
print(nelson[0])
print("")
print(nelson[1])
