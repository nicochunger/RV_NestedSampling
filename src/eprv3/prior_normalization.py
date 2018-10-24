import numpy as np
import pandas as pd
import priors

datasets = ['0001', '0002', '0003', '0004', '0005', '0006']
narrow_priors = {
    1: [(39.8107, 44.6684), (11.4815, 12.8825), (10.0, 10.7452)],
    2: [(15.4882, 16.2181), (14.7911, 17.0608), (158.489, 251.189)],
    3: [(81.2831, 107.152), (38.0189, 42.658), (16.5959, 17.5792)],
    4: [(138.038, 204.174), (15.1356, 16.5959), (398.107, 1000.0)],
    5: [(29.5121, 32.3594), (10.7152, 11.4815), (18.197, 19.9526)],
    6: [(79.4328, 141.254), (31.6228, 50.1187), (316.228, 398.107)]}

narrowZ_notnorm = pd.read_csv('eprv3_narrowZ.txt')
narrowZ_norm = pd.read_csv('eprv3_narrowZ.txt')

input_dict = {'harps': {'period': [0, 1, ['Jeffreys', 1.25, 1e4], 0.13]}}

prior_dict = priors.prior_constructor(input_dict, {})
prior = prior_dict['harps_period']

# x = np.linspace(1.25, 1e4, 1000)
# pdf = prior.pdf(x)

for i, dset in enumerate(datasets):
    integrales = np.ones(4)
    for n in range(3):
        x = np.linspace(narrow_priors[i+1][n][0],
                        narrow_priors[i+1][n][1], 1000)
        y = prior.pdf(x)
        integrales[n] = np.trapz(y, x)
        narrowZ_norm[dset][n+1] = narrowZ_notnorm[dset][n+1] - \
            np.log10(np.prod(integrales))

narrowZ_norm.to_csv('eprv3_narrowZ_norm.txt', index=False)

print('\nNo normalizado:')
print(narrowZ_notnorm)

print('\nNormalizado:')
print(narrowZ_norm)
