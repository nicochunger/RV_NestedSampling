import pickle
import matplotlib.pyplot as plt
from PyAstronomy.pyasl import foldAt
import numpy as np
import pandas as pd
from pprint import pprint

datapath = '/home/nunger/tesis/codigo/data/HD40307/'
data = pd.read_csv(datapath+'HD40307_HARPS03(DRS-3-5)_night_binning.rdb',
                   sep='\t', comment='#', skiprows=[1, ])

t = data['rjd']
y = data['vrad']
# Correct offset
y = y - 31334.45
ey = data['svrad']
rhk = data['rhk']


# # Plot of raw RV data
# plt.figure(0)
# plt.errorbar(t, y, yerr=ey, fmt='k.')

# Load evidence results
output = pickle.load(open('output.p', 'rb'))
posterior = output.posterior
samples = posterior.samples

# Get the medians for each parameter
paramnames = posterior.getParamNames().list()
medians = dict(zip(paramnames, np.median(samples, axis=0)))
pprint(medians)

# Phase folding
nplanets = int((len(paramnames)-2)/5)
phases = np.zeros((len(t), nplanets))
for i in range(nplanets):
    phases[:, i] = foldAt(t, medians[f'planet{i+1}_period'])

# Sort with respect to phase
# First, get the order of indices ...
sortIndi = np.argsort(phases, axis=0)
# ... and, second, rearrange the arrays.
y_folded = np.zeros((len(t), nplanets))
for i in range(nplanets):
    phases[:, i] = phases[:, i][sortIndi[:, i]]
    y_folded[:, i] = y[sortIndi[:, i]]

    # Plot the result
    plt.figure(i+1)
    plt.plot(phases[:, i], y_folded[:, i], '.')
    plt.title(f"Planet {i}. Period: {medians[f'planet{i+1}_period']}")
plt.show()
