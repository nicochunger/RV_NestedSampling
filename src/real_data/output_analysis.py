import pickle
import matplotlib.pyplot as plt
from PyAstronomy.pyasl import foldAt
import pandas as pd
import numpy as np
from pprint import pprint

datapath = '/home/nunger/tesis/codigo/data/HD40307/'
data = pd.read_csv(datapath+'HD40307_HARPS03(DRS-3-5)_night_binning.rdb',
                   sep='\t', comment='#', skiprows=[1, ])

t = data['rjd']
y = data['vrad']
# Correct offset
y = y - 31334.45
ey = data['svrad']
# rhk = data['rhk']

# # Plot of raw RV data
# plt.figure(0)
# plt.errorbar(t, y, yerr=ey, fmt='k.')

# Load evidence results
output = pickle.load(open('output.p', 'rb'))
print(output.logZ)
posterior = output.posterior
samples = posterior.samples

# Get the medians for each parameter
paramnames = posterior.getParamNames().list()
medians = dict(zip(paramnames, np.median(samples, axis=0)))
pprint(medians)

# Count the number of planets in model
nplanets = 0
for i in range(10):
    if f'planet{i}_period' in paramnames:
        nplanets += 1

# Construct dictionary with posterior for each parameter
post = {}
for i in range(len(paramnames)):
    post.update({paramnames[i]: samples[:, i]})

# Plot posterior for the period of each planet
for i in range(nplanets):
    plt.figure(i)
    # plt.hist(post[f'planet{i+1}_period'], bins='fd')
    plt.plot(post[f'planet{i+1}_period'], '.')
    plt.title(f"Planet {i+1}. Period: {medians[f'planet{i+1}_period']}")
plt.show()


# # Phase folding
# phases = np.zeros((len(t), nplanets))
# for i in range(nplanets):
#     phases[:, i] = foldAt(t, medians[f'planet{i+1}_period'])

# # Sort with respect to phase
# # First, get the order of indices ...
# sortIndi = np.argsort(phases, axis=0)
# # ... and, second, rearrange the arrays.
# y_folded = np.zeros((len(t), nplanets))
# for i in range(nplanets):
#     phases[:, i] = phases[:, i][sortIndi[:, i]]
#     y_folded[:, i] = y[sortIndi[:, i]]

#     # Plot the result
#     plt.figure(i+1)
#     plt.plot(phases[:, i], y_folded[:, i], '.')
#     plt.title(f"Planet {i}. Period: {medians[f'planet{i+1}_period']}")
# plt.show()
