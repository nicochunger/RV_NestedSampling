import pickle
import matplotlib.pyplot as plt
# from PyAstronomy.pyasl import foldAt
import pandas as pd
import numpy as np
import numpy.ma as ma
from pprint import pprint
from scipy.optimize import curve_fit
import argparse
import os

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fit', action='store_true',
                    help='If posterior should be fit with a gaussian.')
args_params = parser.parse_args()

# Load evidence results
output = pickle.load(open('output.p', 'rb'))
print(f'Evidence (logZ) = {output.logZ}')
# Change direcory of posterior
output.base_dir = os.path.dirname(os.path.abspath(__file__))
posterior = output.posterior
samples = posterior.samples

# Get the medians for each parameter
paramnames = posterior.getParamNames().list()
medians = dict(zip(paramnames, zip(np.round(np.median(
    samples, axis=0), 5), np.round(np.std(samples, axis=0), 7))))
pprint(medians)

# --------------- POSTERIORS -----------------------------
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
fig, ax = plt.subplots(nplanets, 1)
ax = np.atleast_1d(ax)  # To support 1 planet models

# Separate samples for the 4th planet for 51 days and 131 days
period_post = post['planet1_period']
p51 = np.where((period_post >= 51.1) & (period_post <= 52))[0]
print(len(p51))
p131 = np.where((period_post >= 130.6) & (period_post <= 131.5))[0]
print(len(p131))

parnames51 = [par+'51' for par in paramnames]
medians_51 = dict(zip(parnames51, zip(np.round(np.median(
    samples[p51, :], axis=0), 5), np.round(np.std(samples[p51, :], axis=0), 7))))

parnames131 = [par+'131' for par in paramnames]
medians_131 = dict(zip(parnames131, zip(np.round(np.median(
    samples[p131, :], axis=0), 5), np.round(np.std(samples[p131, :], axis=0), 7))))

medians2 = {}
medians2.update(medians_51)
medians2.update(medians_131)

pprint(medians2)

# print('\nValues at 51 days:')
# pprint(medians_51)

# print('\nValues at 131 days')
# pprint(medians_131)
