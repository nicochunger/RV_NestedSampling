import numpy as np
import scipy.special as sp
from pprint import pprint
import pickle
import os

# Load evidence results
output = pickle.load(open('output.p', 'rb'))
print('Evidence (logZ) = {}'.format(output.logZ))
# Change direcory of posterior
output.base_dir = os.path.dirname(os.path.abspath(__file__))
posterior = output.posterior
samples = posterior.samples

# Get the medians for each parameter
paramnames = posterior.getParamNames().list()
medians = dict(zip(paramnames, zip(np.round(np.median(
    samples, axis=0), 5), np.round(np.std(samples, axis=0), 7))))

logZs = np.array(output.logZs)
planet1_medians = {}
for i in range(1, 50):
    cluster_posterior = output.cluster_posterior(i)
    cluster_samples = cluster_posterior.samples
    medians = dict(zip(paramnames, zip(np.round(np.median(
        cluster_samples, axis=0), 3), np.round(np.std(cluster_samples, axis=0), 3))))
    planet1_medians[i] = [medians['planet1_period'], logZs[-i]]

pprint(planet1_medians)

# Find which clusters correspong to 51 days and 131 days
idx51 = []
idx131 = []
idx098 = []
for key, value in planet1_medians.items():
    if value[0][0] >= 50 and value[0][0] <= 53:
        idx51.append(-1*key)
    if value[0][0] >= 130 and value[0][0] <= 133:
        idx131.append(-1*key)
    if value[0][0] >= 0.8 and value[0][0] <= 1.2:
        idx098.append(-1*key)

# Find evidences for those peaks and add them in log
logZ51 = sp.logsumexp(np.take(logZs, idx51))
print('(Log)Evidence for 51 days peak: {}'.format(logZ51))
logZ131 = sp.logsumexp(np.take(logZs, idx131))
print('(Log)Evidence for 131 days peak: {}'.format(logZ131))
logZ098 = sp.logsumexp(np.take(logZs, idx098))
print('(Log)Evidence for 131 days peak: {}'.format(logZ098))
