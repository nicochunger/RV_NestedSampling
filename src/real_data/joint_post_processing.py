#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python imports
import os
import subprocess
import argparse
import pickle

# Dependencies
import numpy as np
import matplotlib.pyplot as plt

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-k', type=int, default=3,
                    help='Maximum number of planets to analyse.')
parser.add_argument('-prior', type=str, default='u',
                    help='Prior to use for number of planets. Default: Uniform')
args_params = parser.parse_args()

# Maximum number of planets
k = args_params.k

# Loading of all posteriors
dirname = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), 'posteriors/')
# Get name of folders
runs = subprocess.check_output(
    'ls -d '+dirname+'*/', shell=True).decode('utf-8').replace(dirname, '').split('/\n')
runs.remove('')

posteriors = []  # All posteriors and evidences are going to be stored here
for i in range(k+1):
    # Dictionary to save all relevant data of this planet model
    postdict = {'k': i}
    folderpath = runs[i]  # Folder with posterios
    # Load PolyChord output
    output = pickle.load(open(dirname+folderpath+'/output.p', 'rb'))
    output.base_dir = dirname+folderpath  # Change direcory of posteriors
    posterior = output.posterior
    postdict['output'] = output
    postdict['posterior'] = posterior
    postdict['samples'] = posterior.samples
    postdict['Nsamples'] = len(posterior.samples)
    postdict['parnames'] = posterior.getParamNames().list()
    postdict['logZ'] = output.logZ
    print(f'ln(Evidence) for k={i}: {output.logZ}')
    posteriors.append(postdict)

# Prior distribution for planet models
if args_params.prior == 'u':
    pk = 1./(k+1)
# Probability of there being k planets given the data
pky = np.zeros(k+1)
for i, post in enumerate(posteriors):
    pky[i] = np.exp(post['logZ']) * pk
pky = pky/np.sum(pky)  # Normalization

Tobs = 3060.81  # Total time of observation
factor = 1
# With the 1/p0 + 1/Tobs criterion
# Doing it this way there is no overlap and the right edge of one interval
# coincides with the left edge of the next one
# pn = 1.25
# center_points = []
# while pn < (Tobs*factor)/2:
#     center_points.append(pn)
#     pn = 1./((1/pn)-(2/(factor*Tobs)))

width = 4/Tobs
N = int((1-1e-4)/width)
grid_freq = np.linspace(1e-4, 1, N)
grid_period = 1/grid_freq[::-1]

joint_probability = np.zeros(len(grid_period))
for i in range(len(grid_period)-1):
    suma = 0
    # Count number of samples which have at least one of the periods in the
    # given range
    left_edge = grid_period[i]
    right_edge = grid_period[i+1]
    # print([left_edge, pt, right_edge])
    for n, postdict in enumerate(posteriors):
        # Get Indeces for periods
        idx = []
        for par in postdict['parnames']:
            if 'period' in par:
                idx.append(postdict['parnames'].index(par))

        # Get all periods for this configuration
        periods = np.zeros([postdict['Nsamples'], n])
        truthtable = np.ones([postdict['Nsamples'], n], dtype=bool)
        for j in range(n):
            periods[:, j] = postdict['samples'][:, idx[j]]
            # Set True for values that are within the interval and False otherwise
            truthtable[:, j] = (periods[:, j] >= left_edge) & (
                periods[:, j] < right_edge)
        # Reduce truthtable to one array with a bitwise or function along axis 1
        reduced = np.bitwise_or.reduce(truthtable, axis=1)
        # Count how many samples have at least one of the periods in the interval
        samp_in_interval = np.count_nonzero(reduced)

        suma += (samp_in_interval/postdict['Nsamples']) * pky[n]

    # Add result to joint probability array
    joint_probability[i] = suma

    # Print Progress
    progress = np.round(100*(i/len(grid_period)), 1)
    print(f'Progress: {progress}%', end='\r')

print(sum(joint_probability))

pino0 = joint_probability.nonzero()
grid_period = np.array(grid_period)

# logp = 1-joint_probability
# plt.plot(grid_period, logp, 'k.-')
# plt.yscale('log')
plt.plot(grid_period, joint_probability, 'k.-')
# plt.hist(joint_probability, bins=grid_period, histtype='step')

plt.xscale('log')
plt.xlabel('Period')
plt.ylabel('Probability')
plt.title(f'kmax = {k}')
plt.show()
