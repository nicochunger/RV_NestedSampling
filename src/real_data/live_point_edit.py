#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

# Add path of model files
# import sys
# path = '/codigo/src/real_data/'
# sys.path.append(path)

# Dependencies
from model_rvrd_kepler import lnlike, lnprior, preprocess
import config
import numpy as np
import argparse
from pprint import pprint
import subprocess

# Read last run to edit
chains_path = '/scratch/nunger/polychord_chains/'
# chains_path = '/media/nunger/Windows/Nico/Facu/Tesis/polychord_chains/'
runs = subprocess.check_output(
    'ls -d '+chains_path+'*/', shell=True).decode('utf-8').replace(chains_path, '').split('/\n')
dates = []
for run in runs:
    dates.append(run[-9:])
# dates.remove('dump')
# dates.remove('Cluster')
dates.sort()
for run in runs:
    # Search for the most recent run
    if dates[-1] in run:
        prev_run = run
# Extract the number of planets analyzed in previous run
nplanets = int(prev_run.split('_')[1][0])

# Assign modelpath
modelpath = 'configfiles/hd40307_model_vizier_cluster.py'

# Generate dictionaries
parnames, datadict, priordict, fixedpardict = config.read_config(
    modelpath, nplanets, False)

print('\n Parameter names and order:')
print(parnames)
print('\n')

covdict = preprocess(datadict)[0]  # Covariance dictionary

# Number of dimensions is the number of parameters in the model
nDims = len(parnames)
nDerived = 0


def logLikelihood(theta):
    """ Log Likelihood for the planet model. """

    result = lnlike(theta, parnames, fixedpardict, datadict, covdict)

    return result, []


def prior(hypercube):
    """ Priors for each parameter. """
    theta = np.zeros(nDims)
    for i, x in enumerate(hypercube):
        theta[i] = priordict[parnames[i]].ppf(x)

    return theta


def inv_prior(parameters):
    """ Priors for each parameter. """
    hypercube = np.zeros(nDims)
    for i, x in enumerate(parameters):
        hypercube[i] = priordict[parnames[i]].cdf(x)

    return hypercube


def live_point_parser(live):
    """ Recieves a string of one live point from the resume file
    and return 3 arrays. The array for the hypercube, for the parameters,
    and likelihood/s. """
    live = live[2:]  # Remove leading spaces
    # Split string
    points = live.split(' ')
    # Remove empty strings
    while '' in points:
        points.remove('')
    # Allocate arrays
    hypercube = np.array([float(i) for i in points[:nDims]])
    parameters = np.array([float(i) for i in points[nDims:nDims*2]])
    likelihoods = np.array([float(i) for i in points[-2:]])
    likelihoods[0] = likelihoods[1]*2

    return hypercube, parameters, likelihoods


def new_live_point(hypercube, parameters, likelihoods):
    """ 
    Takes the new hypercube and likelihoods and forms them
    together to a string.
    """
    string = ''
    # Append hypercube values
    for h in hypercube:
        # If number is negative there is only one space before
        if h < 0:
            string += ' ' + format_string(h)
        # If positive there are two spaces
        else:
            string += '  ' + format_string(h)

    # Append parameter values
    for p in parameters:
        if p < 0:
            string += ' ' + format_string(p)
        else:
            string += '  ' + format_string(p)

    # Append Likelihoods
    string += ' ' + format_string(likelihoods[0])
    string += ' ' + format_string(likelihoods[1])

    return string


def format_string(n):
    p = 0
    while n <= -1 or n >= 1:
        n = n / 10.0
        p += 1
    while n < 1e-1 and n > -1e-1:
        n = n * 10.0
        p -= 1
    return "{:.15f}E{:+04d}".format(n, p)


harps = [
    31334.412,  # Offset
    1.70696,    # Jitter
    5.24698,    # Slope
    0.03765,    # Linear Drift
    0.00484,    # Quadratic Drift
    -0.00161,   # Cubic Drift
    22.2,       # Drift scaling factor
]

planet1 = [
    2.35612,    # K1
    9.6211,     # Period
    0.11619,    # Eccentricity
    3.5423,     # Oemga
    2.21374     # Mean Anomaly
]

planet2 = [
    2.46472,    # K1
    20.41452,   # Period
    0.05073,    # Eccentricity
    1.3456,     # Oemga
    2.63202     # Mean Anomaly
]

planet3 = [
    1.76562,    # K1
    4.3114,     # Period
    0.07617,    # Eccentricity
    2.47893,    # Oemga
    1.02086     # Mean Anomaly
]

planet4 = [
    0.79496,    # K1
    51.66486,   # Period
    0.09336,    # Eccentricity
    2.15823,    # Oemga
    2.37245     # Mean Anomaly
]

planets = [planet1, planet2, planet3, planet4]

# Construct ideal params array
ideal_params = []
for par in harps:
    ideal_params.append(par)
for n in range(nplanets):
    for i in planets[n]:
        ideal_params.append(i)

# Convert to numpy array
ideal_params = np.array(ideal_params)

new_logL = logLikelihood(ideal_params)[0]
new_logLs = [-1e30, new_logL]
new_hypercube = inv_prior(ideal_params)

# Generate string of new live point
new_live = new_live_point(new_hypercube, ideal_params, new_logLs)


dirname = chains_path + prev_run + '/hd40307_k{}.resume'.format(nplanets)

# Read and edit resume file
with open(dirname, 'r') as resume_file:
    resume = resume_file.readlines()
    for i in range(len(resume)):
        if 'Minimum loglikelihood positions' in resume[i]:
            min_lhood = int(resume[i+1])
            continue
        if '=== live points ===' in resume[i]:
            if min_lhood == 0:
                print('It was the first point.')
                resume[i+3] = new_live + '\n'
            else:
                print('\nPrevioud point:')
                print(resume[i+2])
                resume[i+2] = new_live + '\n'
                print('\nNew point:')
                print(resume[i+2])
            break

# Write edits
with open(dirname, 'w') as file:
    file.writelines(resume)
