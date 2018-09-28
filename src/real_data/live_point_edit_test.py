#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

# Add path of model files
import sys
path = '/home/nunger/tesis/codigo/src/real_data/'
sys.path.append(path)

# Dependencies
from model_rvrd_kepler import lnlike, lnprior, preprocess
import config
import numpy as np
import argparse
from pprint import pprint

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-n', type=int, default=1,
                    help='Number of planets in model.')

args_params = parser.parse_args()

# Get number of planets in model
nplanets = args_params.n  # Number of Planets in the model

# Assign modelpath
modelpath = path + 'configfiles/hd40307_model.py'

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
    31334.41254,  # Offset
    1.70696,  # Jitter
    1.4646,  # Slope
    0.03492,  # Linear Drift
    0.00193,  # Quadratic Drift
    -0.00124,  # Cubic Drift
    14.29913,  # Drift scaling factor
]

planet1 = [
    2.35612,  # K1
    9.6211,  # Period
    0.11619,  # Eccentricity
    3.5423,  # Oemga
    2.21374  # Mean Anomaly
]

planet2 = [
    2.46472,  # K1
    20.41452,  # Period
    0.05073,  # Eccentricity
    1.3456,  # Oemga
    2.63202  # Mean Anomaly
]

planet3 = [
    1.76562,  # K1
    4.3114,  # Period
    0.07617,  # Eccentricity
    2.47893,  # Oemga
    1.02086  # Mean Anomaly
]

planet4 = [
    0.79496,  # K1
    51.66486,  # Period
    0.09336,  # Eccentricity
    2.15823,  # Oemga
    2.37245  # Mean Anomaly
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
new_logLs = [2*new_logL, new_logL]
new_hypercube = inv_prior(ideal_params)

new_live = new_live_point(new_hypercube, ideal_params, new_logLs)

print(new_live)
