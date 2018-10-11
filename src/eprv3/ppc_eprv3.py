#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

# Add path of config files
import sys
path = '/home/nunger/tesis/codigo/src/eprv3/'
sys.path.append(path)

# Dependencies
# Standard libraries
import time
import datetime
import argparse
import os

# Related Third party imports
import numpy as np
import pandas as pd
import pickle

# Local libraries
from model_eprv3_kepler import lnlike, lnprior, preprocess
import config

# PolyChord imports
import PyPolyChord as PPC
from PyPolyChord.settings import PolyChordSettings

# # Remove stack size limit
# import resource
# resource.setrlimit(resource.RLIMIT_STACK,
#                    (resource.RLIM_INFINITY, resource.RLIM_INFINITY))

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-n', type=int, default=1,
                    help='Number of planets in model.')
parser.add_argument('-nlive', type=int, default=25,
                    help='Number of live points in nested sampling. As: nlive*ndim')
parser.add_argument('-nrep', type=int, default=3,
                    help='Number of repeats in slice sampling. As: nrep*ndim')
parser.add_argument('-prec', type=float, default=0.01,
                    help='Precision criterion for termination.')
parser.add_argument('-dfile', type=int, default=1,
                    help='Which dataset to analyze.')

parser.add_argument('--noclust', action='store_true',
                    help='If clustering will be deactivated.')
parser.add_argument('--narrow', action='store_true',
                    help='If the narrowpriors will be used.')
parser.add_argument('--save', action='store_true',
                    help='If the run should be saved.')
args_params = parser.parse_args()

datafile = args_params.dfile  # Data set to analyze
assert datafile in range(
    1, 7), "Incorrect datafile. Has to an integer from 1 to 6."

# Generate dictionaries
nplanets = args_params.n  # Number of Planets in the model
modelpath = 'eprv3rv01_model.py'
parnames, datadict, priordict, fixedpardict = config.read_config(
    path + modelpath, nplanets, args_params.dfile, args_params.narrow)
covdict = preprocess(datadict)[0]  # Covariance dictionary

nDims = 2 + (nplanets * 5)  # Number of parameters to fit
assert nDims == len(
    parnames), "Number of parameters and dimensions don't match"
nDerived = 0


def logLikelihood(theta):
    """ Log Likelihood for the planet model. """

    result = lnlike(theta, parnames, fixedpardict, datadict, covdict)

    return result, []


def prior(hypercube):
    """ Priors for each parameter. """
    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = priordict[parnames[i]].ppf(x)

    return theta


dirname = os.path.dirname(os.path.abspath(__file__))
timecode = time.strftime("%m%d_%H%M")
folder_path = f'000{datafile}_{nplanets}a_' + timecode
if args_params.narrow:
    folder_path = f'000{datafile}_{nplanets}b_' + timecode

# Define PolyChord settings
settings = PolyChordSettings(nDims, nDerived, )
settings.do_clustering = not args_params.noclust
settings.nlive = nDims * args_params.nlive
settings.num_repeats = nDims * args_params.nrep

settings.base_dir = 'chains/'
if args_params.save:
    # Save all the files from this run
    settings.base_dir = dirname+'/saved_runs/'+folder_path
print(settings.base_dir)
settings.file_root = 'eprv3'
settings.read_resume = False
settings.write_resume = False
settings.precision_criterion = args_params.prec

# Save Parameter names list
try:
    pickle.dump(parnames, open(settings.base_dir+'/parnames.p', 'wb'))
except:
    os.makedirs(settings.base_dir)
    pickle.dump(parnames, open(settings.base_dir+'/parnames.p', 'wb'))

# Initialize start time to measure run time
start = time.time()

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)

# Parameter names
latexnames = [r'\sigma_J', r'C']
for j in range(nplanets):
    latexnames.extend(
        [fr'K_{j}', fr'P_{j}', fr'e_{j}', fr'\omega_{j}', fr'M_{j}'])
paramnames = [(x, latexnames[i]) for i, x in enumerate(parnames)]
output.make_paramnames_files(paramnames)

end = time.time()  # End time
Dt = end - start
print(f'\nTotal run time was: {datetime.timedelta(seconds=int(Dt))}')
print(f'\nlog10(Z) = {output.logZ*0.43429} \n')  # Log10 of the evidence

# Save output as a pickle file
if args_params.save:
    pickle_file = settings.base_dir + '/output.p'
    pickle.dump(output, open(pickle_file, "wb"))

# Save evidence and other relevant data
results = {}
results['run_time'] = Dt
results['logZ'] = output.logZ
results['logZerr'] = output.logZerr
results['log10Z'] = output.logZ * np.log10(np.e)  # Total evidence in log_10
results['nlive'] = settings.nlive  # Number of live points
results['prec'] = settings.precision_criterion  # Precision crtierion
medians = np.median(output.posterior.samples, axis=0)
for i in range(nDims):
    results[parnames[i]] = medians[i]

# Convert to pandas DataFrame
results = pd.DataFrame(results, index=[0])
# Order the parameters
order = ['run_time', 'logZ', 'logZerr', 'log10Z', 'nlive', 'prec']
for par in parnames:
    order.append(par)
results = results[order]

# Name of data file
filename = f'results/000{datafile}/results000{datafile}_{nplanets}a_pd.txt'
if args_params.narrow:
    filename = f'results/000{datafile}/results000{datafile}_{nplanets}b_pd.txt'

try:
    # Append results to file
    f = pd.read_csv(filename, sep='\t')
    f = f.append(results)
    f.to_csv(filename, sep='\t', index=False, float_format='%8.5f')
except:
    # File does not exist, must create it first
    results.to_csv(filename, sep='\t', index=False, float_format='%8.5f')


# # Plotting
# if nDims < 8:
#     try:
#         import getdist.plots
#         import matplotlib.pyplot as plt
#         posterior = output.posterior
#         g = getdist.plots.getSubplotPlotter()
#         g.triangle_plot(posterior, filled=True)
#         plt.show()
#     except ImportError:
#         print("Install matplotlib and getdist for plotting examples")
