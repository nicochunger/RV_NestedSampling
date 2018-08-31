#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

# Dependencies
from model_rvrd_kepler import lnlike, lnprior, preprocess
import config
import numpy as np
import time
import datetime
import argparse
import pickle
import os
import subprocess

# PolyChord imports
import PyPolyChord as PPC
from PyPolyChord.settings import PolyChordSettings

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-n', type=int, default=1,
                    help='Number of planets in model.')
parser.add_argument('-nlive', type=int, default=25,
                    help='Number of live points in nested sampling. As: nlive*ndim')
parser.add_argument('-nrep', type=int, default=5,
                    help='Number of repeats in slice sampling. As: nrep*ndim')
parser.add_argument('-prec', type=float, default=0.01,
                    help='Precision criterion for termination.')

# TODO Change this parameter for the new real data
# parser.add_argument('-inst', type=str, default=1,
#                     help='Which dataset to analyze.')

parser.add_argument('--clust', action='store_true',
                    help='If clustering will be activated.')

parser.add_argument('--resume', action='store_true',
                    help='If the previous run should be resumed.')

parser.add_argument('--save', action='store_true',
                    help='If the posterior samples should be saved.')

args_params = parser.parse_args()

# Initialize start time to measure run time
start = time.time()

# Get number of planets in model
nplanets = args_params.n  # Number of Planets in the model

# Change number of plantes if resume is true
dirname = '/media/nunger/Windows/Nico/Facu/Tesis/polychord_chains/'
prev_run = subprocess.check_output(
    'ls -d '+dirname+'*/', shell=True).decode('utf-8').replace(dirname, '').split('/\n')[-2]
if args_params.resume:
    # Extract the number of planets analyzed in previous run
    nplanets = int(prev_run.split('_')[1][0])

# Assign modelpath
modelpath = f'configfiles/hd40307_k{nplanets}.py'


# Generate dictionaries
parnames, datadict, priordict, fixedpardict = config.read_config(modelpath)
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


# Filepath for the polychord data and storing the samples
timecode = time.strftime("%m%d_%H%M")
folder_path = f'hd40307_{nplanets}a_' + timecode
if not args_params.save:
    # Save the samples in a dump folder if the data shouldn't be saved
    # This overwrites any data saved before of the same model
    folder_path = 'dump'

# Define PolyChord settings
settings = PolyChordSettings(nDims, nDerived, )
settings.do_clustering = args_params.clust
settings.nlive = nDims * args_params.nlive
settings.base_dir = dirname+folder_path
settings.file_root = modelpath[12:-3]
settings.num_repeats = nDims * args_params.nrep
settings.precision_criterion = args_params.prec
settings.read_resume = False

# Change settings if resume is true
if args_params.resume:
    settings.read_resume = args_params.resume
    settings.base_dir = dirname+prev_run

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

# Save output data as a pickle file
if not args_params.save:
    pickle_file = settings.base_dir + '/output.p'
    pickle.dump(output, open(pickle_file, "wb"))

# Save evidence and other relevant data
result = np.zeros(6)
result[0] = Dt  # Run time
result[1] = output.logZ  # Total evidence in log_e
result[2] = output.logZerr  # Error for the evidence
result[3] = output.logZ * np.log10(np.e)  # Total evidence in log_10
result[4] = settings.nlive  # Number of live points
result[5] = settings.precision_criterion  # Precision crtierion
result = np.reshape(result, (1, 6))

header = 'run_time logZ logZerr log10Z nlive prec'
filename = f'results_{nplanets}a.txt'

try:
    # Append results to file
    f = np.loadtxt(dirname+filename)
    if len(np.shape(f)) == 1:
        f = np.reshape(f, (1, 6))
    results = np.append(f, result, axis=0)
    np.savetxt(dirname+filename, results, header=header, fmt='%.8e')
except:
    # File does not exist, must create it first
    np.savetxt(dirname+filename, result, header=header, fmt='%.8e')