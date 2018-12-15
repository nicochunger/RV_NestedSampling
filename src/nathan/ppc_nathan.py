#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

# Dependencies
from model_nathan_kepler import lnlike, lnprior, preprocess
import config
import numpy as np
import pandas as pd
import time
import datetime
import argparse
import pickle
import subprocess
import os

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# PolyChord imports
import PyPolyChord as PPC
from PyPolyChord.settings import PolyChordSettings

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-model', type=int, default=1,
                    help='Model to be used')
parser.add_argument('-dfile', type=int, default=1,
                    help='Datafile to be used')
parser.add_argument('-nlive', type=int, default=25,
                    help='Number of live points in nested sampling. As: nlive*ndim')
parser.add_argument('-nrep', type=int, default=3,
                    help='Number of repeats in slice sampling. As: nrep*ndim')
parser.add_argument('-prec', type=float, default=0.01,
                    help='Precision criterion for termination.')
# parser.add_argument('-narrow', type=float, default=0,
#                     help='Wether to use the narrow priors for the periods.')

parser.add_argument('--noclust', action='store_false',
                    help='If clustering should be turned off.')

parser.add_argument('--resume', action='store_true',
                    help='If the previous run should be resumed.')

parser.add_argument('--save', action='store_true',
                    help='If the posterior samples should be saved.')

parser.add_argument('--cluster', action='store_true',
                    help='If the run is going to be run in the cluster.')

args_params = parser.parse_args()

# Initialize start time to measure run time
start = time.time()

# Assign modelpath
model = args_params.model
filepath = os.path.dirname(__file__)
modelpath = os.path.join(filepath, 'models/nathan_newmodel{}.py'.format(model))

# Generate dictionaries
datafile = args_params.dfile
data_path = os.path.join(filepath, 'data/data_2gen/')
if args_params.cluster:
    data_path = '/home/spectro/nunger/codigo/src/nathan/data/data_2gen/'
data_files = subprocess.check_output(
    'ls {}'.format(data_path), shell=True).decode('utf-8').split('\n')
data_files.remove('')  # Remove last empty item
assert datafile in range(
    1, 201), "Incorrect datafile has to be one of {}".format(range(1, 101))
parnames, datadict, priordict, fixedpardict = config.read_config(
    modelpath, data_path, datafile)

if rank == 0:
    print('Datafile: {}'.format(data_files[datafile-1]))
    print('\n Parameter names and order:')
    print(parnames)

# Number of dimensions is the number of parameters in the model
nDims = len(parnames)
nDerived = 0


def logLikelihood(theta):
    """ Log Likelihood for the planet model. """

    result = lnlike(theta, parnames, fixedpardict, datadict)

    return result, []


def prior(hypercube):
    """ Priors for each parameter. """
    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = priordict[parnames[i]].ppf(x)

    return theta


# Filepath for the polychord data and storing the samples
dirname = filepath
if args_params.cluster:
    dirname = '/scratch/nunger/nathan/'
timecode = time.strftime("%m%d_%H%M")
folder_path = 'nathan_model{}_'.format(model) + timecode
if not args_params.save:
    # Save the samples in a dump folder if the data shouldn't be saved
    # This overwrites any data saved before of the same model
    folder_path = 'dump'

# Define PolyChord settings
settings = PolyChordSettings(nDims, nDerived, )
settings.do_clustering = args_params.noclust
settings.nlive = nDims * args_params.nlive
settings.base_dir = os.path.join(dirname, folder_path)
settings.file_root = 'nathan_model{}'.format(model)
settings.num_repeats = nDims * args_params.nrep
settings.precision_criterion = args_params.prec
settings.read_resume = False

# Change settings if resume is true
# if args_params.resume:
#     settings.read_resume = args_params.resume
#     settings.base_dir = dirname+prev_run

# Save Parameter names list
try:
    pickle.dump(parnames, open(settings.base_dir+'/parnames.p', 'wb'))
except:
    os.makedirs(settings.base_dir)
    pickle.dump(parnames, open(settings.base_dir+'/parnames.p', 'wb'))

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)

# Parameter names
# latexnames = [r'\sigma_J', r'C']
# for j in range(nplanets):
#     latexnames.extend(
#         [fr'K_{j}', fr'P_{j}', fr'e_{j}', fr'\omega_{j}', fr'M_{j}'])
# paramnames = [(x, latexnames[i]) for i, x in enumerate(parnames)]
paramnames = [(x, x) for x in parnames]

output.make_paramnames_files(paramnames)

end = time.time()  # End time
Dt = end - start
if rank == 0:
    print('\nTotal run time was: {}'.format(
        datetime.timedelta(seconds=int(Dt))))
    # Log10 of the evidence
    print('\nlog10(Z) = {} \n'.format(output.logZ*0.43429))

    # Save output data as a pickle file
    pickle_file = settings.base_dir + '/output.p'
    pickle.dump(output, open(pickle_file, "wb"))

    # Save evidence and other relevant data
    results = {}
    results['run_time'] = Dt
    results['logZ'] = output.logZ
    results['logZerr'] = output.logZerr
    results['log10Z'] = output.logZ * \
        np.log10(np.e)  # Total evidence in log_10
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

    print('\nParameters:')
    print(results)

    # Name of data file
    filename = os.path.join(dirname, 'results_2gen/results_{}_model{}.txt'.format(
        data_files[datafile-1][5:-4], model))

    try:
        # Append results to file
        f = pd.read_csv(filename, sep='\t')
        f = f.append(results)
        f = f[order]
        f.to_csv(filename, sep='\t', index=False, float_format='%8.5f')
    except:
        # File does not exist, must create it first
        results.to_csv(filename, sep='\t', index=False, float_format='%8.5f')
