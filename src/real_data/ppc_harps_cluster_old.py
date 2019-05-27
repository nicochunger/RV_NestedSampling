#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

# Python Packages
import time
import datetime
import argparse
import subprocess
import os
import random

# Dependencies
import numpy as np
import pandas as pd
import pickle
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# print("I'm rank {} of {}".format(rank, size))

# Own modules
from model_rvrd_kepler import lnlike, preprocess
import config

# PolyChord imports
import PyPolyChord as PPC
from PyPolyChord.settings import PolyChordSettings

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
parser.add_argument('-narrow', type=float, default=0,
                    help='Wether to use the narrow priors for the periods.')

parser.add_argument('--noclust', action='store_false',
                    help='If clustering should be turned off.')

parser.add_argument('--resume', action='store_true',
                    help='If the previous run should be resumed.')

parser.add_argument('--save', action='store_true',
                    help='If the posterior samples should be saved.')

parser.add_argument('--cluster', action='store_true',
                    help='If it is going to be run on the cluster.')

parser.add_argument('-nano', type=int, default=1,
                    help='Nano seconds used as seed for current job.')

args_params = parser.parse_args()

# Initialize start time to measure run time
start = time.time()

# Get number of planets in model
nplanets = args_params.n  # Number of Planets in the model

# Change number of plantes if resume is true
dirname = '/media/nunger/Windows/Nico/Facu/Tesis/polychord_chains/definitive'
if args_params.cluster:
    dirname = '/scratch/nunger/polychord_chains/definitive'
if args_params.resume:
    runs = subprocess.check_output(
        'ls -d '+dirname+'*/', shell=True).decode('utf-8').replace(dirname, '').split('/\n')
    dates = []
    for run in runs:
        dates.append(run[-9:])
    dates.sort()
    for run in runs:
        if dates[-1] in run:
            prev_run = run
    print('Corrida a resumir: {}'.format(prev_run))
    # Extract the number of planets analyzed in previous run
    nplanets = int(prev_run.split('_')[1][0])

# Assign modelpath
filepath = os.path.dirname(__file__)
modelpath = os.path.join(
    filepath, 'configfiles/hd40307_model_vizier_cluster.py')

# Generate dictionaries
parnames, datadict, priordict, fixedpardict = config.read_config(
    modelpath, nplanets, args_params.narrow)


# Filepath for the polychord data and storing the samples
timecode = time.strftime("%m%d")
folder_path = 'hd40307_{}a_'.format(nplanets) + timecode
if args_params.narrow:
    folder_path = 'hd40307_{}b_'.format(nplanets) + timecode
if not args_params.save:
    # Save the samples in a dump folder if the data shouldn't be saved
    # This overwrites any data saved before of the same model
    folder_path = 'dump'

random.seed(args_params.nano)
rnd = "%04d" % random.randint(0, 9999)
if folder_path != "dump":
    # Use time in nano seconds as seed to create a unique folder code name
    folder_path += '_' + rnd
    base_dir = os.path.join(dirname, folder_path)
else:
    base_dir = os.path.join(dirname, folder_path)

if rank == 0:
    print('\n Parameter names and order:')
    print(parnames)
    print(base_dir)
    # Save Parameter names list
    try:
        pickle.dump(parnames, open(base_dir+'/parnames.p', 'wb'))
    except:
        os.makedirs(base_dir)
        pickle.dump(parnames, open(base_dir+'/parnames.p', 'wb'))

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
    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = priordict[parnames[i]].ppf(x)

    return theta


# Define PolyChord settings
settings = PolyChordSettings(nDims, nDerived, )
settings.do_clustering = args_params.noclust
settings.nlive = nDims * args_params.nlive
settings.base_dir = base_dir
settings.file_root = 'hd40307_k{}'.format(nplanets)  # modelpath[12:-3]
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

    ########################################

    # Save evidence and other relevant data
    results = {}
    results['id'] = timecode + '_' + rnd
    results['run_time'] = Dt
    results['logZ'] = output.logZ
    results['logZerr'] = output.logZerr
    # results['log10Z'] = output.logZ * \
    #     np.log10(np.e)  # Total evidence in log_10
    results['ndims'] = int(nDims)
    results['nlive'] = int(settings.nlive)  # Number of live points
    results['nreap'] = int(settings.num_repeats)
    results['prec'] = settings.precision_criterion  # Precision crtierion
    # medians = np.median(output.posterior.samples, axis=0)
    # for i in range(nDims):
    #     results[parnames[i]] = medians[i]

    # Convert to pandas DataFrame
    results = pd.DataFrame(results, index=[0])
    # Order the parameters
    order = ['id', 'run_time', 'logZ', 'logZerr',
             'ndims', 'nlive', 'nreap', 'prec']
    # for par in parnames:
    #     order.append(par)
    results = results[order]

    print('\nParameters:')
    print(results)

    # Name of data file
    filename = dirname + '/results_{}a.txt'.format(nplanets)
    if args_params.narrow:
        filename = dirname + '/results_{}b.txt'.format(nplanets)

    try:
        # Append results to file
        f = pd.read_csv(filename, sep='\t')
        f = f.append(results)
        f = f[order]
        f.to_csv(filename, sep='\t', index=False, float_format='%.3f')
    except:
        # File does not exist, must create it first
        results.to_csv(filename, sep='\t', index=False, float_format='%.3f')
