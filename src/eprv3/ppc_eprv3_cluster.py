#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

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
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Local libraries
from model_eprv3_kepler import lnlike, lnprior, preprocess
import config_cluster

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
parser.add_argument('-prec', type=float, default=0.002,
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
    1, 7), "Incorrect datafile. Has to be an integer from 1 to 6."

# Generate dictionaries
nplanets = args_params.n  # Number of Planets in the model

filepath = os.path.dirname(os.path.abspath(__file__))
modelpath = os.path.join(filepath, 'eprv3rv01_model.py')

parnames, datadict, priordict, fixedpardict = config_cluster.read_config(
    modelpath, nplanets, args_params.dfile, args_params.narrow)
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


dirname = '/scratch/nunger/eprv3'
timecode = time.strftime("%m%d_%H%M")
folder_path = '000{}_{}a_'.format(datafile, nplanets) + timecode
if args_params.narrow:
    folder_path = '000{}_{}b_'.format(datafile, nplanets) + timecode

# Define PolyChord settings
settings = PolyChordSettings(nDims, nDerived, )
settings.do_clustering = not args_params.noclust
settings.nlive = args_params.nlive
settings.num_repeats = nDims * args_params.nrep

settings.base_dir = dirname+'/chains/'
if args_params.save:
    # Save all the files from this run
    settings.base_dir = dirname+'/saved_runs/'+folder_path
print(settings.base_dir)
settings.file_root = 'eprv3'
settings.read_resume = False
settings.write_resume = False
settings.precision_criterion = args_params.prec

# Initialize start time to measure run time
start = time.time()

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)

end = time.time()  # End time
Dt = end - start
if rank == 0:
    print('\nTotal run time was: {}'.format(
        datetime.timedelta(seconds=int(Dt))))
    # Log10 of the evidence
    print('\nlog10(Z) = {} \n'.format(output.logZ*0.43429))

    # Parameter names
    latexnames = [r'\sigma_J', r'C']
    for j in range(nplanets):
        latexnames.extend(
            [r'K_{}'.format(j), r'P_{}'.format(j), r'e_{}'.format(j), r'\omega_{}'.format(j), r'M_{}'.format(j)])
    paramnames = [(x, latexnames[i]) for i, x in enumerate(parnames)]
    output.make_paramnames_files(paramnames)

    # Save output as a pickle file
    if args_params.save:
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
    # medians = np.median(output.posterior.samples, axis=0)
    # for i in range(nDims):
    #     results[parnames[i]] = medians[i]

    # Convert to pandas DataFrame
    results = pd.DataFrame(results, index=[0])
    # Order the parameters
    order = ['run_time', 'logZ', 'logZerr', 'log10Z', 'nlive', 'prec']
    # for par in parnames:
    #     order.append(par)
    results = results[order]

    # Name of data file
    filename = dirname + \
        '/results2/000{}/results000{}_{}a.txt'.format(
            datafile, datafile, nplanets)
    if args_params.narrow:
        filename = dirname + \
            '/results2/000{}/results000{}_{}b.txt'.format(
                datafile, datafile, nplanets)

    try:
        # Append results to file
        f = pd.read_csv(filename, sep='\t')
        f = f.append(results)
        f = f[order]
        f.to_csv(filename, sep='\t', index=False, float_format='%8.5f')
        print('Guardado a:' + filename)
    except:
        # File does not exist, must create it first
        results.to_csv(filename, sep='\t', index=False, float_format='%8.5f')
        print('Guardado a:' + filename)
