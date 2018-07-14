#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

# Add path of config files
import sys
path = '/home/nunger/tesis/codigo/run/'
sys.path.append(path)

# Dependencies
from model_eprv3_kepler import lnlike, lnprior, preprocess
import config
import numpy as np
import time
import datetime
import argparse

# PolyChord imports
import PyPolyChord as PPC 
from PyPolyChord.settings import PolyChordSettings

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-n', type=int, default=1)
parser.add_argument('-nlive', type=int, default=25)
parser.add_argument('-nrep', type=int, default=3)
parser.add_argument('-prec', type=float, default=0.01)
parser.add_argument('-dfile', type=int, default=1)
parser.add_argument('--clust', action='store_true')
parser.add_argument('--narrow', action='store_true')
args_params = parser.parse_args()

datafile = args_params.dfile # Data set to analyze
assert datafile in range(1,7), "Incorrect datafile. Has to be between 1 and 6."

# Initialize start time to measure run time
start = time.time()

# Generate dictionaries
nplanets = args_params.n # Number of Planets in the model
modelpath = f'configfiles/eprv3rv01_k{nplanets}.py'
if args_params.narrow:
    splitpath = modelpath.split('.')
    modelpath = splitpath[0] + '_narrowprior.' + splitpath[1]

#---- Edit configfile to set the correct datafile ------
with open(path+modelpath) as f:
    lines = f.readlines()

for i in range(len(lines)):
    if 'datafile' in lines[i]:
        # Edit that line to insert datafile
        lines[i] = lines[i][:49] + str(datafile) + lines[i][50:]

with open(path+modelpath, "w") as f:
    # Write edited file
    f.write(''.join(lines))
#--------------------------------------------------------

rundict, initial_values, datadict, priordict, fixedpardict = config.read_config(path + modelpath)
covdict = preprocess(datadict)[0] # Covariance dictionary
parnames = list(initial_values.keys()) # Parameter names

nDims = 2 + (nplanets * 5) # Number of parameters to fit
assert nDims == len(parnames), "Number of parameters and dimensions don't match"
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
settings.do_clustering = args_params.clust
settings.nlive = nDims * args_params.nlive
settings.file_root = modelpath[12:-3]
settings.read_resume = False
settings.num_repeats = nDims * args_params.nrep
#settings.feedback = 1
settings.precision_criterion = args_params.prec

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)
# Parameter names
latexnames = [r'\sigma_J', r'C']
for j in range(nplanets):
    latexnames.extend([fr'K_{j}', fr'P_{j}', fr'e_{j}', fr'\omega_{j}', fr'M_{j}'])
paramnames = [(x, latexnames[i]) for i, x in enumerate(parnames)]

output.make_paramnames_files(paramnames)

end = time.time() # End time
Dt = end - start
print(f'\nTotal run time was: {datetime.timedelta(seconds=int(Dt))}')
print(f'\nlog10(Z) = {output.logZ*0.43429} \n') # Log10 of the evidence

# Save evidence and other relevant data
result = np.zeros(6+nDims)
result[0] = Dt # Run time
result[1] = output.logZ # Total evidence in log_e
result[2] = output.logZerr # Error for the evidence
result[3] = output.logZ * np.log10(np.e) # Total evidence in log_10
result[4] = settings.nlive # Number of live points
result[5] = settings.precision_criterion # Precision crtierion
result[6:] = np.mean(output.posterior.samples[-500:], axis=0) # Average for each parameter from the posterior
result = np.reshape(result, (1,6+nDims))

header = 'run_time logZ logZerr log10Z nlive prec '
for i in range(nDims):
    header += parnames[i]
    if i < nDims-1:
        header += ' ' # Add space after each parameter, except the last one

dataset = datadict['eprv']['datafile'][-8:-4]
# Name of data file
if 'narrowprior' not in modelpath:
    filename = f'results/results{dataset}_{nplanets}a.txt'
else:
    filename = f'results/results{dataset}_{nplanets}b.txt'

try:
    # Append results to file
    f = np.loadtxt(filename)
    if len(np.shape(f)) == 1:
        f = np.reshape(f, (1, 6+nDims))
    results = np.append(f, result, axis=0)
    np.savetxt(filename,results, header=header, fmt='%.8e')
except:
    # File does not exist, must create it first
    np.savetxt(filename, result, header=header, fmt='%.8e')

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




