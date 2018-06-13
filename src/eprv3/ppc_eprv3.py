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

# PolyChord imports
import PyPolyChord as PPC 
from PyPolyChord.settings import PolyChordSettings

# Initialize start time to measure run time
start = time.time()

# Generate dictionaries
nplanets = 2 # Number of Planets in the model
modelpath = f'configfiles/eprv3rv01_k{nplanets}.py'
rundict, initial_values, datadict, priordict, fixedpardict = config.read_config(path + modelpath)
covdict = preprocess(datadict)[0] # Covariance dictionary
parnames = list(initial_values.keys()) # Parameter names

nDims = 2 + (nplanets * 5) # Number of parameters to fit
assert nDims == len(parnames), 'Error in number of parameters'
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
settings.do_clustering = False
settings.nlive = 25*nDims
settings.file_root = modelpath[12:-3]
settings.read_resume = False
#settings.num_repeats = nDims * 5
settings.feedback = 1
settings.precision_criterion = 0.001

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
result[3] = output.logZ*0.43429 # Total evidence in log_10
result[4] = settings.nlive # Number of live points
result[5] = settings.precision_criterion # Precision crtierion
result[6:] = np.mean(output.posterior.samples[-500:], axis=0) # Average for each parameter from the posterior
result = np.reshape(result, (1,6+nDims))

header = 'run_time logZ logZerr log10Z nlive prec '
for i in range(nDims):
    header += parnames[i]
    if i < nDims-1:
        header += ' ' # Add comma after each parameter, except the last one
try:
    # Append results to file
    f = np.loadtxt(f'results{nplanets}.txt')
    if len(np.shape(f)) == 1:
        f = np.reshape(f, (1, 6+nDims))
    results = np.append(f, result, axis=0)
    np.savetxt(f'results{nplanets}.txt',results, header=header, fmt='%.8e')
except:
    # File does not exist, must create it first
    np.savetxt(f'results{nplanets}.txt', result, header=header, fmt='%.8e')

if nDims < 6:
    try:
        import getdist.plots
        import matplotlib.pyplot as plt
        posterior = output.posterior
        g = getdist.plots.getSubplotPlotter()
        g.triangle_plot(posterior, filled=True)
        plt.show()
    except ImportError:
        print("Install matplotlib and getdist for plotting examples")




