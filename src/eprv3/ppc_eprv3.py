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

# PolyChord imports
import PyPolyChord as PPC 
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior

# Initialize start time to measure run time
start = time.time()

# Generate dictionaries
nplanets = 1 # Number of Planets in the model
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
settings.do_clustering = True
#settings.nlive = 25*nDims
settings.file_root = modelpath[12:-3]
settings.read_resume = True
#settings.num_repeats = nDims * 5

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)
#latexnames = [r'\sigma_J', r'C']
# latexplantes = [r'K_1', r'P_1', r'e_1', r'\omega_1', r'M_1']
# paramnames = [(x, latexnames[i]) for i, x in enumerate(parnames)]
latexnames = [r'\sigma_J', r'C']
for j in range(nplanets):
    latexnames.extend([fr'K_{j}', fr'P_{j}', fr'e_{j}', fr'\omega_{j}', fr'M_{j}'])
paramnames = [(x, latexnames[i]) for i, x in enumerate(parnames)]

output.make_paramnames_files(paramnames)

#print(output.logZs)
#print(output.logZerrs)

end = time.time()
print(f'Total run time was: {end-start}')
print(f'\n log10(Z) = {output.logZ*0.43429} \n') # Log10 of the evidence

try:
    import getdist.plots
    import matplotlib.pyplot as plt
    posterior = output.posterior
    g = getdist.plots.getSubplotPlotter()
    g.triangle_plot(posterior, filled=True)
    plt.show()
except ImportError:
    print("Install matplotlib and getdist for plotting examples")




