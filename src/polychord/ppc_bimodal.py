#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_bimodal.py

import PyPolyChord as PPC
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import time
import datetime

nDims = 1
nDerived = 0

Tmin = 0
Tmax = 10


def PLUS(x, y):
    """ Logarithmic sum. """
    if x > y:
        return x + np.log(1+np.exp(y-x))
    else:
        return y + np.log(1+np.exp(x-y))


def logLikelihood(theta):
    """ log Likelihood of this toy bimodal example. """
    loc1 = 3
    sigma1 = 0.4
    loc2 = 7
    sigma2 = 0.4
    #result = np.log(norm.pdf(theta, loc1, sigma1) + norm.pdf(theta, loc2, sigma2))
    result = PLUS(norm.logpdf(theta, loc1, sigma1),
                  norm.logpdf(theta, loc2, sigma2))

    return float(result), []


def prior(hypercube):
    """ Uniform Prior for [Tmin, Tmax]. """
    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = UniformPrior(Tmin, Tmax)(x)

    return theta


# Define PolyChord settings
settings = PolyChordSettings(nDims, nDerived)
settings.do_clustering = True
settings.nlive = 200
settings.file_root = 'bimodal'
settings.read_resume = False

start = time.time()

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)

# End time
end = time.time()
Dt = end - start
print(f'\nTotal run time was: {datetime.timedelta(seconds=int(Dt))}')

print(f'\nZ = {np.exp(output.logZ)}')

result = np.exp(output.logZ)
filename = 'results.txt'

try:
    # Append results to file
    f = np.loadtxt(filename)
    if len(np.shape(f)) == 0:
        f = np.atleast_1d(f)
    results = np.append(f, np.atleast_1d(result), axis=0)
    np.savetxt(filename, results, header='evidence', fmt='%.8e')
except:
    # File does not exist, must create it first
    np.savetxt(filename, np.atleast_1d(result), header='evidence', fmt='%.8e')

# try:
#     import getdist.plots
#     import matplotlib.pyplot as plt
#     #plt.close('all')
#     posterior = output.posterior
#     g = getdist.plots.getSubplotPlotter()
#     #print("hola")
#     g.triangle_plot(posterior, filled=True)
#     plt.show()
# except ImportError:
#     print("Install matplotlib and getdist for plotting examples")
