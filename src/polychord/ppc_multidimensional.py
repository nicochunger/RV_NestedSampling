#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_multidimensional.py

import PyPolyChord as PPC 
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior
import numpy as np
from scipy.stats import norm, multivariate_normal
from scipy.misc import logsumexp
import matplotlib.pyplot as plt
import time

# Initialize start time to measure run time
start = time.time()

# Number of dimensions and derived parameters
nDims = 5
nDerived = 0

# Min and max values for all parameters (in this simple case)
Tmin = -60
Tmax = 60

def logLikelihood(theta):
    """ log Likelihood of this toy multiimodal example. """

    #mean = [[3,3], [3,12], [12,3], [12,12]]
    nModes = 10
    mean = []
    currDim = -1
    for i in range(nModes):
        mode = np.zeros([nDims])
        if i % 2 == 0:
            currDim += 1
            currDim = currDim % nDims
            
        loop = int(i / (nDims*2)) + 1
        pos = loop * (-1)**i * 10
        mode[currDim] = pos
        mean.append(mode)

    # # Add 0's to the means if nDims > 2 for every new dimension
    # if nDims > 2:
    #     for i in mean:
    #         for _ in range(nDims-2):
    #             i.append(0)

    rvs = []
    for i in range(nModes):
        rvs.append(multivariate_normal(mean[i]).logpdf(theta))

    result = logsumexp(rvs)


    # rv1 = multivariate_normal(mean[0])
    # rv2 = multivariate_normal(mean[1])
    # rv3 = multivariate_normal(mean[2])
    # rv4 = multivariate_normal(mean[3])

    # result1 = np.logaddexp(rv1.logpdf(theta), rv2.logpdf(theta))
    # result2 = np.logaddexp(rv3.logpdf(theta), rv4.logpdf(theta))
    # result = np.logaddexp(result1, result2)
    
    return result, []

def prior(hypercube):
    """ Uniform Prior for [Tmin, Tmax]. """
    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = UniformPrior(Tmin, Tmax)(x)

    return theta

# Define PolyChord settings
settings = PolyChordSettings(nDims, nDerived, )
settings.do_clustering = True
settings.nlive = 25*nDims
settings.file_root = 'multidimensional'
settings.read_resume = False
#settings.num_repeats = nDims * 5

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)

print(output.logZs)
print(output.logZerrs)

end = time.time()
print(f'Total run time was: {end-start}')

try:
    import getdist.plots
    import matplotlib.pyplot as plt
    #plt.close('all')
    posterior = output.posterior
    g = getdist.plots.getSubplotPlotter()
    #print("hola")
    g.triangle_plot(posterior, filled=True)
    plt.show()
except ImportError:
    print("Install matplotlib and getdist for plotting examples")
