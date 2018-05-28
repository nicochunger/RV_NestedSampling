#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_multidimensional.py

import PyPolyChord as PPC 
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior
import numpy as np
from scipy.stats import norm, multivariate_normal
import matplotlib.pyplot as plt

nDims = 2
nDerived = 0

Tmin = 0
Tmax = 30

def PLUS(x, y):
    """ Logarithmic sum. """
    if x > y:
        return x + np.log(1+np.exp(y-x))
    else:
        return y + np.log(1+np.exp(x-y))

def logLikelihood(theta):
    """ log Likelihood of this toy bimodal example. """
    # loc1 = 3
    # mean1 = [loc1] * nDims
    # mean1 = [3, 3]
    # rv1 = multivariate_normal(mean1)
    # #result = rv1.logpdf(theta)

    # # Multimodal
    # loc2 = 12
    # mean2 = [loc2] * nDims
    # rv2 = multivariate_normal(mean2)
    # result = np.logaddexp(rv1.logpdf(theta), rv2.logpdf(theta))

    mean = [[3,3], [3,12], [12,3], [12,12]]
    rv1 = multivariate_normal(mean[0])
    rv2 = multivariate_normal(mean[1])
    rv3 = multivariate_normal(mean[2])
    rv4 = multivariate_normal(mean[3])

    result1 = np.logaddexp(rv1.logpdf(theta), rv2.logpdf(theta))
    result2 = np.logaddexp(rv3.logpdf(theta), rv4.logpdf(theta))
    result = np.logaddexp(result1, result2)
    
    return result, []

def prior(hypercube):
    """ Uniform Prior for [Tmin, Tmax]. """
    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = UniformPrior(Tmin, Tmax)(x)

    return theta

# Define PolyChord settings
settings = PolyChordSettings(nDims, nDerived)
settings.do_clustering = True
settings.nlive = 100
settings.file_root = 'multidimensional'
settings.read_resume = False

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)


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
