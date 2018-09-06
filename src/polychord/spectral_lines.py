#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: spectral_lines.py

import numpy as np
import PyPolyChord as PPC
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior

# Importing data
path = "../../data/spectral_lines.txt"
data = np.loadtxt(path, delimiter=" ")

channel = data[:, 0]
signal = data[:, 1]

Tmin = 0.1
Tmax = 100

nDims = 1
nDerived = 0


def model(T, nu, nu0=37, sigmaL=2.0):
    """ Model 1 for the spectral lines model. """
    return T * np.exp(-(nu - nu0)**2 / (2 * sigmaL**2))


def logLikelihood(theta):
    """ Analytic expression of the log of the Likelihhod for the spectral lines problem."""
    sigma = 1.0
    N = len(data)

    exponent = -1 * np.sum((signal - model(theta, channel))**2) / (2*sigma**2)
    logL = (-N/2.)*np.log(2*np.pi) - N*np.log(sigma) + exponent
    return logL, []


def prior(hypercube):
    " Uniform Prior for [Tmin, Tmax]. "
    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = UniformPrior(Tmin, Tmax)(x)

    return theta


settings = PolyChordSettings(nDims, nDerived)
settings.nlive = 5000
settings.file_root = 'gregory'
settings.do_clustering = True
settings.read_resume = False

output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)

print(f'Z = {np.exp(output.logZ)}')

try:
    import getdist.plots
    import matplotlib.pyplot as plt
    # plt.close('all')
    posterior = output.posterior
    g = getdist.plots.getSubplotPlotter()
    # print("hola")
    g.triangle_plot(posterior, filled=True)
    plt.show()
except ImportError:
    print("Install matplotlib and getdist for plotting examples")
