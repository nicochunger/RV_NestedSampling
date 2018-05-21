#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_bimodal.py

import PyPolyChord as PPC 
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

nDims = 1
nDerived = 0

Tmin = 0
Tmax = 10

def logLikelihood(theta):
    """ log Likelihood of this toy bimodal example. """
    loc1 = 3; sigma1 = 0.4
    loc2 = 7; sigma2 = 0.4
    result = np.log(norm.pdf(theta, loc1, sigma1) + norm.pdf(theta, loc2, sigma2))
    
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
settings.nlive = 100
settings.file_root = 'bimodal'
settings.read_resume = False

# Run PolyChord
output = PPC.run_polychord(logLikelihood, nDims, nDerived, settings, prior)
