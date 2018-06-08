#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ppc_eprv3.py

import sys, os
path = '/home/nunger/tesis/codigo/run'
sys.path.append(path)

from model_eprv3_kepler import lnlike, lnprior, preprocess
import numpy as np

# Import Data

nDims = planets * len(params) # number of planets time the number of parameters for each planet

# Generate dictionaries


def logLikelihood(theta):
    """ Log Likelihood for the planet model. """


    return result, []


def prior(theta):
    """ Priors for each parameter. """
    theta = [0.0] * nDims
    for i, x in enumerate(hypercube):
        theta[i] = UniformPrior(Tmin, Tmax)(x)

    return theta



