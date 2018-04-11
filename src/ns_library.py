#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ns_library.py

import numpy as np
from scipy import interpolate
from scipy.integrate import quad

# IDEA change the priors functions to return a single variable function instead
# of the already evaluated value.

def Uniform(x, xmin, xmax):
    """ This function creates a uniform prior for a certain range of values."""
    assert xmin < xmax, "xmin has to be lower than xmax."
    return 1./(xmax - xmin)

def Jeffreys(x, xmin, xmax):
    """ This function returns the Jeffreys prior for a certain value. """
    assert xmax > 0, "xmax hasta to be greater than 0."
    assert xmax > xmin, "xmax hasta to be greater than xmin."
    def prior(x):
        return 1./(x * np.log(float(xmax)/float(xmin)))
    return prior #1./(x * np.log(float(xmax)/float(xmin)))

def CDF(pdf, x):
    """ Calculates the cummulative distribution function for a certain
    probability density function.

    Parameters
    ----------
    pdf : np_array
          Numerical array with the probability density function.
    x : np_array
        x values of the pdf.

    Returns
    -------
    function
        Interpolation of the inverse of the cummulative distribution function.
    """
    # TODO implement use of analytic pdf
    N = len(pdf)
    res = np.zeros(N) # Array with stored results
    for i in range(N):
        res[i] = np.trapz(pdf[:i], x[:i])

    res /= max(res) # Normalization of cdf
    return interpolate.interp1d(res, x) # Interpolation of inverse of the cdf

def SampleDist(cdf, N=1):
    """ Samples values acording to a probability density function given by the
    inverse of the cdf. N is the number of values to be sampled, default is 1.
    Returns a numpy array.

    Parameters
    ----------
    cdf : function
          Inverse cummulative density function.
    N : int
        Number of samples to be returned. (default is 1)

    Returns
    -------
    ndarray
            List of samples.
    """
    samples = np.random.uniform(size=N) # N random numbers in [0,1)
    dist = cdf(samples) # Final distribution of x values following the cdf
    return dist

def PLUS(x, y):
    """ Logarithmic sum. """
    if x > y:
        return x + np.log(1+np.exp(y-x))
    else:
        return y + np.log(1+np.exp(x-y))

def sampleXi(N, xi):
    """ Samples the next prior mass value. """
    samples = np.random.uniform(size=N)
    return xi * np.max(samples)
