#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ns_library.py

import numpy as np
from scipy import interpolate

def Uniform(x, xmin, xmax):
    """ This function creates a uniform prior for a certain range of values."""
    return 1./np.abs(xmax - xmin)

def Jeffreys(x, xmin, xmax):
    """ This function returns the Jeffreys prior for a certain value. """
    assert xmax > 0, "xmax hasta to be greater than 0."
    assert xmax > xmin, "xmax hasta to be greater than xmin."
    return 1./(x * np.log(xmax/xmin))

def CDF(pdf, x):
    """ Calculates de cummulative distribution function for a certain
    probability density function.
    Returns a function with the inverse of the cdf. """
    N = len(pdf)
    res = np.zeros(N) # Array with stored results
    for i in range(N):
        res[i] = np.trapz(pdf[:i], x[:i])

    res /= max(res)
    return interpolate.interp1d(res, x)

def SampleDist(cdf, x, N=1):
    """ Samples x values acording to a probability density function given by
    cdf. N is the number of values to be sampled, default is 1.
    Returns a numpy array."""
    samples = np.random.uniform(size=N) # N random numbers in [0,1)
    dist = cdf(samples) # Final distribution of x values following the cdf
    return dist

def PLUS(x, y):
    """ Logarithmic sum. """
    if x > y:
        return x + np.log(1+np.exp(y-x))
    elif y > x:
        return y + np.log(1+np.exp(x-y))
