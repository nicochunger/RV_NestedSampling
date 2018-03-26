#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ns_test.py

import numpy as np
import matplotlib.pyplot as plt

path = "/home/nunger/tesis/codigo/data/spectral_lines.txt"
data = np.loadtxt(path, delimiter=" ")

channel = data[:,0]
signal = data[:,1]

Tmin = 0.01
Tmax = 100

# --------------------- Function definitions -----------------------------
def model(T, nu, nu0=37, sigmaL=2):
    """" Model 1 for the spectral lines model. """
    return T * np.exp(-(nu - nu0)**2 / (2 * sigmaL**2))

def Uniform(Vmin, Vmax):
    """ This function creates a uniform prior for a certain range of values."""
    return 1/(Vmax - Vmin)

def Likelihood(T, data, grid, sigma=1):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    N = len(data)
    exponent = sum((data - model(T, grid))**2)
    return (2*np.pi)**(-N/2.) * sigma**(-N) * np.exp(-exponent / (2*sigma**2))

def CDF(pdf, x):
    """ Calculates de cummulative distribution function for a certain
    probability density function. """
    N = len(pdf)
    res = np.zeros(N) # Array with stored results
    for i in range(N):
        res[i] = np.trapz(pdf[:i], x[:i])

    return res

def SampleDist(cdf, x, samples):
    """ Samples x values acording to a probability density function given by
    cdf. samples is an array of values sampled uniformly from [0,1)"""
    npts = len(samples) # Number of samples
    dist = np.zeros(npts) # Final distribution of x values following the cdf
    for i in range(npts):
        idx = (np.abs(cdf-samples[i])).argmin()
        dist[i] = x[idx]
    return dist
# --------------------------------------------------------------------------

npts = 2000
T = np.linspace(Tmin,Tmax,npts)
pdf = [Likelihood(t, signal, channel) for t in T]

cdf = CDF(pdf, T)
cdf /= max(cdf) # Normalization of the CDF

samples = np.random.uniform(size=5000)
dist = SampleDist(cdf, T, samples)

plt.hist(dist, bins=60)
plt.show()

# ------------------- Nested Sampling algorithm ----------------------------
# Definition of variables and objects
N = 10 # Number of active objects
Obj = [] # List of active objects
Samples = [] # All Samples
Lstar = 0 # Current Likelihhod constraint
w = 1 - np.exp(-1./N) # First width in prior mass
H = 0.0 # Information
Z = 0 # Evidence, initially 0
Znew = Z # Updated Evidence
nest = 0 # Current iteration of the Nested Sampling loop
end = 2.0 # End condition for loop

# Initialization of first objects
for i in range(N):
    tsample = SampleDist(cdf, T, np.random.uniform()) # Sample 1 point from pdf
    Obj[i] = [tsample, Likelihood(tsample, signal, channel), 0] # Set objects

# Begin Nested Sampling loop
while nest <= end * N * H:
    # Search for worst Likelihhod withing active objects
    worst = 0
    for i in range(N):
        if(Obj[i][1] - Obj[worst][1]):
            worst = i
    Obj[worst][2] = w * Obj[worst][1] # Weight of object

    # Update Evidence and Information
    Znew = Z + Obj[worst][2]
    H = (Obj[worst][2] / Znew) * np.log(Obj[worst][1]) + (Z/ZNew)*(H+np.log(Z)) - np.log(Znew)
    Z = Znew

    Samples[nest] = Obj[worst]

    #Kill worst object in favour of a new exponent
    Lstar = Obj[worst][1]
    Obj[worst] = 



    nest += 1



# --------------------------------------------------------------------------



# # Plot of data
# plt.plot(channel, signal, 'k--')
# plt.plot(channel, signal, 'k*')
# plt.show()
