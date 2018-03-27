#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ns_test.py

import numpy as np
import matplotlib.pyplot as plt
import sys
import pprint as pp

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
    return 1./(Vmax - Vmin)

def Likelihood(T, data, grid, sigma=1):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    N = len(data)
    exponent = sum((data - model(T, grid))**2)
    lh = (2*np.pi)**(-N/2.) * sigma**(-N) * np.exp(-exponent / (2*sigma**2))
    return Uniform(Tmin, Tmax) * lh

def CDF(pdf, x):
    """ Calculates de cummulative distribution function for a certain
    probability density function. """
    N = len(pdf)
    res = np.zeros(N) # Array with stored results
    for i in range(N):
        res[i] = np.trapz(pdf[:i], x[:i])

    return res

def SampleDist(cdf, x, N):
    """ Samples x values acording to a probability density function given by
    cdf. N is the number of values to be sampled. Returns a numpy array."""
    samples = np.random.uniform(size=N) # N random numbers in [0,1)
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
print "Tiene que dar: {}".format(np.trapz(pdf,T))

cdf = CDF(pdf, T)
cdf /= max(cdf) # Normalization of the CDF

#dist = SampleDist(cdf, T, 5000)

# plt.hist(dist, bins=60)
# plt.show()

# ------------------- Nested Sampling algorithm ----------------------------
# Definition of variables and objects
N = 10 # Number of active objects
N_MAX = 10000 # Maximum samples
Obj = [] # List of active objects
Samples = [] # All Samples
Lstar = 0 # Current Likelihhod constraint
w = 1 - np.exp(-1./N) # First width in prior mass
H = 0.0 # Information
Z = sys.float_info.min #-sys.float_info.max * sys.float_info.epsilon # Evidence, initially 0
Znew = Z # Updated Evidence
nest = 0 # Current iteration of the Nested Sampling loop
end = 2.0 # End condition for loop

# Initialization of first objects
for i in range(N):
    tsample = SampleDist(cdf, T, 1) # Sample 1 point from pdf
    Obj.append([tsample, Likelihood(tsample, signal, channel)]) # Set objects

# Begin Nested Sampling loop
while nest <= 50: #end * N * H:
    # Search for worst Likelihood within the active objects
    worst = 0
    for i in range(N):
        if(Obj[i][1] < Obj[worst][1]):
            worst = i
    currZ = w * Obj[worst][1] # Weight of object

    # Update Evidence and Information
    Znew = Z + currZ
    H = (currZ / Znew) * np.log(Obj[worst][1]) + (Z/Znew)*(H+np.log(Z)) - np.log(Znew)
    Z = Znew
    print "Z = {} \t H = {:.3f} \t n = {}".format(Z, H, nest)

    Samples.append(Obj[worst])

    #Kill worst object in favour of a new object
    Lstar = Obj[worst][1] # Update Likelihood constraint
    if N!=1: # Keep it if N=1
        while 1:
            tsample = SampleDist(cdf, T, 1)
            tL = Likelihood(tsample, signal, channel)
            if tL > Lstar:
                Obj[worst] = [tsample, tL, 0]
                break

    nest += 1
    # Break loop if nest exeeds the maximum value
    if nest >= N_MAX:
        break
    w = np.exp(-float(nest/N)) - np.exp(-float(nest+1)/N) # Shrink weight size

print "\nFinal evidende: {}".format(Z)


# --------------------------------------------------------------------------



# # Plot of data
# plt.plot(channel, signal, 'k--')
# plt.plot(channel, signal, 'k*')
# plt.show()
