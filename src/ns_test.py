#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ns_test.py

import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ns_library as ns

path = "../data/spectral_lines.txt"
data = np.loadtxt(path, delimiter=" ")

channel = data[:,0]
signal = data[:,1]

Tmin = 0.01
Tmax = 100

#np.random.seed(1234)
# --------------------- Function definitions -----------------------------
def model(T, nu, nu0=37, sigmaL=2.0):
    """" Model 1 for the spectral lines model. """
    return T * np.exp(-(nu - nu0)**2 / (2 * sigmaL**2))


def Likelihood(T, data, grid, sigma=1.0):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    N = len(data)
    exponent = np.sum((data - model(T, grid))**2)
    lh = (2*np.pi)**(-N/2.) * sigma**(-N) * np.exp(-exponent / (2*sigma**2))
    #if lh == 0.0:
    #   lh = 1e-308
    return lh

def logLikelihood(T, sigma=1.0):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    N = len(data)
    exponent = -np.sum((signal - model(T, channel))**2) / (2*sigma**2)
    lh = (2*np.pi)**(-N/2.) * sigma**(-N) * np.exp(exponent)
    return np.log(lh)

class ActiveObj:
    def __init__(self, T, cdf):
        self.param = 0
        self.logLhood = 0
        self.cdf = cdf
        self.T = T
        self.logwt = 0

    def Sample(self):
        """ Samples the object. """
        self.param = ns.SampleDist(self.cdf, self.T)
        self.logLhood = logLikelihood(self.param)

    def Evolve(self, Lconstraint):
        """ Evolves the Object to find a new sample given the
        Likelihood constraint. """
        while 1:
            tsample = ns.SampleDist(self.cdf, self.T, 1)
            tL = logLikelihood(tsample)
            if tL > Lconstraint:
                self.param = tsample
                self.logLhood = tL
                break

# --------------------------------------------------------------------------

npts = 5000
T = np.linspace(Tmin,Tmax,npts)
#pdf = [ns.Jeffreys(t, Tmin, Tmax) for t in T]
pdf = [ns.Uniform(t, Tmin, Tmax) for t in T]
#print "Analytic integration: {}\n".format(np.trapz(pdf,T))

# TODO implement analytic CDF for Uniform
cdf = ns.CDF(pdf, T)

# x = np.linspace(0,1,1000)
# dist = cdf(x)
# plt.figure()
# plt.plot(x, dist)
# plt.show()

# ------------------- Nested Sampling algorithm ----------------------------
# Definition of variables and objects
N = 300 # Number of active objects
N_MAX = 10000 # Maximum samples
Obj = [] # List of active objects
Samples = [] # All Samples
logw = np.log(1 - np.exp(-1./N)) # First width in prior mass
xi = [] # Prior mass
H = 0.0 # Information
logZ = -sys.float_info.max # log(Evidence, initially 0)
logZnew = logZ # Updated Evidence
nest = 0 # Current iteration of the Nested Sampling loop
end = 2.0 # End condition for loop

# Initialization of first objects
for i in range(N):
    Obj.append(ActiveObj(T, cdf)) # Creates an Active Object
    Obj[i].Sample() # Samples it

# Begin Nested Sampling loop
while nest <= 7 * N: #end * N * H:
    # Search for worst Likelihood within the active objects
    lhoods = np.zeros([N])
    for i in range(N):
        lhoods[i] = Obj[i].logLhood
    worst = np.argmin(lhoods)
    currZ = logw + Obj[worst].logLhood # Weight of object
    Obj[worst].logwt = currZ

    # Update Evidence and Information
    logZnew = ns.PLUS(logZ, currZ)
    log1 = Obj[worst].logLhood
    # log2 = logZ
    # log3 = logZnew
    #H = (currZ / Znew) * log1 + (logZ/Znew)*(H+log2) - log3
    #H = (currZ / Znew) * np.log(Obj[worst].Lhood) + (Z/Znew)*(H+np.log(Z)) - np.log(Znew)

    H = np.exp(currZ - logZnew) * log1 + np.exp(logZ - logZnew) + (H+logZ) - logZnew

    logZ = logZnew
    print "logZ = {} \t H = {:.3f} \t logL = {} \t n = {}".format(logZ, H, log1, nest)

    # Save all chosen samples
    Samples.append(copy.deepcopy(Obj[worst]))

    #Kill worst object in favour of a new object
    Lstar = Obj[worst].logLhood # Update Likelihood constraint
    Obj[worst].Evolve(Lstar) # Evolve the old sample for a new one

    xi.append(np.exp(-float(nest)/N))
    nest += 1
    #w = np.exp(-float(nest-1)/N) - np.exp(-float(nest)/N) # Shrink weight size
    logw -= 1.0/N

    # Break loop if nest exeeds the maximum value
    if nest >= N_MAX:
        print "Loop exceeded maximum iteration number of {}".format(N_MAX)
        break

# Final correction


print "Iterations: {}".format(nest)
print "Final evidence: {}".format(np.exp(logZ))

# Plotting of solution
lvector = [np.exp(obj.logLhood) for obj in Samples]
lvector = np.array(lvector)

plt.figure()
plt.plot(xi, lvector)
plt.xscale('log')

#plt.show()
# --------------------------------------------------------------------------


# # Plot of data
#plt.figure()
#plt.plot(channel, signal, 'k--*')
#plt.plot(channel, signal, 'k*')
plt.show()
