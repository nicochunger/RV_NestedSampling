#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ns_test.py

import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ns_library as ns

path = "/home/nunger/tesis/codigo/data/spectral_lines.txt"
data = np.loadtxt(path, delimiter=" ")

channel = data[:,0]
signal = data[:,1]

Tmin = 0.01
Tmax = 100

#np.random.seed(1234)
# --------------------- Function definitions -----------------------------
def model(T, nu, nu0=37, sigmaL=2):
    """" Model 1 for the spectral lines model. """
    return T * np.exp(-(nu - nu0)**2 / (2 * sigmaL**2))


def Likelihood(T, data, grid, sigma=1):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    N = len(data)
    exponent = np.sum((data - model(T, grid))**2)
    lh = (2*np.pi)**(-N/2.) * sigma**(-N) * np.exp(-exponent / (2*sigma**2))
    res = lh * ns.Uniform(Tmin, Tmax)
    if res == 0.0:
       res = 1e-308
    return res

class ActiveObj:
    def __init__(self, T, cdf, xdata, ydata):
        self.param = 0
        self.Lhood = 0
        self.xdata = xdata
        self.ydata = ydata
        self.cdf = cdf
        self.T = T
        self.wt = 0

    def Sample(self):
        """ Samples the object. """
        self.param = ns.SampleDist(self.cdf, self.T)
        self.Lhood = Likelihood(self.param, self.ydata, self.xdata)

    def Evolve(self, Lconstraint):
        """ Evolves the Object to find a new sample given the
        Likelihood constraint. """
        while 1:
            tsample = ns.SampleDist(self.cdf, self.T, 1)
            tL = Likelihood(tsample, self.ydata, self.xdata)
            if tL > Lconstraint:
                self.param = tsample
                self.Lhood = tL
                break

# --------------------------------------------------------------------------

npts = 5000
T = np.linspace(Tmin,Tmax,npts)
#pdf = [ns.Jeffreys(t, Tmin, Tmax) for t in T]
pdf = [ns.Uniform(Tmin, Tmax) for t in T]
#print "Analytic integration: {}\n".format(np.trapz(pdf,T))

cdf = ns.CDF(pdf, T)

# x = np.linspace(0,1,1000)
# dist = cdf(x)
# plt.figure()
# plt.plot(x, dist)
# plt.show()

# ------------------- Nested Sampling algorithm ----------------------------
# Definition of variables and objects
N = 100 # Number of active objects
N_MAX = 10000 # Maximum samples
Obj = [] # List of active objects
Samples = [] # All Samples
Lstar = 0 # Current Likelihhod constraint
w = 1 - np.exp(-1./N) # First width in prior mass
xi = [] # Prior mass
H = 0.0 # Information
Z = sys.float_info.min # Evidence, initially 0
Znew = Z # Updated Evidence
nest = 0 # Current iteration of the Nested Sampling loop
end = 2.0 # End condition for loop


# Initialization of first objects
for i in range(N):
    Obj.append(ActiveObj(T, cdf, channel, signal)) # Creates an Active Object
    Obj[i].Sample() # Samples it

# Begin Nested Sampling loop
while nest <= 700: #end * N * H:
    # Search for worst Likelihood within the active objects
    lhoods = [obj.Lhood for obj in Obj]
    worst = np.argmin(lhoods)
    currZ = w * Obj[worst].Lhood # Weight of object
    Obj[worst].wt = currZ

    # Update Evidence and Information
    Znew = Z + currZ
    log1 = np.log(Obj[worst].Lhood)
    log2 = np.log(Z)
    log3 = np.log(Znew)
    H = (currZ / Znew) * log1 + (Z/Znew)*(H+log2) - log3
    #H = (currZ / Znew) * np.log(Obj[worst].Lhood) + (Z/Znew)*(H+np.log(Z)) - np.log(Znew)
    Z = Znew
    print "Z = {} \t H = {:.3f} \t L = {} \t n = {}".format(Z, H, Obj[worst].Lhood, nest)

    # Save all chosen samples
    Samples.append(copy.deepcopy(Obj[worst]))

    #Kill worst object in favour of a new object
    Lstar = Obj[worst].Lhood # Update Likelihood constraint
    Obj[worst].Evolve(Lstar) # Evolve the old sample for a new one

    xi.append(np.exp(-float(nest)/N))
    nest += 1
    w = np.exp(-float(nest)/N) - np.exp(-float(nest+1)/N) # Shrink weight size

    # Break loop if nest exeeds the maximum value
    if nest >= N_MAX:
        print "Loop exceeded maximum iteration number of {}".format(N_MAX)
        break

# Final correction


print "Iterations: {}".format(nest)
print "Final evidence: {}".format(Z)

# Plotting of solution
lvector = [obj.Lhood for obj in Samples]

plt.figure()
plt.plot(xi, lvector)
plt.xscale('log')

plt.show()
# --------------------------------------------------------------------------


# # Plot of data
# plt.plot(channel, signal, 'k--')
# plt.plot(channel, signal, 'k*')
# plt.show()
