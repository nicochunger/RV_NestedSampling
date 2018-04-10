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

Tmin = 0.1
Tmax = 100

#np.random.seed(1234)
# --------------------- Function definitions -----------------------------
def model(T, nu, nu0=37, sigmaL=2.0):
    """ Model 1 for the spectral lines model. """
    return T * np.exp(-(nu - nu0)**2 / (2 * sigmaL**2))


def Likelihood(T, sigma=1.0):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    N = len(data)
    exponent = -np.sum((signal - model(T, channel))**2) / (2*sigma**2)
    lh = (2*np.pi)**(-float(N)/2.0) * sigma**(-float(N)) * np.exp(exponent)
    #if lh == 0.0:
    #   lh = 1e-308
    return lh

def logLikelihood(T, sigma=1.0):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    # N = len(data)
    # exponent = -np.sum((signal - model(T, channel))**2) / (2*sigma**2)
    # lh = (2*np.pi)**(-N/2.) * sigma**(-N) * np.exp(exponent)
    return np.log(Likelihood(T)) #np.log(lh)

class ActiveObj:
    def __init__(self, cdf):
        self.param = 0
        self.logLhood = 0
        self.cdf = cdf
        self.logwt = 0

    def Sample(self):
        """ Samples the object. """
        self.param = ns.SampleDist(self.cdf)
        self.logLhood = logLikelihood(self.param)

    def Evolve(self, Lconstraint):
        """ Evolves the Object to find a new sample given the
        Likelihood constraint. """
        count = 0
        limit = 100000
        while count <= limit:
            tsample = ns.SampleDist(self.cdf)
            tL = logLikelihood(tsample)
            count += 1
            if tL > Lconstraint:
                self.param = tsample
                self.logLhood = tL
                return

        sys.exit("Evolve() couldn't find a sample with greater Likelihood \
                  within the sample limit constraint")

# --------------------------------------------------------------------------

npts = 5000
T = np.linspace(Tmin,Tmax,npts)
#pdf = [ns.Jeffreys(t, Tmin, Tmax) for t in T]
pdf = np.array([ns.Uniform(t, Tmin, Tmax) for t in T])
lanaly = np.array([Likelihood(t) for t in T]) * ns.Uniform(0, Tmin, Tmax)
print "Analytic integration: {}\n".format(np.trapz(lanaly,T))

# TODO implement analytic CDF for Uniform
cdf = ns.CDF(pdf, T)

# x = np.linspace(0,1,1000)
# dist = cdf(x)
# plt.figure()
# plt.plot(x, dist)
# plt.show()

# ------------------- Nested Sampling algorithm ----------------------------
# Definition of variables and objects
N = 800 # Number of active objects
N_MAX = 10000 # Maximum samples
Obj = [] # List of active objects
Samples = [] # All Samples
logw = np.log(1 - np.exp(-1./N)) # First width in prior mass
xi = [1] # Prior mass points (inicially 1)
xi.append(ns.sampleXi(N, xi[-1])) # Calculate next xi to be one step ahead
H = 0.0 # Information
logZ = -sys.float_info.max * sys.float_info.epsilon # log(Evidence, initially 0)
logZnew = logZ # Updated Evidence
nest = 0 # Current iteration of the Nested Sampling loop
end = 2.0 # End condition for loop

# Initialization of first objects
for i in range(N):
    Obj.append(ActiveObj(cdf)) # Creates an Active Object
    Obj[i].Sample() # Samples it

# Begin Nested Sampling loop
# TODO implement better termination condition
tol = 1e-2
termination = False
while not termination: #end * N * H:
    # Search for worst Likelihood within the active objects
    lhoods = np.array([Obj[i].logLhood for i in range(N)])
    worst = np.argmin(lhoods)
    currZ = logw + Obj[worst].logLhood # Weight of object
    Obj[worst].logwt = currZ

    # Update Evidence and Information
    logZnew = ns.PLUS(logZ, currZ)
    log1 = Obj[worst].logLhood
    H = np.exp(currZ - logZnew) * log1 + np.exp(logZ - logZnew) + (H+logZ) - logZnew
    logZ = logZnew

    # Print current data every 10 iteration
    if nest % 10 == 0:
        print "logZ = {} \t H = {:.3f} \t logL = {} \t n = {}".format(logZ, H, log1, nest)

    # Save all chosen samples
    Samples.append(copy.deepcopy(Obj[worst]))

    #Kill worst object in favour of a new object
    Lstar = Obj[worst].logLhood # Update Likelihood constraint
    Obj[worst].Evolve(Lstar) # Evolve the old sample for a new one

    # Update next prior mass value
    # xi is always one step ahead to calculate wi with the trapezoidal rule
    xi.append(ns.sampleXi(N, xi[-1]))
    logw = np.log((xi[-3] - xi[-1])/2)

    # Update termination condition
    termination = xi[nest]*np.exp(np.mean(lhoods)) < tol * np.exp(logZ) #np.exp(logZ)

    # Increment iteration number
    nest += 1

    # Break loop if nest exeeds the maximum value
    if nest >= N_MAX:
        print "Loop exceeded maximum iteration number of {}".format(N_MAX)
        break

# Final correction
logw = -float(nest)/N - np.log(float(N))
final_corr = 0
for obj in Obj:
    obj.logwt = logw + obj.logLhood
    final_corr += obj.logwt
    logZnew = ns.PLUS(logZ, obj.logwt)
    H = np.exp(obj.logwt - logZnew) * obj.logLhood + np.exp(logZ - logZnew) + (H+logZ) - logZnew
    logZ = logZnew
print "Final correction: {}\n".format(np.exp(final_corr))

print "Iterations: {}".format(nest)
print "Final evidence: {}".format(np.exp(logZ))

# Plotting of solution
xi = np.array(xi[:-2])
lvector = np.array([np.exp(obj.logLhood) for obj in Samples])

plt.figure()
plt.plot(xi, lvector)
plt.xscale('log')

#plt.show()
# --------------------------------------------------------------------------


# # Plot of data
#plt.figure()
#plt.plot(channel, signal, 'k--*')
plt.show()
