#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: ns_test.py

import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
import ns_library as ns
import scipy.stats as st

path = "../data/spectral_lines.txt"
data = np.loadtxt(path, delimiter=" ")

channel = data[:,0]
signal = data[:,1]

# Minimun and maximum values for the parameter
Tmin = 0.1
Tmax = 100

# --------------------- Function definitions -----------------------------
def model(T, nu, nu0=37, sigmaL=2.0):
    """ Model 1 for the spectral lines model. """
    return T * np.exp(-(nu - nu0)**2 / (2 * sigmaL**2))


def Likelihood(T, sigma=1.0):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    N = len(data)
    exponent = -np.sum((signal - model(T, channel))**2) / (2*sigma**2)
    lh = (2*np.pi)**(-float(N)/2.0) * sigma**(-float(N)) * np.exp(exponent)
    return lh

def logLikelihood(T, sigma=1.0):
    """ Analytic expression of the Likelihhod for the spectral lines problem."""
    N = len(data)
    exponent = -np.sum((signal - model(T, channel))**2) / (2*sigma**2)
    logL = (-N/2.)*np.log(2*np.pi) - N*np.log(sigma) + exponent
    return logL

class ActiveObj:
    def __init__(self, ppf, logL):
        self.param = 0.
        self.logLhood = 0.
        self.ppf = ppf
        self.logwt = 0.
        self.logLfunc = logL

    def Sample(self):
        """ Samples the object. """
        self.param = ns.SampleDist(self.ppf)
        #self.logLhood = logLikelihood(self.param)
        self.logLhood = self.logLfunc(self.param)

    def Evolve(self, Lconstraint, sampLimit):
        """ Evolves the Object to find a new sample given the
        Likelihood constraint. """
        count = 0
        limit = 100000
        while count <= limit:
            tsample = np.random.uniform(sampLimit[0], sampLimit[1])
            #tsample = ns.SampleDist(self.cdf)
            # tL = logLikelihood(tsample)
            tL = self.logLfunc(tsample)
            count += 1
            if tL > Lconstraint:
                self.param = tsample
                self.logLhood = tL
                return

        sys.exit("Evolve() couldn't find a sample with greater Likelihood \
                  within the sample limit constraint")

def bimodal(loc1, loc2, sigma1, sigma2):
    """ Creates a bimodal percent point function (inverse of cdf) with both
    mean values and standart deviation given as parameter.

    Returns a single variable function with the chosen bimodal distribution."""
    def foo(x):
        return np.log(st.norm.ppf(x, loc1, sigma1) + st.norm.ppf(x, loc2, sigma2))

    return foo

# --------------------------------------------------------------------------

npts = 5000
T = np.linspace(Tmin,Tmax,npts)
#pdf = [ns.Jeffreys(t, Tmin, Tmax) for t in T]
pdf = np.array([ns.Uniform(t, Tmin, Tmax) for t in T])
lanaly = np.array([Likelihood(t) for t in T]) * ns.Uniform(0, Tmin, Tmax)

# TODO implement analytic CDF for Uniform
#cdf = ns.CDF(pdf, T)
ppf = bimodal(3, 7, 0.3, 0.3)

# x = np.linspace(0,1,1000)
# dist = cdf(x)
# plt.figure()
# plt.plot(x, dist)
# plt.show()



# ------------------- Nested Sampling algorithm ----------------------------
def NestedSampling(N, priorPPF, logL, iterations=50, tol=1e-2):
    """ Runs the Nested Sampling algorithm for a certain amount of iterations
    and Calculates the resulting evidence each time.

    Parameters
    ----------
    N: int
        Number of live points to be used in the Nested Sampling algortihm.
    iterations: int (default = 50)
        How many times the algorithm is run. Saving all the results.
    tol: float (default = 1e-2)
        Tolerance for the termination condition. Nested Sampling loop ends when
        the estimated remaining evidence is less than the current calculates
        evidence time the tolerance.

    Returns
    -------
    evidences: ndarray
        Contains all the evidences calculated in each iteration.
    xi: ndarray
        Prior mass values for the last iteration.
    lvector: ndarray
        Contributions to the Likelihood of each xi for the last iteration.
    """
    evidences = []
    for i in range(iterations):
        # Initial sampling limits
        # sampLimit = [Tmin, Tmax]

        # Definition of variables and objects
        N_MAX = 50000 # Maximum samples
        Obj = [] # List of active objects
        Samples = [] # All Samples
        xi = [1] # Prior mass points (inicially 1)
        xi.append(ns.sampleXi(N, xi[-1])) # Calculate next xi to be one step ahead
        logw = np.log(1 - xi[-1]) # First width in prior mass
        H = 0.0 # Information
        logZ = -sys.float_info.max # log(Evidence, initially 0)
        logZnew = logZ # Updated Evidence
        nest = 0 # Current iteration of the Nested Sampling loop

        # Initialization of first objects
        for i in range(N):
            Obj.append(ActiveObj(priorPPF, logL)) # Creates an Active Object
            Obj[i].Sample() # Samples it

        # ------Begin Nested Sampling loop---------
        termination = False
        while not termination: #end * N * H:
            # Search for worst Likelihood within the active objects
            lhoods = np.array([obj.logLhood for obj in Obj])
            worst = np.argmin(lhoods)
            currZ = logw + Obj[worst].logLhood # Weight of object
            Obj[worst].logwt = currZ

            # Update Evidence and Information
            logZnew = ns.PLUS(logZ, currZ)
            log1 = Obj[worst].logLhood
            H = np.exp(currZ - logZnew) * log1 + np.exp(logZ - logZnew) + (H+logZ) - logZnew
            logZ = logZnew

            # # Print current data every 10 iteration
            # if nest % 10 == 0:
            #     print "logZ = {} \t logL = {} \t n = {}".format(logZ, log1, nest)

            # Save all chosen samples
            Samples.append(copy.deepcopy(Obj[worst]))

            # Adjust sampling limits
            params = [obj.param for obj in Obj]
            print params
            print lhoods
            sampLimit = [min(params), max(params)]
            # idx = np.abs(sampLimit - Obj[worst].param).argmin()
            # sampLimit[idx] = Obj[worst].param

            #Kill worst object in favour of a new object
            Lstar = Obj[worst].logLhood # Update Likelihood constraint
            Obj[worst].Evolve(Lstar, sampLimit) # Evolve the old sample for a new one

            # Update next prior mass value
            # xi is always one step ahead to calculate wi with the trapezoidal rule
            xi.append(ns.sampleXi(N, xi[-1]))
            logw = np.log((xi[-3] - xi[-1])/2)

            # Update termination condition
            termination = xi[nest]*np.exp(np.max(lhoods)) < tol * np.exp(logZ)

            # Increment iteration number
            nest += 1

            # Break loop if nest exeeds the maximum value
            if nest >= N_MAX:
                print("Loop exceeded maximum iteration number of {}".format(N_MAX))
                break
        #--------End of Nested Sampling loop----------------

        # Final correction
        logw = -float(nest)/N - np.log(float(N))
        final_corr = -sys.float_info.max
        for obj in Obj:
            obj.logwt = logw + obj.logLhood
            final_corr = ns.PLUS(final_corr, obj.logwt)
            logZnew = ns.PLUS(logZ, obj.logwt)
            H = np.exp(obj.logwt - logZnew) * obj.logLhood + np.exp(logZ - logZnew) + (H+logZ) - logZnew
            logZ = logZnew

        evidences.append(np.exp(logZ))

        # Plotting of solution
        xi = np.array(xi[:-2], dtype=float)
        lvector = np.array([np.exp(obj.logLhood) for obj in Samples], dtype=float)

        xsamples = np.array([obj.param for obj in Samples])

        # plt.plot(xsamples, '.')

        # plt.figure()
        # plt.plot(xi, lvector)
        # plt.xscale('log')

        #plt.show()
        # --------------------------------------------------------------------------


        # # Plot of data
        #plt.figure()
        #plt.plot(channel, signal, 'k--*')
        # plt.show()
    print("Average evidence for N={} and {} iterations: {}".format(N,iterations,np.mean(evidences)))
    return np.array(evidences), xi, lvector
