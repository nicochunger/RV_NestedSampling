#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: estimator_evolution.py

import numpy as np
import matplotlib.pyplot as plt
from ns_test import logLikelihood, NestedSampling
from ns_library import Uniform, CDF
import time
import datetime

# Save start time to calculate the total run time
start_time = time.time()

# Prior information
Tmin = 0.1
Tmax = 100
npts = 5000
T = np.linspace(Tmin,Tmax,npts)
pdf = pdf = np.array([Uniform(t, Tmin, Tmax) for t in T])
cdf = CDF(pdf, T)

# Number of live points to be analized
Nmin = 7000
Nmax = 7001
step = 2
N = np.arange(Nmin, Nmax, step)
Nsize = len(N)

# Number of iterations of the Nested Sampling Algorithm for each N
niter = 50

# Preallocation of arrays
results = [] # List of arrays with all the evidences
prior_masses = [] # List of arrays with the prior masses of the last iterations
likelihoods = [] # List of arrays with the likelihoods of the last iterations

# Main loop
for i in range(Nsize):
    result, priormass, lh = NestedSampling(N[i], cdf, logLikelihood, iterations=niter)
    results.append(result)
    prior_masses.append(priormass)
    likelihoods.append(lh)

# Calculation of run time
stop_time = time.time()
delta_time = stop_time - start_time
print("Run time: " + str(datetime.timedelta(seconds=int(delta_time))))
# print("Run time: {}:{}:{}".format(int(delta_time/3600),int((delta_time%3600)/60),int(delta_time%60)))

# Save data to txt files
N_text = "{}_{}_{}".format(Nmin, Nmax, step)
path = "../../graficos/nested_sampling_test/"
results_save_loc = path + "results_" + N_text + ".txt"
np.savetxt(results_save_loc, results)

for i in range(Nsize):
    xi_save_loc = path + "xi_" + str(N[i]) + ".txt"
    np.savetxt(xi_save_loc, prior_masses[i])

    lhoods_save_loc = path + "lhoods_" + str(N[i]) + ".txt"
    np.savetxt(lhoods_save_loc, likelihoods[i])


# # Plotting
# plt.figure(1)
# plt.boxplot(results, labels=N, showfliers=False)
# plt.xlabel("Live points")
# plt.ylabel("Evidence")
# plt.title("Evidence distribution for different number of live points.")

# plt.figure(2)
# for i in range(1,Nsize,2):
#     plt.plot(prior_masses[i], likelihoods[i], label="N = "+str(N[i]))
# plt.legend()
# plt.xscale('log')
# plt.xlabel('Prior mass')
# plt.ylabel('Likelihood')

# plt.show()
