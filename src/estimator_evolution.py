#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: estimator_evolution.py

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from ns_test import NestedSampling
import time

start_time = time.time()


# Number of live points to be analized
Nmin = 300
Nmax = 501
step = 200
N = np.arange(Nmin, Nmax, step)
Nsize = len(N)

niter = 500

# Preallocation of arrays
results = [] # List of arrays with all the evidences
prior_masses = [] # List of arrays with the prior masses of the last iterations
likelihoods = [] # List of arrays with the likelihoods of the last iterations

# Main loop
for i in range(len(N)):
    result, priormass, lh = NestedSampling(N[i], iterations=niter)
    results.append(result)
    prior_masses.append(priormass)
    likelihoods.append(lh)

stop_time = time.time()
delta_time = stop_time - start_time
print("Run time: {}:{}".format(int(delta_time/60),int(delta_time%60)))

# Save data to txt files
N_text = "{}_{}_{}".format(Nmin, Nmax, step)
path = "../../graficos/nested_sampling_test/iter{}/".format(niter)
results_save_loc = path + "results_" + N_text + ".txt"
np.savetxt(results_save_loc, results)

for i in range(Nsize):
    xi_save_loc = path + "xi_" + str(N[i]) + ".txt"
    np.savetxt(xi_save_loc, prior_masses[i])

    lhoods_save_loc = path + "lhoods_" + str(N[i]) + ".txt"
    np.savetxt(lhoods_save_loc, likelihoods[i])


# Plotting
plt.figure(1)
plt.boxplot(results, labels=N, showfliers=False)
plt.xlabel("Live points")
plt.ylabel("Evidence")
plt.title("Evidence distribution for different number of live points.")

plt.figure(2)
for i in range(1,Nsize,2):
    plt.plot(prior_masses[i], likelihoods[i], label="N = "+str(N[i]))
plt.legend()
plt.xscale('log')
plt.xlabel('Prior mass')
plt.ylabel('Likelihood')

plt.show()
