#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: estimator_evolution.py

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from ns_test import NestedSampling
import time

start_time = time.time()

results = [] # List of arrays with all the evidences
prior_masses = [] # List of arrays with the prior masses of the last iterations
likelihoods = [] # List of arrays with the likelihoods of the last iterations

# Number of live points to be analized
N = np.arange(1100, 2001, 100)

# Main loop
for i in range(len(N)):
    result, priormass, lh = NestedSampling(N[i])
    results.append(result)
    prior_masses.append(priormass)
    likelihoods.append(lh)

print("Run time: {} seconds".format(time.time() - start_time))

# Plotting
plt.boxplot(results, labels=N)
plt.xlabel("Live points")
plt.ylabel("Evidence")
plt.title("Evidence distribution for different number of live points.")
plt.show()
