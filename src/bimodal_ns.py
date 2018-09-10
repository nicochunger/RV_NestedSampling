#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: bimodal_ns.py

from ns_library import Uniform, CDF
from ns_test import bimodal, NestedSampling
import numpy as np
import matplotlib.pyplot as plt
import time


# Generate a bimodal likelihood
likelihood = bimodal(3, 0.4, 7, 0.4)  # Note: this is the log of the likelihood
xmin = 0
xmax = 10
x = np.linspace(xmin, xmax, 1000)
y = np.exp(likelihood(x))

# # Plot the likelihood
# plt.plot(x, y)
# plt.show()

# Prior
pdf = np.array([Uniform(t, xmin, xmax) for t in x])
prior = CDF(pdf, x)

# ------ Calculation of analytic result ------
numeric_vec = y*Uniform(1, xmin, xmax)
evidence = np.trapz(numeric_vec, x)

print("Analytic result = {}".format(evidence))
# --------------------------------------------

start_time = time.time()
# Main result
result = NestedSampling(200, prior, likelihood, 100)

# Print only the evidences
# print(result[0])

stop_time = time.time()
delta_time = stop_time - start_time
print("Run time: {}:{} minutes or {:.3f} seconds".format(
    int(delta_time/60), int(delta_time % 60), delta_time))
