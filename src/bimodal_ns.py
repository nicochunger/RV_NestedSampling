#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: bimodal_ns.py

from ns_library import *
from ns_test import *
import numpy as np
import matplotlib.pyplot as plt

# FIXME the likelihood gives nan on almost all values.
likelihood = bimodal(3, 0.3, 7, 0.3)
xmin = 0
xmax = 10
x = np.linspace(xmin, xmax, 1000)
pdf = np.array([Uniform(t, xmin, xmax) for t in x])
prior = CDF(pdf, x)

numeric_vec = likelihood(x)/Uniform(1, xmin, xmax)
evidence = np.trapz(numeric_vec, x)

#
print("Analytic result = {}".format(evidence))

result = NestedSampling(100, prior, likelihood, 10)

print(result[0])
