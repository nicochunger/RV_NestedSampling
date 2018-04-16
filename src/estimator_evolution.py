#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: estimator_evolution.py

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from ns_test import ns_algorithm

results = []

N = np.arange(100, 1001, 100)

for i in range(len(N)):
    result = ns_algorithm(N[i])
    results.append(result)

#print(results)

plt.boxplot(results, labels=N)
plt.xlabel("Live points")
plt.ylabel("Evidence")
plt.title("Evidence distribution for different number of live points.")
plt.show()
