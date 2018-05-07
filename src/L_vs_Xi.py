#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file: L_vs_Xi.py

import numpy as np
from ns_test import Likelihood, logLikelihood
from scipy.special import erfinv
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt



# Likelihood test
T = np.linspace(0.1, 100, 200)
Lhood = Likelihood(T)
integral = np.trapz(Lhood, T)
#print(integral)

L_norm = Lhood #/ integral
int_norm = np.trapz(L_norm, T)
#print(int_norm)

# Gaussian model
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

p0 = [1.0e-38, 1.5, 1.]

coeff, var_matrix = curve_fit(gauss, T, L_norm, p0=p0) # Fit data
mu = coeff[1]
sigma = coeff[2]
print(coeff)

# Prior mass
xi = np.logspace(-2, 0, 200)

# Inverse error function
u = erfinv(xi)

x = 2*u*sigma + mu
L = Likelihood(x)

plt.plot(xi, L, '.')
plt.xscale('log')
plt.show()


