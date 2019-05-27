import pickle
import matplotlib.pyplot as plt
# from PyAstronomy.pyasl import foldAt
import pandas as pd
import numpy as np
import numpy.ma as ma
from pprint import pprint
import os
import getdist
from uncertainties import ufloat
from uncertainties.umath import *

# Load evidence results
output = pickle.load(open('output.p', 'rb'))
print(f'Evidence (logZ) = {output.logZ:.2f} +/- {output.logZerr:.2f}')
# Change direcory of posterior
output.base_dir = os.path.dirname(os.path.abspath(__file__))
posterior = output.posterior
# posterior.removeBurn(0.5)

# Construst "real" posterior
idxs = []
for i, x in enumerate(posterior.weights):
    if np.random.random() < x:
        idxs.append(i)

samples = posterior.samples[idxs]
print(f'Nr. of samples in posterior: {len(samples)}')

# Get the medians for each parameter
paramnames = posterior.getParamNames().list()
medians = dict(zip(paramnames, zip(np.round(np.median(
    samples, axis=0), 7), np.round(np.std(samples, axis=0), 7))))
names = ['c', 'd', 'b', 'f', 'g']
pprint(medians)

# --------------- POSTERIORS -----------------------------
# Count the number of planets in model
nplanets = 0
for i in range(10):
    if f'planet{i}_period' in paramnames:
        nplanets += 1

# Construct dictionary with posterior for each parameter
post = {}
for i in range(len(paramnames)):
    post.update({paramnames[i]: samples[:, i]})


def gauss(x, *p):
    # Gaussian model
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# fig2, ax2 = plt.subplots(4, 1)
# ax2[0].hist(post['drift1_lin']*1e4,
#             bins='auto', histtype='step', normed=True)
# ax2[1].hist(post['drift1_qua']*1e4,
#             bins='auto', histtype='step', normed=True)
# ax2[2].hist(post['drift1_cub']*1e4,
#             bins='auto', histtype='step', normed=True)
# ax2[3].hist(post['drift1_unitconstant'],
#             bins='auto', histtype='step', normed=True)
# plt.show()


# Plot posterior for the period of each planet
fig, ax = plt.subplots(1, nplanets, figsize=(12, 3))
ax = np.atleast_1d(ax)  # To support 1 planet models


def int_hist(hist, bin_edges):
    res = 0
    lim95 = 0
    counter1 = 1
    lim99 = 0
    counter2 = 1
    for i in range(len(hist)):
        width = bin_edges[i+1] - bin_edges[i]
        res += hist[i] * width
        if res > 0.95 and counter1:
            lim95 = bin_edges[i]
            counter1 = 0
        if res > 0.99 and counter2:
            lim99 = bin_edges[i]
            counter2 = 0

    return res, lim95, lim99


# MINIM MASS
def min_mass(K, period, ecc, mstar):
    """ 
    Calculates de minimum mass using known parameters for the planet and the star
    K [m/s]: radial velocity amplitud 
    period [days]: orbital period of planet
    ecc: eccentricity of planet
    mstar [Msun]: Mass of star in Sun masses
    """
    # Unit conversions
    period *= 86400 # days to seconds
    mstar *= 1.9891e30 # Sun masses to kg

    msini = K*(mstar)**(2/3) * (period/(2*np.pi*6.67e-11))**(1/3) * sqrt(1-ecc**2)

    # Convert to Earth masses
    msini /= 5.972e24
    return msini


def semi_mayor_axis(P, mstar):
    """ Calculates the semimayor axis for a planet.
    P [days]: orbital period of planet
    mstar [Msun]: Mass of star in Solar Masses """
    # Unit conversion
    P *= 86400 # days to seconds
    mstar *= 1.9891e30 # Msun to kg
    G = 6.67e-11

    sma = ((P**2 * G * mstar) / (4*np.pi**2))**(1/3)

    # In AU
    sma /= 1.496e+11

    return sma

# MAIN LOOP
for i in range(nplanets):
    period_post = post[f'planet{i+1}_period']
    median = medians[f'planet{i+1}_period'][0]

    # Histogram
    # hist, bin_edges, patches = ax[i].hist(
    #     period_post, label='Posterior', bins='fd', histtype='step', normed=True)
    hist, bin_edges, patches = ax[i].hist(
        period_post, label='Posterior', bins=40, histtype='step', normed=True)

    res, lim95, lim99 = int_hist(hist, bin_edges)
    # print(
    #     f"Median ({median:5.2f}): lim95 = {lim95:5f}; lim99 = {lim99:5f}; res = {res:5f}")

    # Graph posterior only for planet f (0.978)
    # if medians[f'planet{i+1}_period'][0] < 1:
    #     plt.figure()
    #     plt.hist(
    #         period_post, label='Posterior', bins="auto", histtype='step', normed=True)
    #     plt.xlabel('Período [d]')
    #     plt.ylabel('PDF')
    #     plt.title(f"Mediana = {medians[f'planet{i+1}_period'][0]:5.5f} días")

    # Time series of posterior
    # ax[i].plot(period_post, '.')

    ax[i].set_title(
        f"Mediana = {median:5.3f} días")
    ax[i].set_xlabel('Período [d]')
    if i == 0:
        ax[i].set_ylabel('PDF')

    # # Minimum mass calculation
    K = ufloat(medians[f'planet{i+1}_k1'][0], medians[f'planet{i+1}_k1'][1])
    period = ufloat(median, medians[f'planet{i+1}_period'][1])
    ecc = ufloat(medians[f'planet{i+1}_ecc'][0], medians[f'planet{i+1}_ecc'][1])
    mstar = ufloat(0.75, 0.035)
    print(f"Planet {period:.2f}, m*sin(i) = {min_mass(K, period, ecc, mstar):5f} Mearth;  a = {semi_mayor_axis(period, mstar)} AU")




plt.tight_layout()
plt.show()
# -------------------------------------------------------

# ---------- PHASE FOLDING -----------------
# phases = np.zeros((len(t), nplanets))
# for i in range(nplanets):
#     phases[:, i] = foldAt(t, medians[f'planet{i+1}_period'])

# # Sort with respect to phase
# # First, get the order of indices ...
# sortIndi = np.argsort(phases, axis=0)
# # ... and, second, rearrange the arrays.
# y_folded = np.zeros((len(t), nplanets))
# for i in range(nplanets):
#     phases[:, i] = phases[:, i][sortIndi[:, i]]
#     y_folded[:, i] = y[sortIndi[:, i]]

#     # Plot the result
#     plt.figure(i+1)
#     plt.plot(phases[:, i], y_folded[:, i], '.')
#     plt.title(f"Planet {i}. Period: {medians[f'planet{i+1}_period']}")
# plt.show()
# # -------------------------------------------
