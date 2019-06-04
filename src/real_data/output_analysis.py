import pickle
import matplotlib.pyplot as plt
from PyAstronomy.pyasl import foldAt
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
pprint(medians)

# --------------- POSTERIORS -----------------------------
# Count the number of planets in model
nplanets = 0
for i in range(10):
    if f'planet{i}_period' in paramnames:
        nplanets += 1

# Identify planet period with name
names = {4: ('b', 0), 10: ('c', 1), 20: ('d', 2), 1: ('f', 3), 131: ('f', 3)}

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
    period *= 86400  # days to seconds
    mstar *= 1.9891e30  # Sun masses to kg

    msini = K*(mstar)**(2/3) * \
        (period/(2*np.pi*6.67e-11))**(1/3) * sqrt(1-ecc**2)

    # Convert to Earth masses
    msini /= 5.972e24
    return msini


def semi_mayor_axis(P, mstar):
    """ Calculates the semimayor axis for a planet.
    P [days]: orbital period of planet
    mstar [Msun]: Mass of star in Solar Masses """
    # Unit conversion
    P *= 86400  # days to seconds
    mstar *= 1.9891e30  # Msun to kg
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
    # K = ufloat(medians[f'planet{i+1}_k1'][0], medians[f'planet{i+1}_k1'][1])
    # period = ufloat(median, medians[f'planet{i+1}_period'][1])
    # ecc = ufloat(medians[f'planet{i+1}_ecc'][0], medians[f'planet{i+1}_ecc'][1])
    # mstar = ufloat(0.75, 0.035)
    # print(f"Planet {period:.2f}, m*sin(i) = {min_mass(K, period, ecc, mstar):5f} Mearth;  a = {semi_mayor_axis(period, mstar)} AU")


# plt.tight_layout()
# plt.show()
# -------------------------------------------------------

# ---------- PHASE FOLDING -----------------
# import data
home = os.getenv('HOME')
H03path = os.path.join(
    home, 'tesis/codigo/data/HD40307/HD40307_HARPS03(DRS-3-5).rdb')
H15path = os.path.join(
    home, 'tesis/codigo/data/HD40307/HD40307_HARPS15(DRS-3-5).rdb')

H03data = pd.read_csv(H03path, sep='\t', comment='#', skiprows=[1, ])
H15data = pd.read_csv(H15path, sep='\t', comment='#', skiprows=[1, ])
H03data['vrad'] = H03data['vrad'] - medians['harps_offset'][0]
if 'harps15_offset' in medians.keys():
    H15data['vrad'] = H15data['vrad'] - medians['harps15_offset'][0]
    data = pd.concat([H03data, H15data])  # Join both datasets
    data = data.reset_index(drop=True)  # Reset index
else:
    data = H03data

# Create paramater dictionary with the medians
pardict = medians
for key in pardict.keys():
    pardict[key] = medians[key][0]
pardict.update({'drift1_v0': 0})
pardict.update({'drift1_tref': 54521.6791})
for i in range(nplanets):
    pardict.update({f'planet{i+1}_epoch': 54521.6791})

# Import model to make prediction
import model_rvrd_kepler as mod
time_array = np.linspace(data['rjd'].iloc[0], data['rjd'].iloc[-1], 50000)
prediction = mod.model(pardict, time_array)

# Plot data with model
plt.figure(3)
plt.errorbar(data['rjd'], data['vrad'],
             yerr=data['svrad'], fmt='.r', label='Datos', elinewidth=1)
plt.plot(time_array, prediction, 'k', label='Modelo')
plt.legend()
plt.xlabel("Tiempo (BJD - 2450000)")
plt.ylabel(r"$\Delta$Velocidad radial (m/s)")
# plt.show()


def activity(time):
    tt = (time - pardict['drift1_tref'])/365.25
    if 'drift1_quar' not in pardict.keys():
        quar = 0
    else:
        quar = pardict['drift1_quar']
    return pardict['drift1_unitconstant'] * \
        (pardict['drift1_v0'] + pardict['drift1_lin']*tt + pardict['drift1_qua'] *
         tt**2 + pardict['drift1_cub']*tt**3 + quar*tt**4)


# Corect data for longterm activity
data_activity = activity(data['rjd'])
prediction_activity = activity(time_array)
data_noact = data['vrad'] - data_activity
prediction_noact = prediction - prediction_activity

fig, ax = plt.subplots(1, nplanets, figsize=(13, 5), sharey=True)
for i in range(1, nplanets+1):
    period = pardict[f'planet{i}_period']
    corrected_data = np.array(data_noact)
    for j in range(1, nplanets+1):
        if j == i:
            continue
        model_planet = mod.modelk(
            pardict, np.array(data['rjd']), planet=str(j))
        corrected_data -= model_planet
        # plt.plot(data['rjd'], corrected_data, 'k.')
    prediction = mod.modelk(pardict, np.array(data['rjd']), planet=str(i))

    # Calculate error with jitter
    err03 = np.sqrt(data['svrad'][:239]**2 + pardict['harps_jitter']**2)

    # Calculate reference point with maximum
    t_ref = np.array(data['rjd'])[np.argmax(prediction)]
    phases = foldAt(np.array(data['rjd']), period, T0=t_ref)
    sortIndi = np.argsort(phases)  # Sort by phase

    idx = names[round(period, 0)][1]  # Get index for each planet
    name = names[round(period, 0)][0]  # Get name for each planet
    # Plot data and prediction in corresponding subplot
    try:
        ax[idx].errorbar(phases[:239], corrected_data[:239], yerr=err03,
                         fmt='.', color="xkcd:scarlet", elinewidth=1, barsabove=True)
        err15 = np.sqrt(data['svrad'][239:]**2 + pardict['harps15_jitter']**2)
        ax[idx].errorbar(phases[239:], corrected_data[239:], yerr=err15,
                         fmt='.', color="xkcd:emerald", elinewidth=1, barsabove=True)
        phases = phases[sortIndi]
    except:
        phases = phases[sortIndi]
        corrected_data = corrected_data[sortIndi]
        ax[idx].errorbar(phases[:239], corrected_data[:239], yerr=err03,
                         fmt='r.', elinewidth=1, barsabove=True)
    prediction = prediction[sortIndi]
    ax[idx].plot(phases, prediction, 'b')
    ax[idx].set_title(f"Planeta {name}; Período: {period:.3f} d")
    ax[idx].set_xlabel("Fase orbital")
    if idx == 0:
        ax[idx].set_ylabel("Velocidad radial [m/s]")

plt.show()
# # -------------------------------------------
