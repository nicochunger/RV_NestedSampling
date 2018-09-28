import pickle
import matplotlib.pyplot as plt
from PyAstronomy.pyasl import foldAt
import pandas as pd
import numpy as np
from pprint import pprint
from scipy.optimize import curve_fit
import argparse

# Read arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fit', action='store_true',
                    help='If posterior should be fit with a gaussian.')
args_params = parser.parse_args()

datapath = '/home/nunger/tesis/codigo/data/HD40307/'
data = pd.read_csv(datapath+'HD40307_HARPS03(DRS-3-5)_night_binning.rdb',
                   sep='\t', comment='#', skiprows=[1, ])

t = data['rjd']
y = data['vrad']
# Correct offset
y = y - 31334.45
ey = data['svrad']
# rhk = data['rhk']

# # Plot of raw RV data
# plt.figure(0)
# plt.errorbar(t, y, yerr=ey, fmt='k.')

# Load evidence results
output = pickle.load(open('output.p', 'rb'))
print(output.logZ)
posterior = output.posterior
samples = posterior.samples

# Get the medians for each parameter
paramnames = posterior.getParamNames().list()
medians = dict(zip(paramnames, zip(np.round(np.median(
    samples, axis=0), 5), np.round(np.std(samples, axis=0), 7))))
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

# # Gaussian model
# def gauss(x, *p):
#     A, mu, sigma = p
#     return A*np.exp(-(x-mu)**2/(2.*sigma**2))

fig2, ax2 = plt.subplots(4, 1)
ax2[0].hist(post['drift1_lin']*1e4,
            bins='auto', histtype='step', normed=True)
ax2[1].hist(post['drift1_qua']*1e4,
            bins='auto', histtype='step', normed=True)
ax2[2].hist(post['drift1_cub']*1e4,
            bins='auto', histtype='step', normed=True)
ax2[3].hist(post['drift1_unitconstant'],
            bins='auto', histtype='step', normed=True)
plt.show()

# # Plot posterior for the period of each planet
# fig, ax = plt.subplots(nplanets, 1)
# ax = np.atleast_1d(ax)  # To support 1 planet models

# for i in range(nplanets):
#     # Histogram
#     hist, bin_edges, patches = ax[i].hist(
#         post[f'planet{i+1}_period'], label='Posterior', bins='auto', histtype='step', normed=True)

#     # Time series of posterior
#     # ax[i].plot(post[f'planet{i+1}_period'], '.')

#     ax[i].set_title(
#         f"Planet {names[i]}. P = {medians[f'planet{i+1}_period'][0]:5.2f}")
#     ax[i].set_xlabel('Period [d]')
#     ax[i].set_ylabel('PDF')

#     if args_params.fit:
#         # Gaussian fit
#         # Calculate bin centers
#         bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
#         p0 = [10, medians[f'planet{i+1}_period']
#               [0], medians[f'planet{i+1}_period'][1]]
#         coeff, var_matrix = curve_fit(
#             gauss, bin_centres, hist, p0=p0)  # Fit data
#         x = np.linspace(min(post[f'planet{i+1}_period']),
#                         max(post[f'planet{i+1}_period']), 2000)
#         hist_fit = gauss(x, *coeff)
#         ax[i].plot(x, hist_fit, 'r-', label='Gaussian fit')

#         # Information box of fit parameters
#         textstr = '\n'.join((
#             'Fitted values:',
#             r'$A=%.5f$' % (coeff[0], ),
#             r'$\mu=%.5f$' % (coeff[1], ),
#             r'$\sigma=%.5f$' % (abs(coeff[2]), )))
#         # these are matplotlib.patch.Patch properties
#         props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#         # place a text box in upper left in axes coords
#         ax[i].text(0.05, 0.95, textstr, transform=ax[i].transAxes,
#                    fontsize=11, verticalalignment='top', bbox=props)
#         ax[i].legend()
# plt.show()
# # -------------------------------------------------------

# # ---------- PHASE FOLDING -----------------
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
