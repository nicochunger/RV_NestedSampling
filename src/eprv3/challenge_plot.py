import subprocess
from pprint import pprint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Import all methods names
dirname = 'paper_results/'
methods = subprocess.check_output(
    'ls -d '+dirname+'*/', shell=True).decode('utf-8').replace(dirname, '').split('/\n')
methods.remove('')
order = [-4, 1, 0, 3, 4, 2, 5, -3, -2, 6, 7, 8, 9, 10, 11, -1]
methods = [methods[i] for i in order]

datasets = ['0001', '0002', '0003', '0004', '0005', '0006']
columns = ['comptime', 'nplanets', 'mode_log10Z',
           'median_log10Z', '-2sigma', '-1sigma', '+1sigma', '+2sigma']


def get_evidences(nplanets, prior):
    """ 
    Recieves the number of planets and the type of prior and generates a pandas
    DataFrame with all the evidences for all methods in all datafiles for that
    particular configuration.
    nplanets (int): Number of planets
    prior (str): Type of prior. 'a' for Broad and 'b' for Narrow
    """
    dframe = pd.DataFrame()
    for dset in datasets:
        evidences = []
        for method in methods:
            try:
                data = pd.read_csv(
                    dirname+method+f'/evidences_{dset}{prior}.txt', names=columns)
                evidences.append(data['median_log10Z'].loc[nplanets])
            except:
                evidences.append(None)
        dframe[dset] = evidences

    dframe.set_index([methods], inplace=True)
    return dframe


# ----------------------- BROAD -----------------------------------

# Generate all DataFrames with the evidences for each planet model
# Broad Prior
broad_Zs = [0] * 4  # Empty list of length 4
for i in range(4):
    Zs = get_evidences(i, 'a')
    med = [Zs.median()[dset] for dset in datasets]
    broad_Zs[i] = (Zs, med)

# Plots for Broad Prior
x = np.arange(1, 7)
offset = np.linspace(-0.45, 0.45, len(methods))
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
ax = ax.flatten()
fig.suptitle('Broad Prior', fontsize=20)
plt.subplots_adjust(top=0.94, bottom=0.05, right=0.98, left=0.07,
                    hspace=0.15, wspace=0.15)
for i in range(4):
    for m, method in enumerate(methods):
        if method == 'PolyChord-nlive25nDims-iter50':
            ax[i].plot(x+offset[m], broad_Zs[i][0].loc[method] -
                       broad_Zs[i][1], 'r.', ms=8, label=method)
        elif method == 'PolyChord-nlive2000-iter6':
            ax[i].plot(x+offset[m], broad_Zs[i][0].loc[method] -
                       broad_Zs[i][1], 'g.', ms=8, label=method)
        elif method == 'multirun-multinest-nlive2000-eff0.3':
            ax[i].plot(x+offset[m], broad_Zs[i][0].loc[method] -
                       broad_Zs[i][1], 'k.', ms=8, label=method)
        else:
            ax[i].plot(x+offset[m], broad_Zs[i][0].loc[method] -
                       broad_Zs[i][1], '.', color='xkcd:sky blue', label='_nolegend_')
    # Horizontal line at 0
    ax[i].hlines(0, 0.5, 6.5, linestyles='--', linewidth=1)
    # Horizontal lines to separate datsets
    ax[i].vlines(np.arange(0.5, 6, 1), -5, 5,
                 linestyles='-', linewidth=1)
    props = dict(boxstyle='square', facecolor='white',
                 alpha=1, linewidth=1)
    ax[i].text(0.08, 0.08, f'{i}-planet model', transform=ax[i].transAxes,
               fontsize=17, verticalalignment='bottom', bbox=props)
    ax[i].set_xlim([0.5, 6.5])
    ax[i].set_ylim([-5, 5])
    ax[i].set_xlabel('Dataset Number')
    ax[i].set_ylabel(r'log(Z) - <log(Z)>')
    # Put legend only in first subplot
    if i == 0:
        ax[i].legend(loc=1, fontsize=12)
# -------------------------------------------------------------

# # --------------------------- NARROW --------------------------

# # Generate all DataFrames with the evidences for each planet model
# # Narrow Prior
# narrow_Zs = [0] * 4  # Empty list of length 4
# for i in range(4):
#     Zs = get_evidences(i, 'b')
#     med = [Zs.median()[dset] for dset in datasets]
#     narrow_Zs[i] = (Zs, med)

# # Plots for Narrow Prior
# x = np.arange(1, 7)
# offset = np.linspace(-0.45, 0.45, len(methods))
# fig, ax = plt.subplots(2, 2, figsize=(12, 8))
# ax = ax.flatten()
# fig.suptitle('Narrow Prior', fontsize=20)
# plt.subplots_adjust(top = 0.94, bottom = 0.05, right = 0.98, left = 0.07,
#             hspace = 0.15, wspace = 0.15)
# for i in range(4):
#     for m, method in enumerate(methods):
#         if method == 'PolyChord':
#             ax[i].plot(x+offset[m], narrow_Zs[i][0].loc[method] -
#                        narrow_Zs[i][1], 'r.', label=method)
#         else:
#             ax[i].plot(x+offset[m], narrow_Zs[i][0].loc[method] -
#                        narrow_Zs[i][1], 'b.', label=method)
#     # Horizontal line at 0
#     ax[i].hlines(0, 0.5, 6.5, linestyles='--', linewidth=1)
#     # Horizontal lines to separate datsets
#     ax[i].vlines(np.arange(0.5, 6, 1), -5, 5,
#                  linestyles='-', linewidth=1)
#     props = dict(boxstyle='square', facecolor='white',
#                  alpha=1, linewidth=1)
#     ax[i].text(0.08, 0.08, f'{i}-planet model', transform=ax[i].transAxes,
#                fontsize=17, verticalalignment='bottom', bbox=props)
#     ax[i].set_xlim([0.5, 6.5])
#     ax[i].set_ylim([-5, 5])
#     ax[i].set_xlabel('Dataset Number')
#     ax[i].set_ylabel(r'log(Z) - <log(Z)>')
#     if i == 0:
#         ax[i].legend(loc=4, fontsize=5)
# # ---------------------------------------------------------------

plt.show()
