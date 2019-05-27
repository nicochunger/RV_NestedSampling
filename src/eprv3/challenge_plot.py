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
    error_frame = pd.DataFrame()
    for dset in datasets:
        evidences = []
        errors = []
        for method in methods:
            try:
                data = pd.read_csv(
                    dirname+method+f'/evidences_{dset}{prior}.txt', names=columns)
                evidences.append(data['median_log10Z'].loc[nplanets])
                sigma = (data['+1sigma'].loc[nplanets] -
                         data['-1sigma'].loc[nplanets])/2
                errors.append(sigma)
            except:
                evidences.append(None)
                errors.append(None)
        dframe[dset] = evidences
        error_frame[dset] = errors

    dframe.set_index([methods], inplace=True)
    error_frame.set_index([methods], inplace=True)
    return dframe, error_frame


def plot_results(prior):
    """
    Plots all the results from all the methods for the specified type of prior.
    """
    # Generate all DataFrames with the evidences for each planet model
    Zs_med = [0] * 4  # Empty list of length 4
    stds = [0] * 4
    for i in range(4):
        Zs, Zerrs = get_evidences(i, prior)
        med = [Zs.median()[dset] for dset in datasets]
        Zs_med[i] = (Zs, med, Zerrs)
        std = np.array([Zs.std()[dset] for dset in datasets])
        stds[i] = np.mean(std)

    print(f'Desviacion estandar para cada modelo: {stds}')
    mean_diff = []
    # print(Zs_med[0][0])
    # print(Zs_med[0][1])

    # Plots
    x = np.arange(1, 7)
    offset = np.linspace(-0.45, 0.45, len(methods))
    fig, ax = plt.subplots(4, 1, figsize=(7, 10))
    ax = ax.flatten()
    # if prior == 'a':
    #     fig.suptitle('Broad Prior', fontsize=20)
    # elif prior == 'b':
    #     fig.suptitle('Narrow Prior', fontsize=20)
    plt.subplots_adjust(top=0.94, bottom=0.05, right=0.98, left=0.07,
                        hspace=0.15, wspace=0.15)
    for i in range(4):
        # print(f'STD model {i}: {10**Zs_med[i][0].std()["0005"]}')
        diff = np.zeros([3, 6])
        std_method = np.zeros([3, 6])
        for m, method in enumerate(methods):
            if method == 'PolyChord-nlive25nDims-iter50':
                ax[i].errorbar(x+offset[m], Zs_med[i][0].loc[method] -
                               Zs_med[i][1], yerr=Zs_med[i][2].loc[method], fmt='X', color='xkcd:reddish orange', ms=6, label=method)
                diff[0, :] = np.abs(
                    Zs_med[i][0].loc[method].values - Zs_med[i][1])
                std_method[0, :] = Zs_med[i][2].loc[method]
            elif method == 'PolyChord-nlive2000-iter6':
                ax[i].errorbar(x+offset[m], Zs_med[i][0].loc[method] -
                               Zs_med[i][1], yerr=Zs_med[i][2].loc[method], fmt='b*', color='xkcd:marine blue', ms=6, label=method)
                diff[1, :] = np.abs(
                    Zs_med[i][0].loc[method].values - Zs_med[i][1])
                std_method[1, :] = Zs_med[i][2].loc[method]
            elif method == 'multirun-multinest-nlive2000-eff0.3':
                ax[i].errorbar(x+offset[m], Zs_med[i][0].loc[method] -
                               Zs_med[i][1], yerr=Zs_med[i][2].loc[method], fmt='kD', ms=4, label=method)
                diff[2, :] = np.abs(
                    Zs_med[i][0].loc[method].values - Zs_med[i][1])
                std_method[2, :] = Zs_med[i][2].loc[method]
            else:
                ax[i].errorbar(x+offset[m], Zs_med[i][0].loc[method] -
                               Zs_med[i][1], yerr=Zs_med[i][2].loc[method], fmt='.', color='xkcd:steel', label='_nolegend_')
        # print(diff)
        diff_norm = diff/std_method
        min_max = np.array(
            [np.min(diff_norm, axis=1), np.max(diff_norm, axis=1)])
        print(f"Modelo {i} planetas:\n")
        print(np.transpose(min_max))
        # print(np.transpose(min_max)/std_method)
        # print(std_method)

        print("\n")

        # Horizontal line at 0
        ax[i].hlines(0, 0.5, 6.5, linestyles='--', linewidth=1)
        # Horizontal lines to separate datsets
        ax[i].vlines(np.arange(0.5, 6, 1), -5, 5,
                     linestyles='-', linewidth=1)
        props = dict(boxstyle='square', facecolor='white',
                     alpha=1, linewidth=1)

        for dset in range(6):
            logZ_str = r"$\left<\log(\mathcal{{Z}})\right > =$"
            ax[i].annotate(logZ_str + f"\n{Zs_med[i][1][dset]:.2f}", xy=(
                dset+1, 4), horizontalalignment="center", verticalalignment="center")

        if i == 1:
            ax[i].text(0.08, 0.08, f'{i} planeta', transform=ax[i].transAxes,
                       fontsize=17, verticalalignment='bottom', bbox=props)
        else:
            ax[i].text(0.08, 0.08, f'{i} planetas', transform=ax[i].transAxes,
                       fontsize=15, verticalalignment='bottom', bbox=props)
        ax[i].set_xlim([0.5, 6.5])
        ax[i].set_ylim([-5, 5])
        if i == 3:
            ax[i].set_xlabel('Set de datos')
        ax[i].set_ylabel(
            r'$\log(\mathcal{Z}) - \left< \log(\mathcal{Z}) \right>$')
        # Put legend only in first subplot
        if i == 0:
            ax[i].legend(loc=4, fontsize=10)
    # print(
    #     f'Mean difference of polychord with mean: {np.transpose(np.array(mean_diff))}')
    # print(np.transpose(np.array(mean_diff))/stds)
# ------------------------------------------------------------------


# Plot Borad Priors
plot_results('a')

# Plot Narrow Priors
# plot_results('b')

# Show plots
# plt.show()
