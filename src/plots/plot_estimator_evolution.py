import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import iqr
from scipy.optimize import curve_fit

# Importing data
path = "/home/nunger/tesis/graficos/nested_sampling_test/"

results = np.loadtxt(path + "results_200-7000.txt")

# Construction of live points array
N = np.loadtxt(path + "N_200-7000.txt")
N = N.astype(int)

# Delete small N values
for i in range(8, -1, -2):
    results = np.delete(results, i, 0)
    N = np.delete(N, i)

results = np.transpose(results) # Transpose of results for boxplot

# ----------------------- BoxPlot --------------------------------------
plt.figure(0, figsize=(10, 6))
plt.boxplot(results, positions=N, widths=200, showfliers=True, whis=[5, 95])
plt.xlabel("Numer of live points", fontsize='13')
plt.ylabel("Evidence", fontsize='13')
plt.title("Evidence distribution for different numer of live points", fontsize='15')
plt.xlim([min(N)-300, max(N)+300])

analytic = 1.131e-38
plt.axhline(y=analytic, color='xkcd:grey', label="Analytic result")
plt.legend()
# -----------------------------------------------------------------------

# ---------------- Standard deviations ---------------------------------
plt.figure(1)
stds = np.std(results, axis=0)
plt.plot(N, stds, 'k.-')
plt.xlabel("Live points", fontsize='13')
plt.ylabel("Standard deviation", fontsize='13')
plt.title("Standard deviation for each numer of live points", fontsize='15')
# ---------------------------------------------------------------------

# -------------- Mean and median evidences ----------------------------
plt.figure(2)
avg = np.mean(results, axis=0)
med = np.median(results, axis=0)
plt.plot(N, avg, 'k.-', label="Mean evidence")
plt.plot(N, med, 'b.-', label="Median evidence")
plt.axhline(y=analytic, color='xkcd:grey', label="Analytic result")
plt.xlabel("Live points")
plt.ylabel("Evidence")
plt.legend()
# ----------------------------------------------------------------------

# ---------------------- L vs Xi --------------------------------------
# N = 200
priormass200 = np.loadtxt(path + "xi_200.txt")
lhoods200 = np.loadtxt(path + "lhoods_200.txt")
# N = 7000
priormass7000 = np.loadtxt(path + "xi_7000.txt")
lhoods7000 = np.loadtxt(path + "lhoods_7000.txt")

plt.figure(3)
plt.plot(priormass200, lhoods200, label="N = 200")
plt.plot(priormass7000, lhoods7000, label="N = 7000")
plt.legend()
plt.xscale('log')
plt.xlabel('Prior mass')
plt.ylabel('Likelihood')
# ---------------------------------------------------------------------

# ------------------ Distribution of N --------------------------------
results500 = np.loadtxt(path + "iter500/results_300_501_200.txt")
n1000 = np.loadtxt(path + "iter500/results_1000_1001_1000.txt")
n1000 = np.reshape(n1000, (1, 500))
results500 = np.append(results500, n1000, axis=0)
N500 = [300, 500, 1000]

# Gaussian model
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

p0 = [10, 1e-38, 1e-38]

# Plot
f, ax = plt.subplots(len(N500), 1, sharex=True)
for i in range(len(N500)):
    hist, bin_edges, patches = ax[i].hist(results500[i], bins='fd', label="Histogram of data")
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2 # Calculate bin centers
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0) # Fit data
    print("Parameters [A, mu, sigma]: " + str(coeff))
    # Standard deviation of each coefficient
    print("Errors for parameters: " + str(np.sqrt(np.diag(var_matrix))))
    print("\n")
    x = np.linspace(min(results500[i]), max(results500[i]), 200)
    hist_fit = gauss(x, *coeff)
    ax[i].axvline(x=1.131e-38, color='xkcd:grey', label="Analytic result")
    ax[i].plot(x, hist_fit, label='Gaussian fit')
    ax[i].set_title("Distribution of results for N = {}".format(N500[i]))
    
plt.xlabel("Evidence", fontsize=14)
ax[1].set_ylabel("Ocurrences", fontsize=14)
ax[0].legend()
# -----------------------------------------------------------------------

plt.show()