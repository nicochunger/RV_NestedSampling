import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import iqr
from scipy.optimize import curve_fit

# Importing data
path = "/home/nunger/tesis/graficos/nested_sampling_test/"

results = np.loadtxt(path + "results_200_2001_200.txt")
results = np.append(results, np.loadtxt(path + "results_2500_3001_500.txt"), axis=0)
results = np.append(results, np.loadtxt(path + "results_4000_5001_1000.txt"), axis=0)
results = np.transpose(results)

# Construction of live points array
N = np.arange(200, 2001, 200)
N = np.append(N, [2500, 3000, 4000, 5000])

# # BoxPlot
# plt.figure(0, figsize=(15, 8))
# plt.boxplot(results, positions=N, widths=150, showfliers=True, whis=[5, 95])
# plt.xlabel("Numer of live points")
# plt.ylabel("Evidence")
# plt.title("Evidence distribution for different numer of live points")
# plt.xlim([-100, 5100])

# analytic = 1.131e-38
# plt.axhline(y=analytic, color='xkcd:grey', label="Analytic result")
# plt.legend()

# # Standard deviations
# plt.figure(1)
# stds = np.std(results, axis=0)
# plt.plot(N, stds, 'k.-')
# plt.xlabel("Live points")
# plt.ylabel("Standard deviation")

# Distribution of N
results500 = np.loadtxt(path + "iter500/results_300_501_200.txt")
n1000 = np.loadtxt(path + "iter500/results_1000_1001_1000.txt")
n1000 = np.reshape(n1000, (1, 500))
results500 = np.append(results500, n1000, axis=0)
N500 = [300, 500, 1000]
# Gaussian model
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

p0 = [4, 1.0e-38, 1.0e-38]

# Plot
f, ax = plt.subplots(len(N500), 1, sharex=True)
for i in range(len(N500)):
    hist, bin_edges, patches = ax[i].hist(results500[i], bins='fd', label="Histogram of data") #nBins[i])
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2 # Calculate bin centers
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0) # Fit data
    print("Parameters [A, mu, sigma]: " + str(coeff))
    print("Errors for parameters: " + str(np.sqrt(np.diag(var_matrix)))) # Standard deviation of each coefficient
    print("\n")
    x = np.linspace(min(results500[i]), max(results500[i]), 200)
    hist_fit = gauss(x, *coeff)
    ax[i].axvline(x=1.131e-38, color='xkcd:grey', label="Analytic result")
    ax[i].plot(x, hist_fit, label='Gaussian fit')
    ax[i].set_title("Distribution of results for N = {}".format(N500[i]))
    
plt.xlabel("Evidence", fontsize=14)
ax[1].set_ylabel("Ocurrences", fontsize=14)
ax[0].legend()
plt.show()