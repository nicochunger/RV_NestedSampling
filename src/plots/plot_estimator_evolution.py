import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import iqr

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

# Use Freedman Draconis rule to choose bin size
IQR = iqr(results500[0]) # Interquartile range
binSize = 2 * (IQR / len(results500[0])**(1./3))
nBins = int((np.max(results500[0]) - np.min(results500[0])) / binSize)
# Plot
plt.figure(2)
plt.hist(results500[0], nBins)
plt.xlabel("Evidence")
plt.ylabel("Ocurrences")

plt.show()