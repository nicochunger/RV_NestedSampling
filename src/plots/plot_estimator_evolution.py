import numpy as np 
import matplotlib.pyplot as plt

path = "../../../graficos/nested_sampling_test/"

results = np.loadtxt(path + "results_200_2001_200.txt")
results = np.append(results, np.loadtxt(path + "results_2500_3001_500.txt"), axis=0)
results = np.transpose(results)

N = np.arange(200, 2001, 200)
N = np.append(N, [2500, 3000])

plt.figure(0, figsize=(8, 6))
plt.boxplot(results, labels=N, showfliers=False)
plt.xlabel("Numer of live points")
plt.ylabel("Evidence")
plt.title("Evidence distribution for different numer of live points")

analytic = 1.131e-38
plt.axhline(y=analytic, color='xkcd:grey', label="Analytic result")
plt.legend()

plt.show()