# Import all the data

import pandas as pd
import numpy as np

k1 = pd.read_csv("results_1a.txt", sep='\t')  # , skiprows=[1, ])
k2 = pd.read_csv("results_2a.txt", sep='\t')  # , skiprows=[1, ])
k0 = pd.read_csv("results_0a.txt", sep='\t')  # , skiprows=[1, ])
k3 = pd.read_csv("results_3a.txt", sep='\t')  # , skiprows=[1, ])
k4 = pd.read_csv("results_4a.txt", sep='\t')  # , skiprows=[1, ])
k5 = pd.read_csv("results_5a.txt", sep='\t', comment='#')  # , skiprows=[1, ])

ks = [k0, k1, k2, k3, k4, k5]

print("NEW DATA!!")

for i, k in enumerate(ks):
    med = k['logZ'].median()
    std = k['logZ'].std()
    print(f"Evidence for model k{i} = {med:.2f} +/- {std:.2f}")

print("\n")
# Odds Ratios

odds = []
for i in range(len(ks)-1):
    med1 = ks[i]['logZ'].median()
    std1 = ks[i]['logZ'].std()
    med2 = ks[i+1]['logZ'].median()
    std2 = ks[i+1]['logZ'].std()
    odds = np.log10(np.exp(med2 - med1))
    odds_err = np.sqrt(std1**2 + std2**2) / np.log(10)
    print(f"Odds ratio for {i+1}v{i} = {odds:.2f} +/- {odds_err:.2f}")
