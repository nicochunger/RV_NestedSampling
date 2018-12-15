# Program to calculate the new drift for the new data

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

home = os.getenv('HOME')
datapath = os.path.join(home, 'tesis/codigo/data/HD40307')

# Import data for both instruments
harps03 = pd.read_csv(
    datapath+'/HD40307_HARPS03(DRS-3-5).rdb', skiprows=[1, ], sep='\t')

harps15 = pd.read_csv(
    datapath+'/HD40307_HARPS15(DRS-3-5).rdb', skiprows=[1, ], sep='\t')

# ---------- Polynomial Fit -------------------------
# Fit 3rd degree polynomial to the combined activity
# Adjust observation time with reference time
harps03['rjd'] = (harps03['rjd'] - 54521.6791)/365.25
harps15['rjd'] = (harps15['rjd'] - 54521.6791)/365.25
rjd = harps03['rjd'].append(harps15['rjd'])
rhk = harps03['rhk'].append(harps15['rhk'])

fit_pars, cov = np.polyfit(rjd, rhk, 4, cov=True)
errs = np.sqrt(np.diag(cov))
print(fit_pars)
print(errs)
x = np.linspace(harps03['rjd'][0], harps15['rjd'].iloc[-1], 500)
fit = np.polyval(fit_pars, x)

plt.plot(harps03['rjd'], harps03['rhk'], 'k.')
plt.plot(harps15['rjd'], harps15['rhk'], 'y.')
plt.plot(x, fit, 'r')
plt.show()
# -------------------------------------------------


# ------------------ Keplerian Fit -----------------------------
# Because with the new data a 3rd degree polynomial is no longer enough
# We try with a keplerian curve


# ------------------------------------------------------------
