import matplotlib.pyplot as plt
from PyAstronomy.pyasl import foldAt
import pandas as pd
import numpy as np
import os
import sys

# Construct medians dictionary
pardict = {'harps_offset': 31334.664,
           'harps_jitter': 0.545,
           'harps_slope': 5.1,
           'harps15_offset': 31349.31,
           'harps15_jitter': 0.7,
           'harps15_slope': 17.4,
           'drift1_v0': 0,
           'drift1_tref': 54521.6791,
           'drift1_lin': 0.038,
           'drift1_qua': 0.00068,
           'drift1_cub': -0.001449,
           'drift1_quar': 8.7e-5,
           'drift1_unitconstant': 24.5
           }

pardict.update({'planet1_period': 4.31149,
                'planet1_k1': 1.859,
                'planet1_ecc': 0.038,
                'planet1_ma0': 0.7,
                'planet1_omega': 2.17
                })

pardict.update({'planet2_period': 9.62193,
                'planet2_k1': 2.294,
                'planet2_ecc': 0.04,
                'planet2_ma0': 3.8,
                'planet2_omega': 1.76
                })

pardict.update({'planet3_period': 20.4180,
                'planet3_k1': 2.49,
                'planet3_ecc': 0.03,
                'planet3_ma0': 1.26,
                'planet3_omega': 1.99
                })

pardict.update({'planet4_period': 0.978386,
                'planet4_k1': 1.006,
                'planet4_ecc': 0.059,
                'planet4_ma0': 4.99,
                'planet4_omega': 3.78
                })

# Identify planet period with name
nplanets = 4
names = {4: ('b', 0), 10: ('c', 1), 20: ('d', 2), 1: ('f', 3), 131: ('f', 3)}

# ---------- PHASE FOLDING -----------------
# import data
home = os.getenv('HOME')
H03path = os.path.join(
    home, 'tesis/codigo/data/HD40307/HD40307_HARPS03(DRS-3-5).rdb')
H15path = os.path.join(
    home, 'tesis/codigo/data/HD40307/HD40307_HARPS15(DRS-3-5).rdb')

H03data = pd.read_csv(H03path, sep='\t', comment='#', skiprows=[1, ])
H15data = pd.read_csv(H15path, sep='\t', comment='#', skiprows=[1, ])
H03data['vrad'] = H03data['vrad'] - pardict['harps_offset']
H15data['vrad'] = H15data['vrad'] - pardict['harps15_offset']
data = pd.concat([H03data, H15data])  # Join both datasets
data = data.reset_index(drop=True)  # Reset index

# Create paramater dictionary with the medians
for i in range(nplanets):
    pardict.update({f'planet{i+1}_epoch': 54521.6791})

# Import model to make prediction
modelpath = os.path.join(home, 'tesis/codigo/src/real_data')
sys.path.append(os.path.abspath(modelpath))
import model_rvrd_kepler as mod

# time_array = np.linspace(data['rjd'].iloc[0], data['rjd'].iloc[-1], 50000)
# prediction = mod.model(pardict, time_array)
# Plot data with model
# plt.figure(3)
# plt.errorbar(data['rjd'], data['vrad'],
#              yerr=data['svrad'], fmt='.r', label='Datos', elinewidth=1)
# plt.plot(time_array, prediction, 'k', label='Modelo')
# plt.legend()
# plt.xlabel("Tiempo (BJD - 2450000)")
# plt.ylabel(r"$\Delta$Velocidad radial (m/s)")
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
# prediction_activity = activity(time_array)
data_noact = data['vrad'] - data_activity
# prediction_noact = prediction - prediction_activity

fig, ax = plt.subplots(1, nplanets, figsize=(13, 4), sharey=True)
for i in range(1, nplanets+1):
    period = pardict[f'planet{i}_period']
    corrected_data = np.array(data_noact)

    # Substract contribution of the other planets
    for j in range(1, nplanets+1):
        if j == i:
            continue
        model_planet = mod.modelk(
            pardict, np.array(data['rjd']), planet=str(j))
        corrected_data -= model_planet
    prediction = mod.modelk(pardict, np.array(data['rjd']), planet=str(i))

    # Calculate reference point with maximum
    t_ref = np.array(data['rjd'])[np.argmax(prediction)]
    phases = foldAt(np.array(data['rjd']), period, T0=t_ref)
    sortIndi = np.argsort(phases)  # Sort by phase

    idx = names[round(period, 0)][1]  # Get index for each planet
    name = names[round(period, 0)][0]  # Get name for each planet
    # Plot data and prediction in corresponding subplot
    # H03 data in scralet red
    err03 = np.sqrt(data['svrad'][:239]**2 + pardict['harps_jitter']**2)
    ax[idx].errorbar(phases[:239], corrected_data[:239], yerr=err03,
                        fmt='.', color='xkcd:scarlet', elinewidth=1, barsabove=True)
    
    # H15 data in emerald green
    err15 = np.sqrt(data['svrad'][239:]**2 + pardict['harps15_jitter']**2)
    ax[idx].errorbar(phases[239:], corrected_data[239:], yerr=err15,
                        fmt='.', color='xkcd:emerald', elinewidth=1, barsabove=True)
    phases = phases[sortIndi]

    prediction = prediction[sortIndi]
    ax[idx].plot(phases, prediction, 'b')
    ax[idx].set_title(f"Planet {name}; Period: {period:.3f} d", fontsize=14)
    ax[idx].set_xlabel("Orbital phase", fontsize=13)
    ax[idx].tick_params(axis='both', which='major', labelsize=13)
    if idx == 0:
        ax[idx].set_ylabel("Radial velocity [m/s]", fontsize=13)

plt.tight_layout()
plt.show()
# # -------------------------------------------