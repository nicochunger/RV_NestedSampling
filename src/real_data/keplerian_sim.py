import numpy as np
import matplotlib.pyplot as plt
import model_rvrd_kepler as mod

pardict = {}
pardict.update({'planet1_period': 1})
pardict.update({'planet1_k1': 1})
pardict.update({'planet1_ecc': 0.9})
pardict.update({'planet1_omega': np.pi/4})  # np.pi/2})
pardict.update({'planet1_ma0': np.pi})
pardict.update({'planet1_epoch': 0})

time = np.linspace(0, 1, 1000)

M = time*2*np.pi
true_anom = mod.trueanomaly(M, pardict['planet1_ecc'])

model = mod.modelk(pardict, time, "1")

# plt.plot(time, true_anom)

plt.figure(figsize=(4, 3))
plt.plot(time, model)
plt.xlabel("Fase orbital")
plt.ylabel("Amplitud")
plt.title(
    fr"$e$ = {pardict['planet1_ecc']}; $\omega$ = {pardict['planet1_omega']:.2f}")
plt.show()
