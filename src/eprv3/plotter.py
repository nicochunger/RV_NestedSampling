import numpy as np 
import matplotlib.pyplot as plt 

# Plot of RV data
data = np.loadtxt('../../data/rvs_0001.txt')

time = data[:,0]
rv = data[:,1]
error = data[:,2]

plt.errorbar(time, rv, error, fmt='k.', ecolor='xkcd:grey')
plt.hlines(0, -25, 620, linestyles='dashed', colors='xkcd:light grey')
plt.xlim([-25, 620])
plt.xlabel(r'Time $(days)$', fontsize='12')
plt.ylabel(r'Radial velocity $(m/s)$', fontsize='12')
plt.show()