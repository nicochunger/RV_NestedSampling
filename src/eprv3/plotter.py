import numpy as np 
import matplotlib.pyplot as plt 

# # Plot of RV data
# data = np.loadtxt('../../data/rvs_0001.txt')

# time = data[:,0]
# rv = data[:,1]
# error = data[:,2]

# plt.errorbar(time, rv, error, fmt='k.', ecolor='xkcd:grey')
# plt.hlines(0, -25, 620, linestyles='dashed', colors='xkcd:light grey')
# plt.xlim([-25, 620])
# plt.xlabel(r'Time $(days)$', fontsize='12')
# plt.ylabel(r'Radial velocity $(m/s)$', fontsize='12')
# plt.show()

d10 = np.loadtxt('results0001_0.txt')[:,3]
d11 = np.loadtxt('results0001_1.txt')[:,3]

d20 = np.loadtxt('results0002_0.txt')[:,3]
d21 = np.loadtxt('results0002_1.txt')[:,3]

d30 = np.loadtxt('results0003_0.txt')[:,3]
d31 = np.loadtxt('results0003_1.txt')[:,3]

meanZ = np.zeros((2,3))
meanZ[0,0] = np.mean(d10)
meanZ[0,1] = np.mean(d20)
meanZ[0,2] = np.mean(d30)
meanZ[1,0] = np.mean(d11)
meanZ[1,1] = np.mean(d21)
meanZ[1,2] = np.mean(d31)

print(meanZ)
