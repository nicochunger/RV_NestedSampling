import numpy as np


def distribution(x):
    """ Distribution to sample """
    mu = 0
    sigma = 1
    return (1/np.sqrt(2*np.pi*sigma**2)) * np.exp((x-mu)**2/(2*sigma**2))

# iter = 1000
# for i in range(iter):
#     # Choose slice


print(distribution(0))
