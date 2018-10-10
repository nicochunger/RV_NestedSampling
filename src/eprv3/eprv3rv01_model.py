# Model for the EPRV3 Challenge

import os
from numpy import pi

#datapath = os.path.join(os.getenv('HOME'), 'tesis/codigo/data')

# -------------------
# DATA
# -------------------
datadict = {
    'eprv':
    {'datafile':  '',  # os.path.join(datapath, 'rvs_0001.txt')
     'instrument': 'HARPS',
     'mask': 'G2',
     'type': 'RV',
     'jittermodel': 'constant'
     },
}

# ---------------------
#        PLANETS
# --------------------
fp1dict = {'k1': [0, 1, ['ModJeffreys', 1., 999.], 0.5],
           'period': [0, 1, ['Jeffreys', 1.25, 1e4], 0.13],
           'ecc': [0.1, 1, ['TruncatedRayleigh', 0.2, 1.]],
           'omega': [0.1, 1, ['Uniform', 0., 2*pi]],
           'ma0': [0.1, 1, ['Uniform', 0., 2*pi]],
           'epoch': [278.019, 0]
           }

driftdict = {'v0': [0.0, 1, ['Uniform', -1000., 1000.]],
             'lin': [0.0, 0, ['Uniform', -10., 10.]],
             'qua': [0.0, 0],
             'cub': [0.0, 0],
             'tref': [278.019, 0],
             'unitconstant': [1.0, 0],
             }

eprvdict = {'offset': [0.0, 0, ['Uniform', -10, 10]],
            'jitter': [0.0, 1, ['ModJeffreys', 1.0, 99.0]],
            }

narrow_priors = {
    1: [(39.8107, 44.6684), (11.4815, 12.8825), (10.0, 10, 7452)],
    2: [(15.4882, 16.2181), (14.7911, 17.0608), (158.489, 251.189)],
    3: [(81.2831, 107.152), (38.0189, 42.658), (16.5959, 17.5792)],
    4: [(138.038, 204.174), (15.1356, 16.5959), (398.107, 1000.0)],
    5: [(29.5121, 32.3594), (10.7152, 11.4815), (18.197, 19.9526)],
    6: [(79.4328, 141.254), (31.6228, 50.1187), (316.228, 398.107)]
}

configdicts = [datadict, fp1dict, driftdict, eprvdict]
