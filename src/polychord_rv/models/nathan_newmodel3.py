# New models for the eccentric vs ciruclas planets

import os
from numpy import pi, exp

# Source path of the datafile
datapath = os.path.dirname(os.path.abspath(__file__))

# -------------------
# DATA
# -------------------
datadict = {
    'harps':
    {'datafile': os.path.join(datapath, '../data/corot9b_2.dat'),
     'instrument': 'HARPS03',
     'mask': 'G2',
     'type': 'RV',
     'jittermodel': 'constant',
     'sep': ' '
     },
}

# ---------------------
#        PLANETS
# --------------------
fpdict = {'k1': [0.0, 1, ['Jeffreys', exp(-1), exp(9)]],
          'period': [0.0, 1, ['UniformFrequency', 58.05, 62.07]],
          'ecc': [0.0, 0, ['Beta', 0.867, 3.03]],
          'omega': [0.1, 1, ['Uniform', 0., 2*pi]],
          'ma0': [0.1, 1, ['Uniform', 0., 2*pi]],
          # TODO Define Epoch
          'epoch': [5300, 0]
          }

driftdict = {'v0': [0.0, 0, ['Uniform', -1000., 1000.]],
             'lin': [0.0, 0, ['Normal', 0.035253, 8.82e-03]],
             'qua': [0.0, 0, ['Normal', 0.00246247, 2.542e-03]],
             'cub': [0.0, 0, ['Normal', -0.00123254, 5.879e-04]],
             'tref': [5300, 0],
             'unitconstant': [1.0, 0, ['Uniform', 0., 1000.]],
             }

harpsdict = {'offset': [0., 1, ['Uniform', 19000, 21000]],
             'jitter': [0.75, 1, ['Uniform', 0., 100.]],
             }

input_dict = {'harps': harpsdict,
              'drift1': driftdict,
              'planet1': fpdict}

configdicts = [input_dict, datadict]
