# Number of planets has to be specified when running PolyChord

import os
from numpy import pi

# Source path of the datafile
datapath = os.path.dirname(os.path.abspath(__file__))

# -------------------
# DATA
# -------------------
datadict = {
    'harps':
    {'datafile': os.path.join(datapath, 'data/corot9b_2.dat'),  # Datafile path+name
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
fpdict = {'k1': [0.0, 1, ['Uniform', 0., 40.]],
          'period': [0.0, 1, ['UniformFrequency', 58.05, 62.07]],
          'ecc': [0.1, 1, ['Uniform', 0., 1.]],
          'omega': [0.1, 1, ['Uniform', 0., 2*pi]],
          'ma0': [0.1, 1, ['Uniform', 0., 2*pi]],
          'epoch': [5376, 0]
          }

driftdict = {'v0': [0.0, 0, ['Uniform', -1000., 1000.]],
             'lin': [0.0, 0, ['Normal', 0.035253, 8.82e-03]],
             'qua': [0.0, 0, ['Normal', 0.00246247, 2.542e-03]],
             'cub': [0.0, 0, ['Normal', -0.00123254, 5.879e-04]],
             'tref': [5376, 0],
             'unitconstant': [1.0, 0, ['Uniform', 0., 1000.]],
             }

# TODO Fijarse tema de velocidades si estan en km/s
harpsdict = {'offset': [0., 1, ['Uniform', 19000, 21000]],
             'jitter': [0.75, 1, ['Uniform', 0., 100.]],
             'slope': [0.0, 0, ['Uniform', 0., 50.]],
             }

input_dict = {'harps': harpsdict,
              'drift1': driftdict,
              'planet1': fpdict}

configdicts = [input_dict, datadict]
