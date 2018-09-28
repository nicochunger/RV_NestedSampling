# Model for the star HD40307
# Number of planets has to be specified when running PolyChord

# import os
from numpy import pi

# datapath = os.path.join(os.getenv('HOME'), 'tesis/codigo/data/')

# -------------------
# DATA
# -------------------
datadict = {
    'harps':
    {'datafile': 'data/data_a1.txt',
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
          'epoch': [218.79, 0]
          }

driftdict = {'v0': [0.0, 0, ['Uniform', -1000., 1000.]],
             'lin': [0.0, 0, ['Normal', 0.035253, 8.82e-03]],
             'qua': [0.0, 0, ['Normal', 0.00246247, 2.542e-03]],
             'cub': [0.0, 0, ['Normal', -0.00123254, 5.879e-04]],
             'tref': [218.79, 0],
             'unitconstant': [1.0, 0, ['Uniform', 0., 1000.]],
             }

harpsdict = {'offset': [0., 1, ['Uniform', -10, 10]],
             'jitter': [0.75, 1, ['Uniform', 0., 100.]],
             'slope': [0.0, 0, ['Uniform', 0., 50.]],
             }

input_dict = {'harps': harpsdict,
              'drift1': driftdict,
              'planet1': fpdict}

configdicts = [input_dict, datadict]
