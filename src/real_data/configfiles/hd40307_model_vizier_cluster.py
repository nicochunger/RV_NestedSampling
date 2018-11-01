# Model for the star HD40307
# Number of planets has to be specified when running PolyChord

import os
from numpy import pi

home = os.getenv('HOME')
if 'spectro' in home:
    datapath = os.path.join(home, 'codigo/data/')
else:
    datapath = os.path.join(home, 'tesis/codigo/data/')


# -------------------
# DATA
# -------------------
datadict = {
    'harps':
    {'datafile':  os.path.join(datapath, 'HD40307/HD40307_vizier.txt'),
     'instrument': 'HARPS03',
     'mask': 'G2',
     'type': 'RV',
     'jittermodel': 'constant'
     },
}

# ---------------------
#        PLANETS
# --------------------
fpdict = {'k1': [0.0, 1, ['Uniform', 0., 200.]],
          'period': [0.0, 1, ['Jeffreys', 1.0, 1e4]],
          'ecc': [0.1, 1, ['Beta', 0.867, 3.03]],
          'omega': [0.1, 1, ['Uniform', 0., 2*pi]],
          'ma0': [0.1, 1, ['Uniform', 0., 2*pi]],
          'epoch': [54521.6791, 0]
          }

driftdict = {'v0': [0.0, 0, ['Uniform', -1000., 1000.]],
             'lin': [0.0, 1, ['Normal', 0.03738, 8.6e-04]],
             'qua': [0.0, 1, ['Normal', 0.00463, 4.6e-04]],
             'cub': [0.0, 1, ['Normal', -0.00173, 1.0e-04]],
             'tref': [54521.6791, 0],
             'unitconstant': [1.0, 1, ['Uniform', 0., 100.]],
             }

harpsdict = {'offset': [31334.45, 1, ['Uniform', 28996, 33668]],
             'jitter': [0.75, 1, ['Uniform', 0., 50.]],
             'slope': [0.0, 1, ['Uniform', 0., 50.]],
             }

input_dict = {'harps': harpsdict,
              'drift1': driftdict,
              }

planet_periods = [9.6208, 20.418, 4.311, 51.56, 204.5]
periods_std = [0.0009, 0.005, 0.0003, 0.08, 1]

configdicts = [datadict, fpdict, input_dict]
