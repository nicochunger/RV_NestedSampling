# Model for the star HD40307
# Number of planets has to be specified when running PolyChord

import os
from numpy import pi

datapath = os.path.join(os.getenv('HOME'), 'tesis/codigo/data/')

# -------------------
# DATA
# -------------------
datadict = {
    'harps':
    {'datafile':  os.path.join(datapath, 'HD40307/HD40307_HARPS03(DRS-3-5)_night_binning.rdb'),
     'instrument': 'HARPS03',
     'mask': 'G2',
     'type': 'RV',
     'jittermodel': 'constant'
     },
}

# ---------------------
#        PLANETS
# --------------------
fpdict = {'k1': [0.0, 1, ['Uniform', 0., 200.], 0.5],
          'period': [0.0, 1, ['Jeffreys', 1.0, 1e4], 0.13],
          'ecc': [0.1, 1, ['Beta', 0.867, 3.03]],
          'omega': [0.1, 1, ['Uniform', 0., 2*pi]],
          'ma0': [0.1, 1, ['Uniform', 0., 2*pi]],
          'epoch': [54521.6791, 0]
          }

driftdict = {'v0': [0.0, 0, ['Uniform', -1000., 1000.]],
             'lin': [0.0, 1, ['Normal', 0.035253, 8.82e-03]],
             'qua': [0.0, 1, ['Normal', 0.00246247, 2.542e-03]],
             'cub': [0.0, 1, ['Normal', -0.00123254, 5.879e-04]],
             'tref': [54521.6791, 0],
             'unitconstant': [1.0, 1, ['Uniform', 0., 1000.]],
             }

harpsdict = {'offset': [31334.45, 1, ['Uniform', 28996, 33668]],
             'jitter': [0.75, 1, ['ModJeffreys', 1.0, 99.0]],
             'slope': [0.0, 1, ['Uniform', -200., 200.]],
             }

configdicts = [datadict, fpdict, driftdict, harpsdict]
