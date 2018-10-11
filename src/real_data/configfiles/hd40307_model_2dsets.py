# Model for the star HD40307
# Number of planets has to be specified when running PolyChord

import os
from numpy import pi

datapath = os.path.join(os.getenv('HOME'), 'codigo/data/HD40307/Night_Binning')

# -------------------
# DATA
# -------------------
datadict = {
    'harps':
    {'datafile':  os.path.join(datapath, 'HD40307_HARPS03(DRS-3-5).rdb'),
     'instrument': 'HARPS03',
     'mask': 'G2',
     'type': 'RV',
     'jittermodel': 'constant'
     },
    'harps09':
    {
        'datafile': os.path.join(datapath, 'HD40307_HARPS(Pub-2009).rdb'),
        'instrument': 'HARPS(Pub2009)'
    }
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
             'lin': [0.0, 1, ['Uniform', -2, 2]],
             'qua': [0.0, 1, ['Uniform', -2, 2]],
             'cub': [0.0, 1, ['Uniform', -2, 2]],
             'tref': [54521.6791, 0],
             'unitconstant': [1.0, 0, ['Uniform', 0., 100.]],
             }

harpsdict = {'offset': [31334.45, 1, ['Uniform', 28996, 33668]],
             'jitter': [0.75, 1, ['Uniform', 0., 50.]],
             'slope': [0.0, 0, ['Uniform', 0., 50.]],
             }

planet_periods = [9.6208, 20.418, 4.311, 51.56, 204.5]
periods_std = [0.0009, 0.005, 0.0003, 0.08, 1]

configdicts = [datadict, fpdict, driftdict, harpsdict]
