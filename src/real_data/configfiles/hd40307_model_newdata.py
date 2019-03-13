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
    {
        'datafile':  os.path.join(datapath, 'HD40307/HD40307_HARPS03(DRS-3-5).rdb'),
        'instrument': 'HARPS03',
        'mask': 'G2',
        'type': 'RV',
        'jittermodel': 'constant'
    },
    'harps15':
    {
        'datafile': os.path.join(datapath, 'HD40307/HD40307_HARPS15(DRS-3-5).rdb'),
        'instrument': 'HARPS15'
    }
}

# ---------------------
#        PLANETS
# --------------------
fpdict = {'k1': [0.0, 1, ['Uniform', 0., 200.]],
          'period': [0.0, 1, ['UniformFrequency', 0.5, 1e4]],
          'ecc': [0.1, 1, ['Beta', 0.867, 3.03]],
          'omega': [0.1, 1, ['Uniform', 0., 2*pi]],
          'ma0': [0.1, 1, ['Uniform', 0., 2*pi]],
          'epoch': [54521.6791, 0]
          }

driftdict = {'v0': [0.0, 0, ['Uniform', -1000., 1000.]],
             'lin': [0.0, 1, ['Normal', 3.755e-02, 1.04e-03]],
             'qua': [0.0, 1, ['Normal', 4.914e-04, 2.51e-04]],
             'cub': [0.0, 1, ['Normal', -1.477e-03, 7.84e-05]],
             'quar': [0.0, 1, ['Normal', 9.715e-05, 8.16e-06]],
             'tref': [54521.6791, 0],
             'unitconstant': [1.0, 1, ['Uniform', 0., 100.]],
             }

harpsdict = {'offset': [31334.45, 1, ['Uniform', 28996, 33668]],
             'jitter': [0.75, 1, ['Uniform', 0., 50.]],
             'slope': [0.0, 1, ['Uniform', 0., 50.]],
             }

harps15dict = {'offset': [31334.45, 1, ['Uniform', 28996, 33668]],
               'jitter': [0.75, 1, ['Uniform', 0., 50.]],
               'slope': [0.0, 1, ['Uniform', 0., 50.]],
               }

input_dict = {'harps': harpsdict,
              'harps15': harps15dict,
              'drift1': driftdict,
              }

planet_periods = [9.6208, 20.418, 4.311, 51.56, 204.5]
periods_std = [0.0009, 0.005, 0.0003, 0.08, 1]

configdicts = [datadict, fpdict, input_dict]
