target = 'HD40307'
runid = 'k1'

import os
import pickle
import numpy as np

datapath = os.path.join(os.getenv('HOME'), 'tesis/codigo/data/')

# -------------------
# DATA
# -------------------
datadict = {
    'harps':
    {'datafile':  os.path.join(datapath, target+'/HD40307_HARPS03(DRS-3-5)_night_binning.rdb'),
     'instrument': 'HARPS03',
     'mask': 'G2',
     'type': 'RV',
     'jittermodel': 'constant'
     },
}

# ---------------------
#        PLANETS
# --------------------
fp1dict = {'k1': [2.76, 2, ['Uniform', 0., 200.]],
           'period': [42.0615, 2, ['Jeffreys', 1.0, 1e4]],
           'ecc': [0.1, 1, ['Uniform', 0., 1.]],
           'omega': [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'ma0': [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'epoch': [55500, 0]
           }

driftdict = {'v0': [0.0, 0, ['Uniform', -1000., 1000.]],
             'lin': [0.0, 1, ['Uniform', -10., 10.]],
             'qua': [0.0, 1, ['Uniform', -10., 10.]],
             'cub': [0.0, 1, ['Uniform', -10., 10.]],
             'tref': [55500, 0],
             'unitconstant': [1.0, 0, ['Uniform', -100., 100.]],
             }

harpsdict = {'offset': [31334.45, 1, ['Uniform', 31000, 32000]],
             'jitter': [0.75, 0, ['ModJeffreys', 1.0, 99.0]],
             }

objectdict = {'harps': harpsdict,
              'drift1': driftdict,
              'planet1': fp1dict,
              }

configdicts = [objectdict, datadict]
