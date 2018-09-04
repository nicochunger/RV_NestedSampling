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
    'rvrd':
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
fp1dict = {'k1': [2.76, 2, ['ModJeffreys', 1., 999.], 0.5],
           'period': [42.0615, 2, ['Jeffreys', 1.0, 1e4], 0.13],
           'ecc': [0.1, 1, ['Uniform', 0., 1.]],
           'omega': [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'ma0': [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'epoch': [278.019, 0]
           }

fp2dict = {'k1': [2.51, 2, ['ModJeffreys', 1., 999.], 0.5],
           'period': [42.1173, 2, ['Jeffreys', 1.0, 1e4], 0.13],
           'secos': [0.1, 1, ['Uniform', -1., 1.]],
           'sesin': [0.1, 1, ['Uniform', -1., 1.]],
           # 'ecc': [0.1, 1, ['Uniform', 0., 1.]],
           # 'omega' : [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'ma0': [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'epoch': [278.019, 0]
           }

fp3dict = {'k1': [2.51, 2, ['ModJeffreys', 1., 999.], 0.5],
           'period': [42.1173, 2, ['Jeffreys', 1.0, 1e4], 0.13],
           'secos': [0.1, 1, ['Uniform', -1., 1.]],
           'sesin': [0.1, 1, ['Uniform', -1., 1.]],
           # 'ecc': [0.1, 1, ['Uniform', 0., 1.]],
           # 'omega' : [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'ma0': [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'epoch': [278.019, 0]
           }

# TODO fijarme que el fit sea sobre el tiempo ajustado (?
driftdict = {'v0': [4.14420729e+03, 0, ['Uniform', -1000., 1000.]],
             'lin': [-2.27489889e-01, 0, ['Uniform', -10., 10.]],
             'qua': [4.15577868e-06, 0],
             'cub': [-2.52946502e-11, 0],
             'tref': [278.019, 0],
             'unitconstant': [1.0, 1, ['Uniform', -1000., 1000.]],
             }

rvrddict = {'offset': [31334.45, 0, ['Uniform', -32000, 32000]],
            'jitter': [0.0, 1, ['ModJeffreys', 1.0, 99.0]],
            }

objectdict = {'rvrd': rvrddict,
              'drift1': driftdict,
              'planet1': fp1dict,
              'planet2': fp2dict,
              'planet3': fp3dict
              }

configdicts = [objectdict, datadict]
