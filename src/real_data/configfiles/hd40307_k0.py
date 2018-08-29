target = 'HD40307'
runid = 'k1'

import os
import pickle
import numpy as np

datapath = os.path.join(os.getenv('HOME'), 'tesis/codigo/data/')

#-------------------
# DATA
#-------------------
datadict = {
    'rvrd':
    {'datafile':  os.path.join(datapath, target+'/HD40307_HARPS03(DRS-3-5).rdb'),
     'instrument': 'HARPS03',
     'mask': 'G2',
     'type': 'RV',
     'jittermodel': 'constant'
    },
}

#---------------------
# DRIFT + JITTER
#--------------------

driftdict = {'v0': [0.0, 1, ['Uniform', -1000., 1000.]], 
             'lin': [0.0, 0, ['Uniform', -10., 10.]],
             'qua': [0.0, 0],
             'cub': [0.0, 0],
             'tref': [278.019, 0],
             'unitconstant': [1.0, 0],
             }

rvrddict = {'offset': [0.0, 0, ['Uniform', -10, 10]],
             'jitter': [0.0, 1, ['ModJeffreys', 1.0, 99.0]],
             }

objectdict = {'rvrd': rvrddict,
              'drift1': driftdict,
          }

configdicts = [objectdict, datadict]
