target = 'eprv3rv01'
t2 = target.replace('dot', '.')
runid = 'k0'
n = 3
commentstr = ''
submit = True

import os
import pickle
import numpy as np

datapath = os.path.join(os.getenv('HOME'), 'tesis/codigo/data')

#--------------------
# Dictionaries 
#--------------------
rundict = {
    'target': target,
    'runid': runid,
    'nplanets': n,
    'comment': '',
    'sampler': 'emcee',
    'nwalkers': 300,
    'a': 2.0
            }

#-------------------
# DATA
#-------------------
datadict = {
    'eprv':
    {'datafile':  os.path.join(datapath, 'rvs_0001.txt'),
     'instrument': 'HARPS',
     'mask': 'G2',
     'type': 'RV',
     'jittermodel': 'constant'
 },
    }

#---------------------
#        PLANET
#--------------------
fp1dict = {'k1' : [3.5, 1, ['Uniform', 0., 10.]], 
           'period' : [4.6514, 2, ['Uniform', 4.0, 5.2], 3e-3],
           'sesin' : [0.1, 1, ['Uniform', -1., 1.]],
           'secos' : [0.1, 1, ['Uniform', -1., 1.]],
           'ml0' : [0.0, 1, ['Uniform', -6*np.pi, 6*np.pi]],
           'epoch' : [56845.73, 0]
}

fp2dict = {'k1' : [3.5, 1, ['Uniform', 0., 10.]], 
           'period' : [3.06, 2, ['Uniform', 2.5, 3.5], 3e-3],
           'sesin' : [0.1, 1, ['Uniform', -1., 1.]],
           'secos' : [0.1, 1, ['Uniform', -1., 1.]],
           'ml0' : [0.0, 1, ['Uniform', -6*np.pi, 6*np.pi]],
           'epoch' : [56845.73, 0]
}

driftdict = {'v0': [0.0, 1, ['Uniform', -1000., 1000.]], 
             'lin': [0.0, 0, ['Uniform', -10., 10.]],
             'qua': [0.0, 0],
             'cub': [0.0, 0],
             'tref': [56845.73, 0],
             'unitconstant': [1.0, 0],
             }

eprvdict = {'offset': [0.0, 0, ['Uniform', -10, 10]],
             'jitter': [0.0, 1, ['ModJeffreys', 1.0, 99.0]],
             }

objectdict = {'eprv': eprvdict,
              'drift1': driftdict,
              'planet1': fp1dict,
              'planet2': fp2dict,
          }

configdicts = [rundict, objectdict, datadict]
