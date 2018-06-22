target = 'eprv3rv01'
runid = 'k2'
n = 2
submit = True

import os
import numpy as np

datapath = os.path.join(os.getenv('HOME'), 'tesis/codigo/data')

#--------------------
# Dictionaries 
#--------------------
rundict = {
    'target': target,
    'runid': runid,
    'nplanets': n,
    'comment': 'sesinsecos',
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
fp1dict = {'k1' : [2.76, 2, ['ModJeffreys', 1., 999.], 0.5], 
           'period' : [0.0, 1, ['Jeffreys', 39.8107, 44.6684], 0.13],
           'ecos' : [0.1, 1, ['Normal', 0.0, 0.2]],
           'esin' : [0.1, 1, ['Normal', 0.0, 0.2]],
           # 'ecc': [0.1, 1, ['Uniform', 0., 1.]],
           # 'omega' : [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'ma0' : [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'epoch' : [278.019, 0]
}

fp2dict = {'k1' : [2.76, 2, ['ModJeffreys', 1., 999.], 0.5], 
           'period' : [0.0, 1, ['Jeffreys', 11.4815, 12.8825], 0.13],
           'ecos' : [0.1, 1, ['Normal', 0.0, 0.2]],
           'esin' : [0.1, 1, ['Normal', 0.0, 0.2]],
           # 'ecc': [0.1, 1, ['Uniform', 0., 1.]],
           # 'omega' : [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'ma0' : [0.1, 1, ['Uniform', 0., 2*np.pi]],
           'epoch' : [278.019, 0]
}
driftdict = {'v0': [0.0, 1, ['Uniform', -1000., 1000.]], 
             'lin': [0.0, 0, ['Uniform', -10., 10.]],
             'qua': [0.0, 0],
             'cub': [0.0, 0],
             'tref': [278.019, 0],
             'unitconstant': [1.0, 0],
             }

eprvdict = {'offset': [0.0, 0, ['Uniform', -10, 10]],
             'jitter': [0.0, 1, ['ModJeffreys', 1.0, 99.0]],
             }

objectdict = {'eprv': eprvdict,
              'drift1': driftdict,
              'planet1': fp1dict,
              'planet2': fp2dict
          }

configdicts = [rundict, objectdict, datadict]
