import imp
import numpy as np
import pandas as pd
import copy
import os
from pprint import pprint
import priors


def read_config(configfile, nplanets, dfile, narrow):
    """
    Initialises a sampler object using parameters in config.

    :param string config: string indicating path to configuration file.
    :param int nplanets: integer with the number of planets to be used in the model
    :param int dfile: number data set to use
    :param float narrow: Float indicating how wide to take the narrow priors.
                         If it is 0, broad priors will be used instead.
    """

    # Import configuration file as module
    c = imp.load_source('c', configfile)

    # Make copy of all relavant dictionaries and lists
    datadict, fpdict, driftdict, eprvdict = map(dict.copy, c.configdicts)
    narrow_priors = dict.copy(c.narrow_priors)

    # Write correct data file path
    datapath = os.path.join(os.getenv('HOME'), 'tesis/codigo/data')
    datadict['eprv']['datafile'] = datapath + '/rvs_000{}.txt'.format(dfile)

    # Create input_dict in acordance to number of planets in the model
    input_dict = {'eprv': eprvdict, 'drift1': driftdict}
    for i in range(1, nplanets+1):
        input_dict.update({'planet{}'.format(i): copy.deepcopy(fpdict)})
        # If narrow is True change period priors to the narrow values
        if narrow:
            input_dict['planet{}'.format(
                i)]['period'][2][1] = narrow_priors[dfile][i-1][0]
            input_dict['planet{}'.format(
                i)]['period'][2][2] = narrow_priors[dfile][i-1][1]

    pprint(input_dict)
    # Create prior instances
    priordict = priors.prior_constructor(input_dict, {})

    # Build list of parameter names
    parnames, _ = get_parnames(input_dict)

    # Read data from file(s)
    read_data(c.datadict)

    # Fixed parameters
    fixedpardict = get_fixedparvalues(input_dict)

    return parnames, datadict, priordict, fixedpardict


def get_parnames(input_dict):
    parnames = []
    fixparnames = []
    for obj in input_dict:
        for par in input_dict[obj]:
            if input_dict[obj][par][1] > 0:
                parnames.append(obj+'_'+par)
            elif input_dict[obj][par][1] == 0:
                fixparnames.append(obj+'_'+par)
    return parnames, fixparnames


def get_fixedparvalues(input_dict):
    fpdict = {}
    for obj in input_dict:
        for par in input_dict[obj]:
            if input_dict[obj][par][1] == 0:
                fpdict[obj+'_'+par] = input_dict[obj][par][0]
    return fpdict


def read_data(datadict):
    for inst in datadict:
        # Try to get custom separator
        try:
            sep = datadict[inst]['sep']
        except KeyError:
            # sep = '\t'
            sep = ' '

        # Read rdb file
        data = pd.read_csv(datadict[inst]['datafile'], sep=sep,
                           comment='#', header=None, names=['jbl', 'vrad', 'svrad'], usecols=[0, 1, 2])
        datadict[inst]['data'] = data
