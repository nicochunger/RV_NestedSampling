from model_rvrd_kepler import lnlike, lnprior, preprocess
import config
import numpy as np

nplanets = 0

modelpath = 'configfiles/hd40307_dace.py'

parnames, datadict, priordict, fixedpardict = config.read_config(
    modelpath, nplanets, False)

theta = [31334.88, 0.81, -0.16, -0.04]

print(parnames)

likelihood = lnlike(theta, parnames, fixedpardict, datadict, {})

print(f'Log(L) = {-likelihood}')
