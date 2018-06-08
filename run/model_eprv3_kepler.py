import numpy as np
from scipy.linalg import cho_factor
import warnings

from math import pi

import likelihood

def trueanomaly(M, ecc, method='Newton', niterationmax=1e4):
    if not isinstance(M, float):
        E = M
    else:
        E = n.array([M,])

    Eo = M
    ecc = n.where(ecc > 0.99, 0.99, ecc)

    niteration = 0
    while n.linalg.norm(E - Eo, ord=1) > 1e-5 or niteration==0:
        Eo = E

        ff = E - ecc*n.sin(E) - M
        dff = 1 - ecc*n.cos(E)

        if method == 'Newton':
            # Use Newton method
            E = Eo - ff / dff

        elif method == 'Halley':
            # Use Halley's parabolic method
            d2ff = ecc*n.sin(E)
        
            discr = dff **2 - 2 * ff * d2ff

            E = n.where((discr < 0), Eo - dff / d2ff,
                         Eo - 2*ff / (dff + n.sign(dff) * n.sqrt(discr))
                     )

        # Increase iteration number; if above limit, break with exception.
        niteration += 1
        if niteration >= niterationmax:
            raise RuntimeError('Eccentric anomaly comoputation not converged.')
        
    # Compute true anomaly
    nu = 2. * n.arctan2(n.sqrt(1. + ecc) * n.sin(E / 2.),
                        n.sqrt(1. - ecc) * n.cos(E / 2.)
                        )
    return nu

def preprocess(rundict, initdict, datadict, priordict, fixeddict):
    """
    Compute covariance matrix of data
    """
    covdict = {}
    for inst in datadict.keys():
        t = datadict[inst]['data']['jdb'].values
        dt = t[:, None] - t[None, :]

        # Compute red part covariance according to prescription
        alpha2 = 3.0
        le = 50.0
        lp = 0.5
        tau = 20.0
        cov = alpha2 * np.exp(-0.5 * ((np.sin(pi * dt / tau)/lp)**2 +
                                      (dt/le)**2))
        covdict[inst] = cov

    return [covdict,]

    
def modelk(pardict, time, planet='1'):
    """
    Compute Keplerian curve.

    :param array-like param: Keplerian parameters (K, P, sqrt(e)*cos(w),
    sqrt(e)*sin(w), L0, v0, epoch)
    """

    try:
        K_ms = pardict['planet{}_k1'.format(planet)]
    except KeyError:
        K_ms = np.exp(pardict['planet{}_logk1'.format(planet)])

    try:
        P_day = pardict['planet{}_period'.format(planet)]
    except KeyError:
        P_day = np.exp(pardict['planet{}_logperiod'.format(planet)])

    ###
    # SESIN SECOS ML0
    if 'planet{}_secos'.format(planet) in pardict:
        secos = pardict['planet{}_secos'.format(planet)]
        sesin = pardict['planet{}_sesin'.format(planet)]
        ecc = secos**2 + sesin**2
        omega_rad = np.arctan2(sesin, secos)

    elif 'planet{}_ecos'.format(planet) in pardict:
        ecos = pardict['planet{}_ecos'.format(planet)]
        esin = pardict['planet{}_esin'.format(planet)]
        ecc = np.sqrt(ecos**2 + esin**2)
        omega_rad = np.arctan2(esin, ecos)
        
    else:
        try:
            ecc = pardict['planet{}_ecc'.format(planet)]
            omega_rad = pardict['planet{}_omega'.format(planet)]
        except KeyError:
            raise KeyError('Something is wrong with the eccentricity parametrisation')

    # Compute mean anomaly at epoch
    # ml0 = pardict['planet{}_ml0'.format(planet)]
    # ma0 = ml0 - omega_rad

    ###
    ma0_rad = pardict['planet{}_ma0'.format(planet)]
    epoch = pardict['planet{}_epoch'.format(planet)]

    # Compute mean anomaly
    ma = 2*np.pi/P_day * (time - epoch) + ma0_rad

    # Compute true anomaly
    nu = trueanomaly(ma, ecc)

    return K_ms * ( np.cos(nu + omega_rad) + ecc * np.cos(omega_rad) )

def model(pardict, time):

    v0 = pardict['drift1_v0']
    lin = pardict['drift1_lin']
    qua = pardict['drift1_qua']
    cub = pardict['drift1_cub']

    # Find out number of planets
    planets = []
    for i in range(1, 10):
        if ('planet{}_period'.format(i) in pardict) or \
           ('planet{}_logperiod'.format(i) in pardict):
            planets.append(i)

    rv_planet = np.zeros((len(planets) + 1, len(time)))
    for i, planet in enumerate(planets):
        rv_planet[i] = modelk(pardict, time, planet=planet)
        
    # Add secular acceleration
    tt = (time - pardict['drift1_tref'])/365.25
    rv_planet[-1] = lin*tt + qua*tt**2 + cub*tt*3

    return v0 + rv_planet.sum(axis=0)
    
def lnlike(param, parnames, fixedpardict, data, covdict, **kwargs):

    # Cast param in 2D in all cases (parameters run column-wise)
    #param = np.atleast_2d(param)
    pardict = dict((par, param[parnames.index(par)]) for par in parnames)
    
    # Add fixed parameters to pardict
    pardict.update(fixedpardict)

    lnlike = np.zeros(len(data))
    
    for i, instrument in enumerate(data.keys()):
        t = data[instrument]['data']['jdb'].values
        y = data[instrument]['data']['vrad'].values
        ey = data[instrument]['data']['svrad'].values

        # Get pre-computed covariance matrix
        cov = covdict[instrument].copy()

        # Add diagonal with noise.
        jitter  = pardict['{}_jitter'.format(instrument)]
        cov += np.diag(ey**2 + jitter**2)
        
        # Construct model
        try:
            rvm = model(pardict, t)
        except RuntimeError:
            return -np.inf

        # Compute likelihood
        res = y - pardict['{}_offset'.format(instrument)] - rvm
        lnlike[i] = likelihood.lnlike_gaussian(res, cov)

    return np.sum(lnlike)


def lnprior(param, parnames, priordict, warning=False):
    param = np.atleast_1d(param)

    correction = np.zeros(9)
    
    # Check eccentricity constrain
    for pla in range(1, 10):
        try:
            secos = param[parnames.index('planet{}_secos'.format(pla))]
            sesin = param[parnames.index('planet{}_sesin'.format(pla))]

            # If using secos and sesin, we need to add a correcting factor
            # We can use uniforms in secos and sesin, and this keeps
            # uniforms in omega and ecc, but the density has to be corrected
            correction[pla] = np.log(2./pi)
            
        except ValueError:
            continue
        if secos**2 + sesin**2 >= 1:
            return -np.inf
    
    # Compute prior for each element
    lnpriorpdf = np.zeros_like(param)

    for i, par in enumerate(parnames):
        x = param[parnames.index(par)]
        lnpriorpdf[i] = np.log(priordict[par].pdf(x))

        if warning:
            warnings.warn('Warning! Null prior probability '\
                          'for {}.'.format(par))

    # Return final prior, including eccentricity
    return np.sum(lnpriorpdf, axis=0) + np.sum(correction)
