"""
A module implementing different likelihood functions.
"""

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from math import pi, log

def lnlike_gaussian(residuals, K):
    """
    Compute likelihood for normally distributed errors.

    :param np.array residuals: residuals vector of dimension (n,).
    :param np.array K: (n x n) covariance matrix.
    :return float: likelihood value.
    """

    # Cholesky factor of the covariance
    factor, flag = cho_factor(K, overwrite_a=False)

    # Determinant
    lndet = 2 * np.sum(np.log(np.diag(factor)))

    alpha = cho_solve((factor, flag), residuals)

    # Constant contribution
    nobs = np.shape(K)[0]
    cte = nobs * np.log(2*np.pi)

    # TODO check implementation for multi-residuals using np.diag.
    # return -0.5 * (lndet + np.diag(np.dot(residuals.T, alpha)))
    return -0.5 * (lndet + np.dot(residuals.T, alpha) + cte)


def lnlike_gaussian_invk(residuals, invk, lndet):
    """
    Compute likelihood for normally distributed errors.

    :param np.array residuals: residuals vector of dimension (n,).
    :param np.array invk: (n x n) inverse covariance matrix.
    :param float lndet: log determinant of the covariance matrix.

    :return float: likelihood value.
    """
    return -0.5 * (lndet + np.dot(residuals.T, np.dot(invk, residuals)))


def lnlike_gaussian_marginalisedjitter(residuals, xnodes, wnodes, chofactors,
                                       lndeterminants, priorpdf, a=0, b=99.0):
    """
    Compute likelihood for normally distributed errors, with additional
    white noise term marginalised out.

    The marginalisation is doing through Gauss-Legendre quadrature using the
    nodes, weights and already factorised covariance functions provided.
    """
    lnlike = np.zeros(len(xnodes))

    for i, _ in enumerate(xnodes):
        # Compute value of likelihood and prior in nodes.

        # Solve equation to get K^{-1} r
        alpha = cho_solve(chofactors[i], residuals)

        # Compute likelihood
        lnlike[i] = -0.5 * np.dot(residuals.T, alpha)

    lnmarglike = log_sum(np.log(wnodes) + lndeterminants +
                             lnlike + np.log(priorpdf))
    return log(0.5*(b-a)) + lnmarglike
        
        
def log_sum(log_summands):
    a = np.inf
    x = log_summands.copy()
    while a == np.inf or a == -np.inf or np.isnan(a):
        a = x[0] + np.log(1 + np.sum(np.exp(x[1:] - x[0])))
        np.random.shuffle(x)
    return a

