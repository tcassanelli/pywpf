# Implementation of the methods described in
#
# Hongwei Guo, "A Simple Algorithm for Fitting a Gaussian Function"
# IEEE Signal Processing Magazine, September 2011, pp. 134--137
#
# Author: Stefan van der Walt, 2011

from __future__ import division

import numpy as np
import scipy as sp
import scipy.linalg as sl
import warnings

warnings.simplefilter("error")

def gauss(x, A=1, mu=1, sigma=1):
    """
    Evaluate Gaussian.
    
    Parameters
    ----------
    A : float
        Amplitude.
    mu : float
        Mean.
    std : float
        Standard deviation.

    """
    return np.real(A * np.exp(-(x - mu)**2 / (2 * sigma**2)))

def fit_direct(x, y, F=0, weighted=True, _weights=None):
    """Fit a Gaussian to the given data.

    Returns a fit so that y ~ gauss(x, A, mu, sigma)

    Parameters
    ----------
    x : ndarray
        Sampling positions.
    y : ndarray
        Sampled values.
    F : float
        Ignore values of y &lt;= F.
    weighted : bool
        Whether to use weighted least squares.  If True, weigh
        the error function by y, ensuring that small values
        has less influence on the outcome.

    Additional Parameters
    ---------------------
    _weights : ndarray
        Weights used in weighted least squares.  For internal use
        by fit_iterative.

    Returns
    -------
    A : float
        Amplitude.
    mu : float
        Mean.
    std : float
        Standard deviation.

    """
    mask = (y >= F)
    x = x[mask]
    y = y[mask]

    if _weights is None:
        _weights = y
    else:
        _weights = _weights[mask]

    # We do not want to risk working with negative values
    np.clip(y, 1e-10, np.inf, y)

    e = np.ones(len(x))
    if weighted:
        e = e * (_weights**2)
    
    v = (np.sum(np.vander(x, 5) * e[:, None], axis=0))[::-1]
    A = v[sl.hankel([0, 1, 2], [2, 3, 4])]

    ly = e * np.log(y)
    ls = np.sum(ly)
    x_ls = np.sum(ly * x)
    xx_ls = np.sum(ly * x**2)
    B = np.array([ls, x_ls, xx_ls])

    (a, b, c), res, rank, s = np.linalg.lstsq(A, B)

    A = np.exp(a - (b**2 / (4 * c)))
    mu = -b / (2 * c)
    sigma = sp.sqrt(-1 / (2 * c))

    return A, mu, sigma

def fit_iterative(x, y, F=0, weighted=True, N=10):
    """Fit a Gaussian to the given data.

    Returns a fit so that y ~ gauss(x, A, mu, sigma)

    This function iteratively fits using fit_direct.
    
    Parameters
    ----------
    x : ndarray
        Sampling positions.
    y : ndarray
        Sampled values.
    F : float
        Ignore values of y &lt;= F.
    weighted : bool
        Whether to use weighted least squares.  If True, weigh
        the error function by y, ensuring that small values
        has less influence on the outcome.
    N : int
        Number of iterations.

    Returns
    -------
    A : float
        Amplitude.
    mu : float
        Mean.
    std : float
        Standard deviation.

    """
    y_ = y
    for i in range(N):
        p = fit_direct(x, y, weighted=True, _weights=y_)
        A, mu, sigma = p
        y_ = gauss(x, A, mu, sigma)

    return np.real(A), np.real(mu), np.real(sigma)

def log_gauss_param(A, mu, sigma):
    """Give A, mu, sigma, that represent the Gaussian

    gauss(x, A, mu, sigma),

    compute the a, b, c that parameterise the parabola

    log(gauss(x, A, mu, sigma)) ~ a + bx**2 + cx.

    """
    ss = 2 * sigma**2
    return np.log(A) - mu**2 / ss, \
           2 * mu / ss, \
           -1 / ss