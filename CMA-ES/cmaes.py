#!/usr/bin/env python

# Authors: Nicolas Poilvert <nicolas.poilvert@gmail.com>
# License: BSD

"""
Covariance Matrix Adaptation - Evolution Strategy (CMA-ES)
Taken from a tutorial introduction to CMA-ES by Hansen.
"""

__all__ = ['cmaes']

import numpy as np
from numpy.linalg import norm
from numpy.linalg import eigh

# -- Pretty printing options for numpy arrays
np.set_printoptions(
        precision=4,
        threshold=5,
        linewidth=80,
        suppress=True,
        )

def set_parameters(n):

    assert n >= 1

    # -- Default selection and recombinaison parameters (equations 44 and 45 from
    # Hansen's CMA-ES tutorial)
    l = np.int(4 + np.floor(3. * np.log(n)))

    mup = l / 2.
    mu = np.int(np.floor(mup))
    assert mu <= mup <= l

    w = np.log(mup + 0.5) - np.log(np.arange(1, mu + 1))
    w /= w.sum()
    assert w.size == mu
    assert np.abs(1. - np.sum(w)) < 1e-5

    # -- Default step-size control parameters. Eq. 46 from Hansen's
    mueff = (w.sum()) ** 2 / (w ** 2).sum()
    assert mueff >= 1.

    cs = (mueff + 2.) / (n + mueff + 5.)
    ds = 1. + 2. * np.max(0., np.sqrt((mueff - 1.)/(n + 1.)) - 1.) + cs

    # -- Default parameters for Covariance matrix adaptation. Eqs. 47, 48 and 49
    cc = (4. + mueff / n) / (n + 4. + 2. * mueff / n)
    c1 = 2. / ((n + 1.3) ** 2 + mueff)
    cm = 2. * ((mueff - 2. + 1. / mueff)/((n + 2.) ** 2 + mueff))

    return l, mu, w, cs, ds, cc, c1, cm, mueff

def initialization(
        n,
        range_min=-1.,
        range_max=1.
        ):

    # -- Natural input vector coordinate range (this will be problem dependent)
    assert range_min <= range_max

    # -- Core parameters of CMA-ES algorithm
    sigma = 0.3 * (range_max - range_min)
    m = np.random.uniform(low=range_min, high=range_max, size=n)
    ps = np.zeros(n)
    pc = np.zeros(n)
    B = np.eye(n)
    D = np.eye(n)
    BD = np.dot(B, D)
    C = np.dot(BD, BD.T)

    return m, sigma, ps, pc, B, D, C

def new_sample_population(m, B, D, sigma, l):
    """
    returns the new sample of ``l`` n-dimensional sample vectors
    """

    assert m.ndim == 1
    assert B.ndim == 2
    assert B.shape == D.shape
    assert B.shape[0] == B.shape[1]
    assert m.size == B.shape[0]

    n = m.size

    # -- generating new population search points
    z = np.random.normal(size=(n, l))
    y = np.dot(B, np.dot(D, z))
    x = m[:, None] + sigma * y

    return z, x

def update_mean(x, z, f_x, w):

    assert x.ndim == 2
    assert x.shape == z.shape
    assert len(f_x) == x.shape[1]
    assert w.ndim == 1
    assert w.size <= x.shape[1]

    argsorted_idx = np.argsort(f_x)
    mu = w.size

    xs = x[:, argsorted_idx[:mu]]
    zs = z[:, argsorted_idx[:mu]]

    xm = np.dot(xs, w)
    zm = np.dot(zs, w)

    assert xm.shape == (x.shape[0],)
    assert zm.shape == (x.shape[0],)

    return xm, zm, zs

def chin(n):
    """
    Well-behaved approximation to : sqrt(2)*Gamma((n+1)/2)/Gamma(n/2)
    """
    return np.sqrt(n)*(1. - 1./(4.*n) + 1./(21.*n**2))

def cmaes(
        n, func, stopping_func,
        range_min=-1.,
        range_max=1.
        ):
    """
    ``func`` is the function to optimize (here minimize). ``n`` is the dimension
    of the input vectors fed to ``func``. ``stopping_func`` is the function that
    determines if the main optimization loop should terminate. That latter
    function should evaluate to True or False.
    """

    # -- Set parameters to default values given the dimensionality of input
    # vectors
    l, mu, w, cs, ds, cc, c1, cm, mueff = set_parameters(n)

    # -- Initialization
    m, sigma, ps, pc, B, D, C = initialization(n,
            range_min=range_min, range_max=range_max)

    # -- Main loop
    n_func_eval = 0
    while stopping_func(n_func_eval):

        print 'search vector:', m

        # -- Sampling a new population of search points
        z, x = new_sample_population(m, B, D, sigma, l)

        # -- Computing the "fitness" or "cost" for each sample
        f_x = [func(vec) for vec in x.T]
        n_func_eval += len(f_x)

        # -- Update optimal mean search point
        m, zm, zs = update_mean(x, z, f_x, w)

        # -- Cumulation: update evolution path
        ps = (1. - cs)*ps + np.sqrt(cs*(2. - cs)*mueff) * np.dot(B, zm)
        hs = int(
            (norm(ps) / (chin(n) * np.sqrt(1.-(1. - cs)**(2*n_func_eval/l)))) < \
            (1.4 + 2./(n + 1.))
        )
        pc = (1. - cc)*pc + \
             hs*np.sqrt(cc*(2. - cc)*mueff) * np.dot(B, np.dot(D, zm))

        # -- Covariance Matrix Adaptation
        y = np.dot(B, np.dot(D, zs))
        C = (1. - c1 - cm) * C + \
            c1 * (np.outer(pc, pc) + (1. - hs)*cc*(2. - cc)*C) + \
            cm * np.dot(y, np.dot(np.diag(w), y.T))

        # -- Update step size
        sigma *= np.exp(cs/ds * (norm(ps)/chin(n) - 1.))

        # -- Update intermediate matrices B and D from C
        D, B = eigh(C)

        assert D.ndim == 1

        np.putmask(D, D <= 0., 0.)
        D = np.diag(np.sqrt(D))

    return m
