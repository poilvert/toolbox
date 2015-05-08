#!/usr/bin/env python

from cmaes import cmaes
import numpy as np
import cPickle

# -- Test parameters
N_L = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200]
EPS = 1e-3

# -- Generate a random quadratic function taking n-dimensional vectors as input
def get_func(n):
    H = np.random.randn(n, n).astype('f')
    H = np.dot(H, H.T)
    opt_vec = np.arange(n)
    def my_func(x):
        y = x - opt_vec
        return np.dot(y, np.dot(H, y))
    return my_func, opt_vec

# -- A simple stopping criterion function
def stop_factory(opt_vec=None, eps=EPS):

    def stop(g, m, opt_vec=opt_vec, eps=eps):
        assert m.ndim == 1
        assert opt_vec is not None

        max_g = 2e2 * m.size ** 2

        if g > max_g:
            return False

        error = np.abs(m - opt_vec).mean()
        if error > eps:
            return True
        else:
            return False

    return stop

# -- Running convergence tests with dimensionality
convergence_histories = {}
for n in N_L:
    key = '%i' % n
    myfct, opt_vec = get_func(n)
    mystop = stop_factory(opt_vec=opt_vec)
    m_opt, m_history = cmaes(n, myfct, mystop, range_min=0,
                      range_max=0.3*n)
    convergence_histories[key] = m_history

# -- dump data
cPickle.dump(convergence_histories, open('convergence_history.pkl', 'w'))
