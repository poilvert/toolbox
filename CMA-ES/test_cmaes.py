#!/usr/bin/env python

from cmaes import cmaes
import numpy as np

# -- Test parameters
n = 5
max_iterations = 2000

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
def stop(g, the_max=max_iterations):
    if g <= the_max:
        return True
    else:
        return False

def test_cmaes_random_quadratic():

    # -- instantiating the function
    myfct, opt_vec = get_func(n)

    # -- Testing!
    res = cmaes(n, myfct, stop)

    assert np.abs(res - opt_vec).sum() < 1e-3

if __name__ == "__main__":

    def get_stop(max_iterations=1000):
        def mystop(g):
            if g <= max_iterations:
                return True
            else:
                return False
        return mystop

    n = 50
    max_iterations = 2e1*n**2

    myfct, opt_vec = get_func(n)
    stop = get_stop(max_iterations=max_iterations)

    res = cmaes(n, myfct, stop, range_min=0,
            range_max=0.3*n)
