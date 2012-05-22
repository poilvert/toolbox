#!/usr/bin/env python

# Authors: Nicolas Pinto <nicolas.pinto@gmail.com>
#          Nicolas Poilvert <nicolas.poilvert@gmail.com>

# Licence: BSD

"""
A simple program to generate a filter made of a Laplacian of Gaussian
"""

import numpy as np

def generate_log(filter_shape, sigma_x=1.4, sigma_y=1.4):

    filter_shape = tuple(filter_shape)
    sigma_x = np.float32(sigma_x)
    sigma_y = np.float32(sigma_y)

    assert sigma_x > 0.
    assert sigma_y > 0.

    assert len(filter_shape) == 2
    assert filter_shape[0] > 0
    assert filter_shape[1] > 0

    height = float(filter_shape[0])
    width = float(filter_shape[1])

    X, Y = np.mgrid[-height/2.:height/2.:complex(height),
                   -width/2.:width/2.:complex(width)].astype(np.float32)

    sx2 = np.float32(1. / sigma_x ** 2)
    sy2 = np.float32(1. / sigma_y ** 2)

    XX = sx2 * X ** 2
    YY = sy2 * Y ** 2

    out = (sx2 * XX + sy2 * YY - sx2 - sy2) * np.exp(- 0.5 * (XX + YY))

    out -= out.mean()

    return out
