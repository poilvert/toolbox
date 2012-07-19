#!/usr/bin/env python

# Authors: Nicolas Pinto <nicolas.pinto@gmail.com>
#          Nicolas Poilvert <nicolas.poilvert@gmail.com>

# Licence: BSD

"""
A library that contains a bunch of transformation functions directly
usable by ``transform_image``.

The ``transform`` functions should have the following pre and post
conditions :
    - these functions should take 1D arrays of x and y coordinates
    as inputs.
    - these functions should return 1D arrays of x and y output coordinates.
"""


import numpy as np

EPSILON = 1e-5


def translate(x, y, cx=20., cy=20.):

    assert x.ndim == 1
    assert y.ndim == 1
    assert x.dtype == y.dtype

    input_dtype = x.dtype

    x_out = x + cx
    y_out = y + cy

    return x_out.astype(input_dtype), y_out.astype(input_dtype)


def inverse_translate(x, y, cx=20., cy=20.):

    assert x.ndim == 1
    assert y.ndim == 1
    assert x.dtype == y.dtype

    input_dtype = x.dtype

    x_out = x - cx
    y_out = y - cy

    return x_out.astype(input_dtype), y_out.astype(input_dtype)


def inplane_rotate(x, y, xc=0., yc=0., theta=np.pi/4.):

    assert x.ndim == 1
    assert y.ndim == 1
    assert x.dtype == y.dtype

    input_dtype = x.dtype

    xp = x - xc
    yp = y - yc

    xp_out = np.cos(theta) * xp - np.sin(theta) * yp
    yp_out = np.sin(theta) * xp + np.cos(theta) * yp

    x_out = xp_out + xc
    y_out = yp_out + yc

    return x_out.astype(input_dtype), y_out.astype(input_dtype)


def inverse_inplane_rotate(x, y, xc=0., yc=0., theta=np.pi/4.):

    assert x.ndim == 1
    assert y.ndim == 1
    assert x.dtype == y.dtype

    input_dtype = x.dtype

    xp = x - xc
    yp = y - yc

    xp_out = np.cos(-theta) * xp - np.sin(-theta) * yp
    yp_out = np.sin(-theta) * xp + np.cos(-theta) * yp

    x_out = xp_out + xc
    y_out = yp_out + yc

    return x_out.astype(input_dtype), y_out.astype(input_dtype)


def zoom(x, y, xc=0., yc=0., s=2.):

    assert x.ndim == 1
    assert y.ndim == 1
    assert x.dtype == y.dtype

    input_dtype = x.dtype

    xp = x - xc
    yp = y - yc

    xp_out = s * xp
    yp_out = s * yp

    x_out = xp_out + xc
    y_out = yp_out + yc

    return x_out.astype(input_dtype), y_out.astype(input_dtype)


def inverse_zoom(x, y, xc=0., yc=0., s=2.):

    assert x.ndim == 1
    assert y.ndim == 1
    assert x.dtype == y.dtype

    assert np.abs(s) > EPSILON

    input_dtype = x.dtype

    xp = x - xc
    yp = y - yc

    xp_out = (1. / s) * xp
    yp_out = (1. / s) * yp

    x_out = xp_out + xc
    y_out = yp_out + yc

    return x_out.astype(input_dtype), y_out.astype(input_dtype)
