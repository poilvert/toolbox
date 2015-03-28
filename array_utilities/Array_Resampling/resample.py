#!/usr/bin/env python

# Authors: Nicolas Poilvert <nicolas.poilvert@gmail.com>
#          Nicolas Pinto <nicolas.pinto@gmail.com>
# License: BSD

__all__ = ['resample']

import numpy as np
from _resample import _resample_float32

def resample(arr_in, out_shape):
    """Resample a 3-dimensional array to the desired shape.

    Inputs
    ----------
    arr_in: ndarray, shape = [height, width, depth]
        Input array

    out_shape: tuple of int
        Desired output shape after resampling.
        Format = [new_height, new_width, new_depth]

    Returns
    -------
    arr_out: ndarray, shape = `out_shape`
        Resampled input array of shape `out_shape`.
    """
    assert arr_in.ndim == 3
    assert len(out_shape) == arr_in.ndim

    h_in, w_in, d_in = arr_in.shape
    h_out, w_out, d_out = out_shape

    narr = np.ascontiguousarray(arr_in.copy(), dtype='f')

    arr_out = np.empty(out_shape, dtype=narr.dtype)
    _resample_float32(narr, arr_out)

    return arr_out
