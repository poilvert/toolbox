#!/usr/bin/env python

__all__ = ['zmuv']

import numpy as np
from _zmuv import _zmuv_float32

def zmuv(arr_in):
    """Zero-Mean, Unit Variance along axis 1 of a 2 dimensional array

    Inputs
    ----------
    arr_in: ndarray, shape = [height, width]
        Input array

    Returns
    -------
    arr_out: ndarray, shape = [height, width]
    """
    assert arr_in.ndim == 2

    narr = np.ascontiguousarray(arr_in.copy(), dtype='f')
    arr_out = _zmuv_float32(narr)

    return arr_out
