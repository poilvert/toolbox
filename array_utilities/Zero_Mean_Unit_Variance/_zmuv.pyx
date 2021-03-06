# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: cdivision=True

import numpy as np
cimport numpy as np

from libc.math cimport sqrt, pow
from cython.parallel import prange

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t


# -- square root C function
cdef float _sqrt(float x):
        return sqrt(x)


def _zmuv_float32(
    np.ndarray[DTYPE_t, ndim=2] arr
    ):

    cdef Py_ssize_t rows = arr.shape[0]
    cdef Py_ssize_t cols = arr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] mean = np.zeros(rows, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] std = np.zeros(rows, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] std_inv = np.zeros(rows, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] out = np.zeros((rows, cols), dtype=DTYPE)

    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef DTYPE_t scaler = 1. / <DTYPE_t>cols

    # -- computing the vector of mean values
    for i in prange(rows, nogil=True):

        for j in range(cols):
            mean[i] += arr[i, j]
        mean[i] *= scaler

    # -- now we compute the vector of standard deviations
    for i in prange(rows, nogil=True):

        for j in range(cols):
            std[i] += (arr[i, j] - mean[i]) ** 2
        std[i] *= scaler

    for i in range(rows):

        std_inv[i] = 1. / (<DTYPE_t>_sqrt(<float>std[i]))

    # -- computing the final output array (zero mean, unit variance of the input)
    for i in prange(rows, nogil=True):

        for j in range(cols):
            out[i, j] = std_inv[i] * (arr[i, j] - mean[i])

    return out
