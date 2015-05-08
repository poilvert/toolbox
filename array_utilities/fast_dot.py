#!/usr/bin/env python

import numpy as np

def dot(A, B, out=None):
    """
    A drop in replacement for numpy.dot.
    Computes A.B optimized using fblas calls.
   """
    import scipy.linalg as sp
    gemm = sp.get_blas_funcs('gemm', arrays=(A,B))

    if out is None:
        lda, x, y, ldb = A.shape + B.shape
        if x != y:
            raise ValueError("matrices are not aligned")
        dtype = np.max([x.dtype for x in (A, B)])
        out = np.empty((lda, ldb), dtype, order='F')

    if A.flags.c_contiguous and B.flags.c_contiguous:
        gemm(alpha=1., a=A.T, b=B.T, trans_a=True, trans_b=True, c=out, overwrite_c=True)
    if A.flags.c_contiguous and B.flags.f_contiguous:
        gemm(alpha=1., a=A.T, b=B, trans_a=True, c=out, overwrite_c=True)
    if A.flags.f_contiguous and B.flags.c_contiguous:
        gemm(alpha=1., a=A, b=B.T, trans_b=True, c=out, overwrite_c=True)
    if A.flags.f_contiguous and B.flags.f_contiguous:
        gemm(alpha=1., a=A, b=B, c=out, overwrite_c=True)

    return out
