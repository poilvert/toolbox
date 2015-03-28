#!/usr/bin/env python

# Authors: Nicolas Poilvert <nicolas.poilvert@gmail.com>
#          Nicolas Pinto <nicolas.pinto@gmail.com>
# license: BSD

"""
Some utility functions for memory management
"""

__all__ = ['get_homogeneous_block_shape', 'block_coordinates']

import numpy as np

def get_homogeneous_block_shape(arr_shape, filter_shape,
                                mem_limit=1024**3,
                                itemsize=4):
    """
    This routine will compute the shape of a canonical
    block (hyperrectangle of the same dimensionality than
    the input array) in such a way that the memory footprint
    does not go beyond a fixed limit.

    This canonical block is a homogeneous scaled down version
    of the array.

    Parameters
    ----------

    ``arr_shape``: iterable
        shape of the array from which one wants to extract
        canonical blocks

    ``filter_shape``: iterable
        shape of the filter used in the correlation operation

    ``mem_limit``: int
        maximum memory footprint (in Bytes) for a canonical
        block. Default is 1GB

    ``itemsize``: int
        number of Bytes used to store one element of the array.
        Default is 4, corresponding to float32.

    Returns
    -------

    ``blk_shape``: tuple
        shape of the canonical block
    """

    from math import floor

    # -- checks and transforms on inputs
    arr_shape = np.array(arr_shape).astype('f')
    filter_shape = np.array(filter_shape).astype('f')
    mem_limit = float(mem_limit)
    itemsize = float(itemsize)

    assert (arr_shape > 0).all()
    assert (filter_shape > 0).all()

    # -- dimensionality
    N = float(arr_shape.size)
    assert N >= 1.

    # -- homogeneous scaling factor
    alpha = (mem_limit /
            (itemsize * np.prod(filter_shape) *
                        np.prod(arr_shape))) ** (1. / N)
    assert alpha >= 0.

    # -- if alpha is greater than one, then the block can
    #    be chosen to be the entire array
    if alpha >= 1.:
        blk_shape = tuple(arr_shape.astype(int))
    else:
        blk_shape = ()
        for i in xrange(int(N)):
            blk_shape += (max(1, int(floor(alpha * arr_shape[i]))),)

    return blk_shape

def block_coordinates(arr_shape, blk_shape):
    """
    This algorithm will give the coordinates (so the starting
    and stopping indices) of all the blocks of shape ``blk_shape``
    contained in an array of shape ``arr_shape``.

    The blocks are not necessarily of the same shape in the end
    because the block shape is not necessarily commensurate with
    the array shape.

    Parameters
    ----------

    ``arr_shape``: iterable
        shape of the array (which is the result of the
        correlation operation)

    ``blk_shape``: iterable
        shape of the canonical block with which one wants to
        pave the array

    Returns
    -------

    ``start``, ``stop``: 2D arrays
        arrays of size (``len(arr_shape)``, ``n_blks``)
        where ``n_blks`` is the total number of canonical blocks
        that fit into the array
    """

    # -- checks and transforms on inputs
    arr_shape, blk_shape = np.array(arr_shape).astype(int), \
                           np.array(blk_shape).astype(int)

    assert arr_shape.ndim == 1
    assert blk_shape.ndim == 1
    assert arr_shape.size == blk_shape.size

    assert ((arr_shape - blk_shape) >= 0).all()
    assert (arr_shape > 0).all()
    assert (blk_shape > 0).all()

    # -- dimensionality of input array
    N = arr_shape.size

    # -- total number of expected blocks in input array
    n_blks = 1
    for oi, hi in zip(arr_shape, blk_shape):
        n_blks *= (1 + (oi - 1) / hi)

    # -- 1D arrays of starting and stopping block coordinates
    #    for every dimension in the array
    start_list, stop_list = [], []
    for oi, hi in zip(arr_shape, blk_shape):
        start_idx = np.arange(0, oi, hi)
        stop_idx = np.concatenate([np.arange(hi, oi, hi), np.array([oi])])
        start_list += [start_idx]
        stop_list += [stop_idx]

    # -- augmenting the dimensionality of each array to
    #    dimension N
    for idx, (start_arr, stop_arr) in enumerate(zip(start_list, stop_list)):
        start_arr.shape = max(0, idx - 1) * (1,) + (start_arr.size,) \
                + max(0, N - idx) * (1,)
        stop_arr.shape = max(0, idx - 1) * (1,) + (stop_arr.size,) \
                + max(0, N - idx) * (1,)

    # -- broadcasting the array to full N dimensional arrays
    start_list = np.broadcast_arrays(*start_list)
    stop_list = np.broadcast_arrays(*stop_list)

    # -- adding a first extra dimension (for later concatenation along
    #    that dimension)
    for idx, (start_arr, stop_arr) in enumerate(zip(start_list, stop_list)):
        start_list[idx] = start_arr[np.newaxis, ...]
        stop_list[idx] = stop_arr[np.newaxis, ...]

    # -- concatenation along the first (extra) dimension and
    #    reshaping to get a 2D array
    start = np.concatenate(start_list, axis=0).reshape(N, -1)
    stop = np.concatenate(stop_list, axis=0).reshape(N, -1)

    start = np.atleast_2d(start)
    stop = np.atleast_2d(stop)

    assert start.shape == stop.shape
    assert start.shape == (N, n_blks)

    return (start, stop)
