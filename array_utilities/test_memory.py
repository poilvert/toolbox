#!/usr/bin/env python

"""
Test suite for the ``memory`` module
"""

import numpy as np

from memory import get_homogeneous_block_shape
from memory import block_coordinates


def test_1D_arr_small_enough():
    N = 6
    mem_limit = 100
    d = 3
    out = get_homogeneous_block_shape((N,), (d,),
            mem_limit=mem_limit)
    assert out == (N,)


def test_1D_arr_just_right():
    N = 8
    mem_limit = 100
    d = 3
    out = get_homogeneous_block_shape((N,), (d,),
            mem_limit=mem_limit)
    assert out == (N,)


def test_1D_arr_too_large():
    N = 9
    mem_limit = 100
    d = 3
    out = get_homogeneous_block_shape((N,), (d,),
            mem_limit=mem_limit)
    assert out == (8,)


def test_1D_arr_unit_filter():
    N = 15
    mem_limit = 15
    itemsize = 1
    d = 1
    out = get_homogeneous_block_shape((N,), (d,),
            mem_limit=mem_limit, itemsize=itemsize)
    assert out == (N,)


def test_2D_arr():
    N, M = 8, 8
    f = 13
    arr_shape = (N, M)
    filter_shape = (f,)
    d1, d2 = 2, 2
    mem_limit = 4 * d1 * d2 * f
    out = get_homogeneous_block_shape(arr_shape,
                                      filter_shape,
            mem_limit=mem_limit)
    assert out == (d1, d2)


def test_blk_coord_1D():
    arr_shape = (11,)
    blk_shape = (3,)
    ref_start = np.array([[0, 3, 6, 9]])
    ref_stop = np.array([[3, 6, 9, 11]])
    start, stop = block_coordinates(arr_shape, blk_shape)
    assert (start == ref_start).all()
    assert (stop == ref_stop).all()


def test_blk_coord_2D():
    arr_shape = (5,4)
    blk_shape = (3,2)
    ref_start = np.array([[0, 0, 3, 3],
                          [0, 2, 0, 2]])
    ref_stop = np.array([[3, 3, 5, 5],
                         [2, 4, 2, 4]])
    start, stop = block_coordinates(arr_shape, blk_shape)
    assert (start == ref_start).all()
    assert (stop == ref_stop).all()


def test_blk_coord_unit_blk_shape():
    arr_shape = (1,2,3)
    blk_shape = (1,1,1)
    ref_start = np.array([[0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 1, 1],
                          [0, 1, 2, 0, 1, 2]])
    ref_stop = np.array([[1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 2, 2, 2],
                         [1, 2, 3, 1, 2, 3]])
    start, stop = block_coordinates(arr_shape, blk_shape)
    assert (start == ref_start).all()
    assert (stop == ref_stop).all()
