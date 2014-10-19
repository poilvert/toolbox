#!/usr/bin/env python


import numpy as np
from math import pi
from crystal_utilities import get_elementary_cell_volume
from crystal_utilities import get_reciprocal_space_elementary_vectors
from crystal_utilities import get_Euclidean_dot_product
from crystal_utilities import get_ordered_stars


def test_cartesian_cell_volume():
    basis = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.]
        ]).astype('d')
    reference_volume = 1.
    volume = get_elementary_cell_volume(basis)
    assert abs(reference_volume - volume) < 1e-6
    return

def test_reciprocal_cartesian_cell_volume():
    basis = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.]
        ]).astype('d')
    reciprocal_basis = get_reciprocal_space_elementary_vectors(basis)
    reference_volume = (2. * pi) ** 3
    volume = get_elementary_cell_volume(reciprocal_basis)
    assert abs(reference_volume - volume) < 1e-6
    return

def test_reciprocal_ZnS_basis():
    basis = 0.5 * np.array([
        [0., 1., 1.],
        [1., 0., 1.],
        [1., 1., 0.]
        ]).astype('d')
    reference_reciprocal_basis = 2*pi * np.array([
        [-1., 1., 1.],
        [1., -1., 1.],
        [1., 1., -1.]
        ]).astype('d')
    reciprocal_basis = get_reciprocal_space_elementary_vectors(basis)
    assert np.abs(reference_reciprocal_basis - reciprocal_basis).sum() < 1e-6
    return

def test_reciprocal_ZnS_volume():
    basis = 0.5 * np.array([
        [0., 1., 1.],
        [1., 0., 1.],
        [1., 1., 0.]
        ]).astype('d')
    reciprocal_basis = get_reciprocal_space_elementary_vectors(basis)
    reference_volume = 4 * (2*pi) ** 3
    volume = get_elementary_cell_volume(reciprocal_basis)
    assert abs(reference_volume - volume) < 1e-6
    return

def test_reciprocal_Euclidean_dot_product_of_ZnS():
    basis = 0.5 * np.array([
        [0., 1., 1.],
        [1., 0., 1.],
        [1., 1., 0.]
        ]).astype('d')
    reciprocal_basis = get_reciprocal_space_elementary_vectors(basis)
    rs_dot = get_Euclidean_dot_product(reciprocal_basis)
    for i in range(100):
        m1, m2, m3 = np.random.randint(0, 100, size=3)
        expected_norm_sq = (2*pi) ** 2 * (3.*(m1**2+m2**2+m3**2) -
                2.*(m1*m2+m1*m3+m2*m3))
        norm_sq = rs_dot([m1, m2, m3], [m1, m2, m3])
        assert abs(expected_norm_sq - norm_sq) < 1e-6
    return

def test_cartesian_Euclidean_dot_product():
    basis = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.]
        ]).astype('d')
    mydot = get_Euclidean_dot_product(basis)
    for i in range(100):
        m1, m2, m3 = np.random.randint(0, 100, size=3)
        n1, n2, n3 = np.random.randint(0, 100, size=3)
        expected_dot = m1*n1+m2*n2+m3*n3
        actual_dot = mydot([m1, m2, m3], [n1, n2, n3])
        assert abs(expected_dot - actual_dot) < 1e-6
    return

def test_simple_ordered_stars():
    basis = np.array([
        [1, 0, 0],
        [2, 2, 0],
        [0, 0, 1]
        ])
    stars = get_ordered_stars(basis, 0.5)
    ref_stars = np.array([
        [0,  0,  0]
        ])
    assert np.abs(ref_stars - stars).sum() < 1e-6
    stars = get_ordered_stars(basis, 1.0)
    ref_stars = np.array([
        [ 0,  0,  0],
        [-1,  0,  0],
        [ 0,  0, -1],
        [ 0,  0,  1],
        [ 1,  0,  0]
        ])
    assert np.abs(ref_stars - stars).sum() < 1e-6
    stars = get_ordered_stars(basis, 2.0)
    ref_stars = np.array([
        [ 0,  0,  0],
        [-1,  0,  0],
        [ 0,  0, -1],
        [ 0,  0,  1],
        [ 1,  0,  0],
        [-1,  0, -1],
        [-1,  0,  1],
        [ 1,  0, -1],
        [ 1,  0,  1],
        [-2,  0,  0],
        [-2,  1,  0],
        [ 0,  0, -2],
        [ 0,  0,  2],
        [ 2, -1,  0],
        [ 2,  0,  0]
        ])
    assert np.abs(ref_stars - stars).sum() < 1e-6
    return
