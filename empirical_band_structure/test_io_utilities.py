#!/usr/bin/env python

from io_utilities import read_pseudo
import numpy as np
import os


def test_read_pseudo():

    filename = 'test_pseudo.dat'

    with open('test_pseudo.dat', 'w') as fin:
        fin.write('#\n')
        fin.write('3  4  5  0.2  0.1\n')
        fin.write('  4  1  87  0.99  -0.27\n')

    reference_G_vectors = np.array([
        [3, 4,  5],
        [4, 1, 87]
        ]).astype('int')

    reference_coefficients = np.array([
        0.2+0.1j, 0.99-0.27j]).astype('complex')

    G_vectors, coefficients = read_pseudo(filename)

    assert np.abs(reference_G_vectors - G_vectors).sum() < 1e-6
    assert np.abs(reference_coefficients - coefficients).sum() < 1e-6

    os.remove(filename)
    return
