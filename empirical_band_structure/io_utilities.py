#!/usr/bin/env python

from os import path
import numpy as np

# -- pseudo-potential file reading utility
def read_pseudo(filename):
    """
    Input
    =====
    filename
        string giving the path to the pseudo-potential file

    Returns
    =======
    G_vectors
        2D array of shape (N, 3) where N is the number of G vectors for which
        the pseudo-potential has a non-zero Fourier coefficient

    G_coefficients
        1D array of shape (N,)

    Conventions
    ===========
    The expected format for the pseudo-potential file is the following::

        G1  G2  G3  V_G.real  V_G.imag
        ...

    where the first three objects are integer numbers corresponding to the
    integer coordinates of a G vector, ``V_G.real`` corresponds to the real part
    of the pseudo-potential coefficient for that G vector and ``V_G.imag``
    corresponds to its imaginary part.

    Any line that starts with a ``#`` will be considered as a comment.
    """
    assert path.exists(filename)
    assert path.isfile(filename)

    raw_data = np.loadtxt(filename, dtype='d', comments='#')

    assert raw_data.ndim == 2
    assert raw_data.shape[1] == 5

    n_Gs = raw_data.shape[0]

    G_vectors = raw_data[:, :3].astype('int')

    G_coefficients = np.zeros(n_Gs, dtype='complex')
    G_coefficients.real = raw_data[:, 3]
    G_coefficients.imag = raw_data[:, 4]

    return G_vectors, G_coefficients


def write_pseudo(G_vectors, G_coefficients, filename):
    """
    Given an array of G vector coordinates and the corresponding Fourier
    coefficients, this routine simply writes down the pseudo-potential to file
    """
    assert G_vectors.ndim == 2
    assert G_coefficients.ndim == 1
    assert G_vectors.shape[1] == 3
    assert G_vectors.shape[0] == G_coefficients.size
    assert G_coefficients.dtype == 'complex'

    with open(filename, 'w') as fin:
        for G, G_coeff in zip(G_vectors, G_coefficients):
            G0, G1, G2 = G
            coef_real, coef_imag = G_coeff.real, G_coeff.imag
            fin.write("%6i  %6i  %6i  %10.6f  %10.6f\n" % (G0, G1, G2,
                coef_real, coef_imag))
    return


# -- kpoint file reading utility
def read_kpoints(filename):
    """
    Input
    =====
    filename
        string giving the path to the file containing a list of k points

    Returns
    =======
    k_vectors
        2D array of shape (N, 3) where N is the number of k points in the list

    Conventions
    ===========
    The expected format for the file containing the k points is as follows::

        k1  k2  k3
        ...

    where k1, k2 and k3 are the coordinates of a k point in the basis of
    reciprocal-space elementary vectors.

    Any line that starts with a ``#`` will be considered as a comment.
    """
    assert path.exists(filename)
    assert path.isfile(filename)

    k_vectors = np.loadtxt(filename, dtype='d', comments='#')

    assert k_vectors.ndim == 2
    assert k_vectors.shape[1] == 3

    return k_vectors
