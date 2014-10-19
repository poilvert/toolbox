#!/usr/bin/env python

import numpy as np
from numpy.linalg import eigvalsh
from io_utilities import read_pseudo
from io_utilities import read_kpoints
from crystal_utilities import get_reciprocal_space_elementary_vectors
from crystal_utilities import get_Euclidean_dot_product
from crystal_utilities import get_ordered_stars


# -- routines used in main program
def build_potential_matrix(
    G_basis,
    pseudo_Gs,
    pseudo_coefficients,
    ):

    # -- checks on inputs
    assert G_basis.ndim == 2
    assert G_basis.shape[1] == 3
    assert pseudo_Gs.ndim == 2
    assert pseudo_Gs.shape[1] == 3
    assert pseudo_coefficients.ndim == 1
    assert pseudo_coefficients.dtype == 'complex'
    assert pseudo_Gs.shape[0] == pseudo_coefficients.size
    assert G_basis.dtype == pseudo_Gs.dtype

    (n_Gs, _) = G_basis.shape
    potential_matrix = np.zeros((n_Gs, n_Gs), dtype='complex')

    dico = {}
    for G, coeff in zip(pseudo_Gs, pseudo_coefficients):
        key = '%i %i %i' % (G[0], G[1], G[2])
        dico[key] = coeff

    # -- strictly upper triangular part of the potential matrix
    Gi_minus_Gj = G_basis[:, None, :] - G_basis[None, :, :]
    for i in xrange(n_Gs-1):
        for j in xrange(i+1, n_Gs):
            G = Gi_minus_Gj[i, j]
            key = '%i %i %i' % (G[0], G[1], G[2])
            if key in dico:
                potential_matrix[i, j] = dico[key]

    # -- now we build the lower triangular part by symmetry
    potential_matrix += np.conjugate(potential_matrix.T)

    # -- possibly the diagonal
    key0 = '%i %i %i' % (0, 0, 0)
    if key0 in dico:
        for i in xrange(n_Gs):
            potential_matrix[i, i] = dico[key0]

    return potential_matrix


def build_kinetic_matrix(k_vec, G_basis, G_space_dot):

    assert k_vec.ndim == 1
    assert G_basis.ndim == 2
    assert G_basis.shape[1] == 3

    k_sq_norms = 0.5 * np.array(
            [G_space_dot(k_vec+k, k_vec+k) for k in G_basis]
            ).astype('complex')

    return np.diag(k_sq_norms)


# -- computing the band structure
def band_structure(
    direct_space_basis,
    pseudo_fname,
    k_list_fname,
    G_norm_cutoff,
    n_bands,
    ):

    # -- get all necessary objects needed to compute the band structure
    pseudo_Gs, pseudo_coefficients = read_pseudo(pseudo_fname)
    k_vectors = read_kpoints(k_list_fname)
    reciprocal_space_basis = \
        get_reciprocal_space_elementary_vectors(direct_space_basis)
    rs_dot = get_Euclidean_dot_product(reciprocal_space_basis)

    assert pseudo_Gs.shape[0] == pseudo_coefficients.size
    assert k_vectors.shape[1] == 3

    # -- Basis of G vectors used for the matrix diagonalization
    G_basis = get_ordered_stars(reciprocal_space_basis, G_norm_cutoff)

    # -- generate potential matrix, which is k-independent
    potential_matrix = build_potential_matrix(
            G_basis, pseudo_Gs, pseudo_coefficients
            )

    # -- loop over the k points to compute the band energies
    band_energies = []
    for idx, k_vec in enumerate(k_vectors):
        kinetic_matrix = build_kinetic_matrix(k_vec, G_basis, rs_dot)
        hamiltonian = kinetic_matrix + potential_matrix
        eigvals = np.sort(eigvalsh(hamiltonian))[:n_bands]
        band_energies += [eigvals]

    # -- simple check to see if the eigenvalues are real
    band_energies = np.array(band_energies)
    assert np.abs(band_energies.imag).sum() < 1e-3

    return band_energies
