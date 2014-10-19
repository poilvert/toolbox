#!/usr/bin/env python

from math import pi
import numpy as np
from numpy.linalg import det
from numpy.linalg import norm
from numpy.linalg import inv
from scipy.spatial.distance import cdist
from crystal_structures import REPOSITORY_OF_CRYSTAL_STRUCTURES


def get_cubic_canonical_k_path(npts, dump_k_pts=False,
        k_pts_filename='kpoints_cubic.dat'):
    k_points = get_k_space_path(
            'cubic',
            [('Gamma', 'M'), ('M', 'X'), ('X', 'Gamma')],
            npts, dump_k_pts=dump_k_pts, k_pts_filename=k_pts_filename
        )
    return k_points


def get_fcc_canonical_k_path(npts, dump_k_pts=False,
        k_pts_filename='kpoints_fcc.dat'):
    k_points = get_k_space_path(
            'fcc',
            [('L', 'Gamma'), ('Gamma', 'X'), ('X1', 'Gamma')],
            npts, dump_k_pts=dump_k_pts, k_pts_filename=k_pts_filename
        )
    return k_points


def get_bcc_canonical_k_path(npts, dump_k_pts=False,
        k_pts_filename='kpoints_bcc.dat'):
    k_points = get_k_space_path(
            'bcc',
            [('H', 'Gamma'), ('Gamma', 'N'), ('N', 'P')],
            npts, dump_k_pts=dump_k_pts, k_pts_filename=k_pts_filename
        )
    return k_points


def get_hexagonal_canonical_k_path(npts, dump_k_pts=False,
        k_pts_filename='kpoints_hcp.dat'):
    k_points = get_k_space_path(
            'hexagonal',
            [('Gamma', 'K'), ('K', 'H'), ('H', 'A'), ('A', 'Gamma'),
             ('Gamma', 'M'), ('M', 'L'), ('L', 'A')],
            npts, dump_k_pts=dump_k_pts, k_pts_filename=k_pts_filename
        )
    return k_points


def get_k_space_path(crystal_string, list_of_k_pt_couples, n_ks,
                     lattice_constants=None, dump_k_pts=False,
                     k_pts_filename='k_points.dat'):
    """
    This little routine will produce *normalized* paths in k-space suited for
    band structure calculations. It reads a string that gives the crystal
    structure of the system, and a list of k-point *couples* that indicates the
    lines along which one computes k points.

    Note
    ====

    This code uses the conventions described in the repository of crystal
    structures in the file ``crystal_structures`` for its high-symmetry points.

    Moreover the code produces paths that are normalized in terms of Euclidean
    distances in the Brillouin Zone (i.e. longer lines gets more k points).

    Finally **all** k point coordinates are given in the basis of irreducible
    lattice vectors in *reciprocal* space for the given crystal structure.
    """

    def path_length(k_couple, k_pts_dict, adot):
        k0 = np.array(k_pts_dict[k_couple[0]])
        k1 = np.array(k_pts_dict[k_couple[1]])
        return adot(k1-k0, k1-k0)

    def line_points(k_couple, k_pts_dict, n_pts):
        assert n_pts > 1
        k0 = np.array(k_pts_dict[k_couple[0]])
        k1 = np.array(k_pts_dict[k_couple[1]])
        ks = [k0+(k1-k0)*i/(n_pts-1) for i in range(n_pts)]
        return np.array(ks)

    if lattice_constants is None:
        if crystal_string == 'hexagonal':
            constants = {'a': 1., 'c': 1.}
        else:
            constants = {'a': 1.}
    else:
        constants = lattice_constants

    # -- Direct-Space (i.e. Real Space) basis of elementary vectors
    ds_basis = get_direct_space_elementary_vectors(
                        crystal_string,
                        constants
                    )

    # -- Reciprocal-Space basis of elementary vectors
    rs_basis = get_reciprocal_space_elementary_vectors(ds_basis)

    # -- Dot product function to compute distances and angles between vectors in
    # reciprocal-space
    mydot = get_Euclidean_dot_product(rs_basis)

    # -- Coordinates of the High-Symmetry points
    high_sym_pts = \
        REPOSITORY_OF_CRYSTAL_STRUCTURES[crystal_string]['high_symmetry_k_points']

    # -- Computing the total Euclidean length of the selected k-path
    total_path_length = 0.
    for k_couple in list_of_k_pt_couples:
        total_path_length += path_length(k_couple, high_sym_pts, mydot)

    # -- Finally computing the list of k vectors on the selected path
    k_points = []
    for k_couple in list_of_k_pt_couples:
        n_pts = int(n_ks * path_length(k_couple, high_sym_pts, mydot) \
                    / total_path_length)
        k_points += [line_points(k_couple, high_sym_pts, n_pts)]
    k_points = np.concatenate(k_points)

    if dump_k_pts:
        np.savetxt(k_pts_filename, k_points, fmt='%10.6f')
    else:
        return k_points


def get_direct_space_elementary_vectors(crystal_string, lattice_constants):
    """
    Reads a string like 'fcc' and a dictionnary containing the values for the
    conventional lattice parameters (e.g. 'a' in bcc, cubic and fcc, and 'a' and
    'c' for hexagonal). It then returns the basis of irreducible direct space
    lattice vectors.

    Note
    ====

    This code uses the conventions described in the repository of crystal
    structures in the file ``crystal_structures``
    """
    acceptable_crystal_strings = REPOSITORY_OF_CRYSTAL_STRUCTURES.keys()
    assert crystal_string in acceptable_crystal_strings

    crystal_info = REPOSITORY_OF_CRYSTAL_STRUCTURES[crystal_string]
    irreducible_basis_dict = crystal_info['irreducible_basis']
    basis = np.array([irreducible_basis_dict['a1'],
                      irreducible_basis_dict['a2'],
                      irreducible_basis_dict['a3']])
    lattice_constants_names = lattice_constants.keys()
    lattice_constants_names_ref = crystal_info['lattice_constant_names']
    for name in lattice_constants_names:
        assert name in lattice_constants_names_ref

    if crystal_string == 'hexagonal':
        basis[:2, :] *= lattice_constants['a']
        basis[2, :] *= lattice_constants['c']
    else:
        basis *= lattice_constants['a']

    return basis


def get_elementary_cell_volume(direct_space_elementary_vectors):
    """
    Given a matrix that gives the cartesian coordinates of all direct space
    elementary vectors, this program computes the volume of the elementary cell
    spanned by these elementary vectors.
    """
    a_mat = np.array(direct_space_elementary_vectors).astype('d')
    assert a_mat.ndim == 2
    assert a_mat.shape == (3, 3)

    reference_volume = np.dot(a_mat[0, :], np.cross(a_mat[1, :], a_mat[2, :]))
    check_volume = det(a_mat)

    assert abs(reference_volume - check_volume) < 1e-6

    return check_volume


def get_reciprocal_space_elementary_vectors(direct_space_elementary_vectors):
    """
    Given a matrix that gives the cartesian coordinates of all direct space
    elementary vectors, this program computes the basis of reciprocal space
    elementary vectors.

    Input
    =====
    direct space elementary vectors
        2D array of shape (3, 3)

    Returns
    =======
    reciprocal space elementary vectors
        2D array of shape (3, 3)

    Conventions
    ===========
    the format of the input array is assumed to be the following
        a1x  a1y  a1z
        a2x  a2y  a2z
        a3x  a3y  a3z
    if the matrix elements are given in Bohr units, then the returned matrix is
    given in units of inverse Bohr.
    """
    a_mat = np.array(direct_space_elementary_vectors).astype('d')
    assert a_mat.ndim == 2
    assert a_mat.shape == (3, 3)

    b_mat = np.zeros((3, 3), dtype='d')
    V = get_elementary_cell_volume(a_mat)
    prefactor = 2 * pi * (1. / V)
    b_mat[0, :] = prefactor * np.cross(a_mat[1, :], a_mat[2, :])
    b_mat[1, :] = prefactor * np.cross(a_mat[2, :], a_mat[0, :])
    b_mat[2, :] = prefactor * np.cross(a_mat[0, :], a_mat[1, :])

    assert np.abs(
            np.dot(a_mat, b_mat.T) - 2*pi*np.eye(3, dtype='d')
            ).sum() < 1e-6

    return b_mat


def get_Euclidean_dot_product(basis_matrix):
    """
    This routine will produce a "dot product" function that will allow a user to
    compute Euclidean dot products between two vectors, which coordinates are
    given in the basis of vectors from the input tensor.

    Input
    =====
    basis
        2D array of shape (N, N)

    Returns
    =======
    my_dot
        a function that allows the computation of Euclidean dot products of two
        vectors, which coordinates are given in the basis of vectors given by
        the input tensor

    Note
    ====
    if the input tensor elements are given in a defined unit, then the dot
    product values will have the dimensionality of a square unit.
    """
    mat = np.array(basis_matrix).astype('d')
    assert mat.ndim == 2
    assert mat.shape[0] == mat.shape[1]

    metric_tensor = np.dot(mat, mat.T)
    def my_dot(x, y):
        my_x = np.array(x).astype('d')
        my_y = np.array(y).astype('d')
        assert my_x.ndim == 1
        assert my_y.ndim == 1
        assert my_x.size == my_y.size
        assert my_x.size == metric_tensor.shape[0]
        return np.dot(my_x, np.dot(metric_tensor, my_y))
    return my_dot


def get_ordered_stars(basis, Euclidean_norm_cutoff,
        augment_range=0):
    """
    This program will take a basis of elementary vectors, to build a basis of
    lattice vectors (i.e. vectors for which the coordinates in the basis of
    elementary vectors are integers) of increasing Euclidean norm. The maximal
    Euclidean norm allowed is given by the provided cutoff value.
    """
    a_mat = np.array(basis).astype('d')
    R = float(Euclidean_norm_cutoff)

    assert R >= 0.
    assert a_mat.ndim == 2
    assert a_mat.shape == (3, 3)

    mydot = get_Euclidean_dot_product(a_mat)
    Rsq = R ** 2

    # -- computes the area of the parallelogram spanned by vectors v1 and v2
    def get_area(v1, v2):
        v3 = np.cross(v1, v2)
        return norm(v3)

    # -- here we are looking for a tight regular mesh of lattice vectors that
    # contain the sphere of radius R
    V = get_elementary_cell_volume(a_mat)
    A0 = get_area(a_mat[1, :], a_mat[2, :])
    A1 = get_area(a_mat[2, :], a_mat[0, :])
    A2 = get_area(a_mat[0, :], a_mat[1, :])
    R0, R1, R2 = V / A0, V / A1, V / A2
    m0, m1, m2 = int(R / R0) + 1 + augment_range, \
                 int(R / R1) + 1 + augment_range, \
                 int(R / R2) + 1 + augment_range
    l0, l1, l2 = range(-m0, m0+1), range(-m1, m1+1), range(-m2, m2+1)

    raw_vectors = np.array([
            (x, y, z) for x in l0 for y in l1 for z in l2
            if (mydot([x, y, z], [x, y, z]) <= Rsq)
            ]).astype('int')

    # -- sorting the vectors by increasing Euclidean norm
    norms_sq = [mydot(vec, vec) for vec in raw_vectors]
    sorted_idx = np.argsort(norms_sq)
    sorted_vectors = raw_vectors[sorted_idx]

    return sorted_vectors


def remove_duplicate_atoms(atoms_coords, eps=1e-3):
    """
    Identifies duplicate rows (up to eps) and removes them.
    """

    assert atoms_coords.ndim == 2
    assert atoms_coords.shape[1] == 3

    distances = cdist(atoms_coords, atoms_coords)
    mask = (distances <= eps)
    row_group_idx_l = [np.sort(np.flatnonzero(row)) for row in mask]
    idx_to_keep = np.unique(
                [group[0] for group in row_group_idx_l]
            )

    return atoms_coords[idx_to_keep]


def cartoint(basis, atom_basis, eps=1e-3):
    """
    Transforms atomic coordinates from cartesian to intrinsic. Intrinsic means
    coordinates in the basis of row vectors from ``basis``.
    """

    assert basis.shape == (3, 3)
    assert atom_basis.ndim == 2
    assert atom_basis.shape[1] == 3
    assert det(basis) > eps

    P = inv(basis.T)

    return np.dot(P, atom_basis.T).T


def inttocar(basis, atom_basis, eps=1e-3):
    """
    Transforms atomic coordinates from intrinsic to cartesian. Intrinsic means
    coordinates in the basis of row vectors from ``basis``.
    """

    assert basis.shape == (3, 3)
    assert atom_basis.ndim == 2
    assert atom_basis.shape[1] == 3
    assert det(basis) > eps

    P = basis.T

    return np.dot(P, atom_basis.T).T


def transport_to_home_cell(basis, atom_basis,
        coord_type='cartesian', decimals=5):
    """
    Given a basis of atomic coordinates and a basis of basic translation
    vectors, this code will "transport" all the atoms to the "home cell". The
    "home cell" is defined by imposing that the intrinsic atomic coordinates be
    between 0 (included) and 1 (excluded).

    Warning
    =======
    This function uses a rounding of the atomic coordinates to make sure that
    points that should be close together, won't be "transported" to different
    locations in the home cell when the integer part has been removed.
    """

    assert basis.shape == (3, 3)
    assert atom_basis.ndim == 2
    assert atom_basis.shape[1] == 3
    assert coord_type in ['cartesian', 'intrinsic']

    if coord_type == 'cartesian':
        atom_basis = cartoint(basis, atom_basis)

    atom_basis = np.around(atom_basis, decimals=decimals)

    home_cell = atom_basis - np.floor(atom_basis)

    if coord_type == 'cartesian':
        home_cell = inttocar(basis, home_cell)

    return home_cell


def standardize_atomic_positions(basis, atom_basis, coord_type='cartesian'):
    """
    This function does the following::

        1. transports the atoms to the home cell
        2. de-duplicate atoms if necessary
    """

    assert basis.shape == (3, 3)
    assert atom_basis.ndim == 2
    assert atom_basis.shape[1] == 3
    assert coord_type in ['cartesian', 'intrinsic']

    atom_basis = transport_to_home_cell(
            basis, atom_basis, coord_type=coord_type
            )

    atom_basis = remove_duplicate_atoms(atom_basis)

    return atom_basis


def find_neighbor_shells(
        basis, atom_basis, R,
        coord_type='cartesian', decimals=3
        ):
    """
    Find all the nearest neighbor shells within a distance R for all the atoms
    in the "standardized" basis.

    Warning
    =======
    This code returns the nearest neighbor coordinates and the standardized
    basis vector coordinates in a cartesian basis.
    """

    assert basis.shape == (3, 3)
    assert atom_basis.ndim == 2
    assert atom_basis.shape[1] == 3
    assert R >= 0.
    assert coord_type in ['cartesian', 'intrinsic']

    # -- utility
    def _nn_group_idx_l(distance_vec, R0):

        assert distance_vec.ndim == 1

        argsorted_idx = distance_vec.argsort()
        distance_vec = distance_vec[argsorted_idx]

        mask = (distance_vec <= R0)

        argsorted_idx = argsorted_idx[mask]
        distance_vec = distance_vec[mask]

        u_dist, u_idx = np.unique(distance_vec, return_index=True)
        u_idx = np.concatenate([u_idx, [argsorted_idx.size]])

        group_idx_l = []
        for start_idx, stop_idx in zip(u_idx[:-1], u_idx[1:]):
            group_idx_l += [argsorted_idx[start_idx:stop_idx]]

        return group_idx_l, u_dist

    # -- core of the code
    a1, a2, a3 = basis

    atom_basis = standardize_atomic_positions(basis, atom_basis,
                    coord_type=coord_type)

    # -- transform intrinsic coordinates into cartesian to compute Euclidean
    # distances
    if coord_type == 'intrinsic':
        atom_basis = inttocar(basis, atom_basis)

    atom_distances = np.array([norm(vec) for vec in atom_basis])
    Rc = R + atom_distances.max()
    cells = get_ordered_stars(basis, Rc)

    assert cells.ndim == 2
    assert cells.shape[1] == 3

    all_atoms = np.concatenate(
            [
                (atom_basis + x*a1 + y*a2 + z*a3)
                for (x, y, z) in cells
            ]
        )

    assert all_atoms.ndim == 2
    assert all_atoms.shape[1] == 3

    n_basis_atoms = atom_basis.shape[0]
    n_cells = cells.shape[0]

    cell_indices = np.repeat(cells, n_basis_atoms, axis=0)
    atom_indices = np.tile(np.arange(n_basis_atoms), n_cells)

    assert cell_indices.shape == all_atoms.shape
    assert atom_indices.ndim == 1
    assert atom_indices.shape[0] == cell_indices.shape[0]

    distance_matrix = cdist(all_atoms, all_atoms)
    distance_matrix = np.around(distance_matrix, decimals=decimals)

    all_neighbors = []
    for i in xrange(n_basis_atoms):
        atom_neighbors = {}
        atom_group_idx_l, atom_dist_l = _nn_group_idx_l(
                distance_matrix[..., i], R
                )
        assert len(atom_group_idx_l) == len(atom_dist_l)
        for j in xrange(len(atom_group_idx_l)):
            key = '%i' % j
            group_idx = atom_group_idx_l[j]
            dist = atom_dist_l[j]
            atom_neighbors[key] = {
                    "neighbors_cart_coordinates": all_atoms[group_idx],
                    "neighbors_cell_coordinates": cell_indices[group_idx],
                    "neighbors_basis_indices": atom_indices[group_idx],
                    "distance": dist,
                }
        all_neighbors += [atom_neighbors]

    return all_neighbors


def Wigner_Seitz_cell(
        basis, atom_basis,
        coord_type='cartesian', decimals=5
        ):

    assert basis.shape == (3, 3)
    assert atom_basis.ndim == 2
    assert atom_basis.shape[1] == 3
    assert coord_type in ['cartesian', 'intrinsic']

    atom_basis = standardize_atomic_positions(basis, atom_basis,
                    coord_type=coord_type)
    n_atoms = len(atom_basis)

    if coord_type == 'intrinsic':
        atom_basis = inttocar(basis, atom_basis)

    # -- compute nearest lattice vectors
    a1, a2, a3 = basis
    nearest_latt_vecs = np.array(
                [
                    i*a1 + j*a2 + k*a3
                        for i in xrange(-1,2)
                        for j in xrange(-1,2)
                        for k in xrange(-1,2)
                        if (i, j, k) != (0, 0, 0)
                ]
            )

    # -- compute atoms and translationally equivalent ones for all 8 unit cells
    # surrounding the origin of the "standardized home cell"
    cells = [
                (x, y, z)
                    for x in [-1, 0]
                    for y in [-1, 0]
                    for z in [-1, 0]
            ]

    all_atoms = np.concatenate(
            [
                (atom_basis + x*a1 + y*a2 + z*a3)
                for (x, y, z) in cells
            ]
        )

    # -- normalizing the nearest lattice vectors
    nearest_latt_vecs /= (nearest_latt_vecs ** 2).sum(1)[:, None]
    n_nearest_latt_vecs = len(nearest_latt_vecs)

    # -- matrix of projections onto the nearest lattice vectors
    proj_matrix = np.around(
            np.dot(all_atoms, nearest_latt_vecs.T),
            decimals=decimals
            )

    # -- finding out what are the atoms within the Wigner-Seitz cell
    mask = (proj_matrix >= -0.5) * (proj_matrix <= 0.5)
    atoms_in_ws = (mask.sum(1) == n_nearest_latt_vecs)

    # -- extraction of the final atoms in the WS cell
    selected_atom_idx = []
    for i in xrange(n_atoms):
        atom_mask_vals = atoms_in_ws[i::n_atoms]
        selected_atom_idx += [
                i + np.flatnonzero(atom_mask_vals)[0] * n_atoms
                ]
    final_atom_basis = all_atoms[selected_atom_idx]
    assert final_atom_basis.shape == atom_basis.shape

    sorted_distance_idx = np.argsort([norm(vec) for vec in final_atom_basis])

    if coord_type == 'intrinsic':
        final_atom_basis = cartoint(basis, final_atom_basis)

    return final_atom_basis[sorted_distance_idx]
