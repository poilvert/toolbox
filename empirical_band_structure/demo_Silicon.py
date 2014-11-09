#!/usr/bin/env python

from plot_bands import plot_bands
from crystal_utilities import get_direct_space_elementary_vectors
from crystal_utilities import get_fcc_canonical_k_path

AngtoBohr = 1. / 0.529177

crystal_type = 'fcc'
pseudopotential = 'pseudopotentials/Silicon.dat'
lattice_parameter = 5.43
G_norm_cutoff = 2.8
n_kpts = 300
k_pts_filename = 'demo_k_path.dat'
n_bands = 8
n_electrons_per_cell = 8

# -- Basis of FCC elementary direct-space vectors
basis = get_direct_space_elementary_vectors(
            crystal_type,
            {'a': AngtoBohr * lattice_parameter}
        )

# -- Canonical k-space path for band structure calculations in FCC lattices
k_path = get_fcc_canonical_k_path(n_kpts, dump_k_pts=True,
        k_pts_filename=k_pts_filename)

plot_bands(
        basis,
        pseudopotential,
        k_pts_filename,
        G_norm_cutoff,
        n_bands,
        n_electrons_per_cell=n_electrons_per_cell,
        )
