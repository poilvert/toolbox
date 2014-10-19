#!/usr/bin/env python

"""
Crystal structure information obtained from the following paper::

    http://arxiv.org/pdf/1004.2974.pdf
"""

from math import sqrt


REPOSITORY_OF_CRYSTAL_STRUCTURES = {
    'cubic': {
        'lattice_constant_names': ['a'],
        'irreducible_basis': {
            'a1': [1., 0., 0.],
            'a2': [0., 1., 0.],
            'a3': [0., 0., 1.]
        },
        'high_symmetry_k_points': {
            'Gamma': [0., 0., 0.],
            'M': [0.5, 0.5, 0.],
            'R': [0.5, 0.5, 0.5],
            'X': [0., 0.5, 0.,]
        }
    },
    'fcc': {
        'lattice_constant_names': ['a'],
        'irreducible_basis': {
            'a1': [0., 0.5, 0.5],
            'a2': [0.5, 0., 0.5],
            'a3': [0.5, 0.5, 0.]
        },
        'high_symmetry_k_points': {
            'Gamma': [0., 0., 0.],
            'K': [3./8., 3./8., 3./4.],
            'L': [0.5, 0.5, 0.5],
            'U': [5./8., 1./4., 5./8.],
            'W': [0.5, 1./4., 3./4.],
            'X': [0.5, 0., 0.5],
            'X1': [0.5, 0.5, 1.]
        }
    },
    'bcc': {
        'lattice_constant_names': ['a'],
        'irreducible_basis': {
            'a1': [-0.5, 0.5, 0.5],
            'a2': [0.5, -0.5, 0.5],
            'a3': [0.5, 0.5, -0.5]
        },
        'high_symmetry_k_points': {
            'Gamma': [0., 0., 0.],
            'H': [0.5, -0.5, 0.5],
            'P': [1./4., 1./4., 1./4.],
            'N': [0., 0., 0.5]
        }
    },
    'hexagonal': {
        'lattice_constant_names': ['a', 'c'],
        'irreducible_basis': {
            'a1': [0.5, -sqrt(3.)/2., 0.],
            'a2': [0.5, sqrt(3.)/2., 0.],
            'a3': [0., 0., 1.]
        },
        'high_symmetry_k_points': {
            'Gamma': [0., 0., 0.],
            'A': [0., 0., 0.5],
            'H': [1./3., 1./3., 1./2.],
            'K': [1./3., 1./3., 0.],
            'L': [0.5, 0., 0.5],
            'M': [0.5, 0., 0.]
        }
    }
}
