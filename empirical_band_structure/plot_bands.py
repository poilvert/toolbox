#!/usr/bin/env python

import numpy as np
from compute_band_structure import band_structure
try:
    import pylab
    pylab_present = True
except:
    pylab_present = False


def plot_bands(
    direct_space_basis,
    pseudo_fname,
    k_list_fname,
    G_norm_cutoff,
    n_bands,
    n_electrons_per_cell=None,
    dump_band_energies=False,
    output_filename='band_energies.dat'
    ):

    # -- energy conversion
    HatoeV = 27.21138386

    # -- band structure calculation
    band_energies = band_structure(
        direct_space_basis,
        pseudo_fname,
        k_list_fname,
        G_norm_cutoff,
        n_bands,
        )

    # -- useful parameters for plotting
    n_k_pts, n_bands = band_energies.shape[:2]
    x_grid = np.linspace(0., 1., n_k_pts)

    # -- possibly adjust zero of energy to the Fermi level defined as the energy
    # value halfway between the highest occupied level and the lowest unoccupied
    # level (i.e. midgap if the system is an insulator/semiconductor)
    if n_electrons_per_cell is not None:
        n_electrons = n_k_pts * n_electrons_per_cell
        eigenvalues = np.ravel(band_energies.real)
        sorted_eigenvalues = np.sort(
                np.concatenate(
                            [eigenvalues, eigenvalues]
                          )
                      )
        VBM = sorted_eigenvalues[n_electrons - 1]
        CBM = sorted_eigenvalues[n_electrons]
        fermi_level = 0.5 * (VBM + CBM)
    else:
        fermi_level = 0.
    band_energies -= fermi_level

    # -- dump band energies if asked for it
    if dump_band_energies:
        np.savetxt(output_filename, band_energies)

    # -- plotting the bands
    if pylab_present:
        for idx in range(n_bands):
            pylab.plot(
                    x_grid, HatoeV * band_energies[:, idx].real,
                    '-b', linewidth=3, alpha=0.8
                    )
        pylab.xlabel('Path in k-space', fontsize=16)
        pylab.xticks((0., 1.), ('', ''))
        pylab.ylabel('Band Energies (eV)', fontsize=16)
        pylab.show()
        pylab.clf()
    else:
        print "pylab is missing. Cannot plot bands"
    return
