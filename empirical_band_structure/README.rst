Purpose of these set of programs
================================

The following Python programs will allow one to compute any band structure for any crystal structure as long as the user provides a pseudopotential file containing the Fourier coefficients expressed in the basis of the reciprocal-space vectors.

Quick demo
==========

I provide a little script called `demo_Silicon.py` that should, by its content, explain how to use the different routines in `compute_band_structure.py`, `crystal_utilities.py` and `plot_bands.py` to compute and plot band energies.

Provided pseudopotentials
=========================

The provided pseudopotentials in the `pseudopotentials/` directory are the non-zero Fourier components of the Empirical Pseudopotentials of Silicon, Germanium and Gallium Arsenide as computed from the Pseudopotential *form factors* of Cardona et al. (see table 2.21 of Chapter 2) of **Fundamentals of Semiconductors**, Fourth Edition. I also added a "fake" pseudopotential in order to plot "Free electron" bands in different crystal structures.

Format for pseudopotential files
================================

The set of programs expect pseudopotential files in the following format::

    m  n  p  V_G.real  V_G.imag
    ...

where `m`, `n` and `p` are integer numbers (positive, negative of zero). `V_G.real` and `V_G.imag` are respectively the real part and imaginart part of the Fourier coefficient of the pseudopotential attached to the reciprocal-space vector `G`. As such we assume that the real-space potential is given by::

    V(x) = Sum[V_G * exp(iG.x)]

It is up to the user to make sure that the Fourier coefficients are expressed in *Hartree* atomic units.
