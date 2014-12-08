Large Scale Quantum Conductance Package (LSQC)
==============================================

Purpose
-------

This set of Python programs is to be used as a post-processing utility to the Wannier90_ code.
The intent is to provide the user with a high-level abstraction for manipulating Hamiltonian
matrices provided by the Wannier90_ code. In turn, those matrices can be used to build a very
large Hamiltonian matrix and hence compute the Quantum Conductance of large scale systems.

For more information, please look at the documentation.

Documentation
-------------

You can access the HTML documentation (built with Sphinx_) right here_.

Programs
--------

You can freely access the programs and utilities by looking into the **programs/** and **utilities/**
directories.

Recquirements
-------------

Python >= 2.5
Numpy_ >= 1.0

And for one of the programs:
Mayavi2_

Credits
-------

I would like to thank **Nilton Volpato** for providing such a compact and easy-to-use module for
managing progress bars in Python.

.. _here: http://poilvert.github.com/LSQC/index.html
.. _Wannier90: http://www.wannier.org
.. _Sphinx: http://sphinx.pocoo.org/
.. _Numpy: http://numpy.scipy.org/
.. _Mayavi2: http://code.enthought.com/projects/mayavi/#Mayavi2

