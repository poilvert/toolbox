How to use the CCSD(T) extraction code
======================================

This code expects that the user has performed Hartree-Fock (HF) and Coupled-Cluster (CCSD(T)) calculations using a sequence of Dunning type basis sets. Typically, one would want to have performed HF calculations at the Double, Triple and Quadrupole levels (DTQ), and then CCSD(T) calculations at the Double and Triple levels (DT). The first two lowest levels of HF have to match the two levels of CCSD(T).

This code also assumes that the output file names for the CCSD(T) and HF calculations are the same except that the former has the string "CCSDT" in it while the latter substitutes it for "HF".

Generating Gaussian09 input files
=================================

The *gaussian utilities* module contains routines to generate input files for single point energy calculations and geometry optimization calculations.
