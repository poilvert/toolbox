#!/usr/bin/env python

import numpy as np
from math import log, exp


# -- Utilities for CCSD(T) CBS extrapolations
def three_point_extrapolation(energy_l, X_l=[3, 4, 5]):
    """
    Analytical form of the extrapolation formula::

        E(X) = E_CBS + A.exp(-a*X)

    where X is an integer (X=3 for a triple-Zeta basis set, X=4 for a
    quadruple-zeta basis set and X=5 for a quintuple-zeta basis set).

    And where E is usually an energy like the HF energy.
    """

    E1, E2, E3 = energy_l
    assert E3 <= E2 <= E1

    X1, X2, X3 = X_l
    a = -log((E3 - E2) / (E2 - E1))
    alpha = exp(-a)
    B = (E3 - E2) / (alpha**X3 - alpha**X2)
    E_CBS = 1./3.*(E1 - B*alpha**X1 + \
                   E2 - B*alpha**X2 + \
                   E3 - B*alpha**X3)

    assert B >= 0.

    return E_CBS


def two_point_extrapolation(energy_l, X_l=[3, 4]):
    """
    Analytical form of the extrapolation formula::

        F(X) = F_CBS - B.X**(-3)

    where X is an integer (X=3 for a triple-Zeta basis set, X=4 for a
    quadruple-zeta basis set and X=5 for a quintuple-zeta basis set).

    And where F is any relative quantity like the correlation energy.
    """

    F1, F2 = energy_l
    X1, X2 = X_l

    F_CBS = (X2**3*F2 - X1**3*F1)/(X2**3 - X1**3)

    return F_CBS


def get_CCSDT_CBS_energy(HF_energy_l, CCSDT_energy_l, HF_X_l=[3, 4, 5],
        CCSDT_X_l=[3, 4], use_first_two=True):
    """
    This routine allows one to extract the CBS (Complete Basis Set) limit of the
    CCSD(T) energy.
    WARNING: the code uses the first two values of the HF energy to substract
    from the CCSD(T) energies and establish the correlation energies. If this is
    False, the code will use the last two values. Make sure that the X values
    for both HF and CCSD(T) are appropriately set.
    """

    assert len(HF_energy_l) == len(HF_X_l)
    assert len(CCSDT_energy_l) == len(CCSDT_X_l)

    # -- getting the CBS HF energy first
    HF_CBS = three_point_extrapolation(HF_energy_l, X_l=HF_X_l)

    # -- then we generate the series of correlation energies
    if use_first_two:
        CCSDT_minus_HF_energy_l = np.array(CCSDT_energy_l) - \
                                  np.array(HF_energy_l[:2])
    else:
        CCSDT_minus_HF_energy_l = np.array(CCSDT_energy_l) - \
                                  np.array(HF_energy_l[-2:])

    # -- computing the CCSD(T) CBS energy
    CCSDT_minus_HF_CBS = two_point_extrapolation(
                              CCSDT_minus_HF_energy_l, X_l=CCSDT_X_l)
    E_CCSDT = CCSDT_minus_HF_CBS + HF_CBS

    return E_CCSDT
