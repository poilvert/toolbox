#!/usr/bin/env python

import sys
from os import path
import optparse
import numpy as np
from extrapolation import get_CCSDT_CBS_energy


__all__ = ['CCSDT_CBS_energy']


# -- how and where to find the energy value
HF_KEYWORD = 'SCF Done'
HF_E_POSITION = 4
CCSDT_KEYWORD = 'CCSD(T)= '
CCSDT_E_POSITION = 1

# -- how to find which basis sets were used
BASIS_SET_KEYWORD = 'Standard basis:'
BASIS_SET_POSITION = 2


# -- Internal utilities
def _convert_energy_string_to_float(energy_string):
    energy_string = str(energy_string)
    energy_string = energy_string.upper()
    if 'D' in energy_string:
        if 'D-' in energy_string:
            items = energy_string.split('D-')
            sign = -1
        elif 'D+' in energy_string:
            items = energy_string.split('D+')
            sign = 1
        else:
            items = energy_string.split('D')
            sign = 1
        assert len(items) == 2
        significant, exponent = float(items[0]), int(items[1])
        energy = significant * 10 ** (sign * exponent)
    else:
        energy = float(energy_string)

    return energy


def _get_energies(
        lines,
        key=HF_KEYWORD,
        energy_position=HF_E_POSITION,
        ):
    energy_values = []
    for line in lines:
        line = line.strip()
        if key in line:
            energy_string = line.split()[energy_position]
            energy = _convert_energy_string_to_float(energy_string)
            energy_values += [energy]

    return energy_values


def _get_basis_sets(
        lines,
        key=BASIS_SET_KEYWORD,
        position=BASIS_SET_POSITION,
        ):
    basis_set_l = []
    for line in lines:
        line = line.strip()
        if key in line:
            basis_set = line.split()[position]
            basis_set_l += [basis_set]

    return basis_set_l


def _get_ordered_basis_sets(
        basis_set_string_l,
        top_N=2
        ):
    assert len(basis_set_string_l) >= top_N
    assert top_N >= 1
    unsorted_zetas = []
    for basis_set_string in basis_set_string_l:
        basis_set_string = basis_set_string.upper()
        if 'DZ' in basis_set_string:
            unsorted_zetas += [2]
        elif 'TZ' in basis_set_string:
            unsorted_zetas += [3]
        elif 'QZ' in basis_set_string:
            unsorted_zetas += [4]
        elif '5Z' in basis_set_string:
            unsorted_zetas += [5]
        elif '6Z' in basis_set_string:
            unsorted_zetas += [6]
        else:
            print 'cannot extract Zeta value for basis set %s' % (
                                basis_set_string,)
            sys.exit(1)

    assert len(unsorted_zetas) >= top_N
    sorted_idx = np.argsort(unsorted_zetas)
    sorted_zetas = np.array(unsorted_zetas)[sorted_idx]

    return sorted_idx[-top_N:], sorted_zetas[-top_N:]


# -- Actual utilities
def CCSDT_CBS_energy(
        CCSDT_fname,
        CCSDT_keyword='CCSDT',
        HF_keyword='HF'
        ):

    assert path.exists(CCSDT_fname)
    assert path.isfile(CCSDT_fname)

    print 'looking into: %s' % CCSDT_fname

    dirname = path.dirname(CCSDT_fname)
    basename = path.basename(CCSDT_fname)

    HF_fname = path.join(dirname, basename.replace(CCSDT_keyword, HF_keyword))

    print 'looking into: %s' % HF_fname

    assert path.exists(HF_fname)
    assert path.isfile(HF_fname)

    with open(CCSDT_fname, 'r') as CCSDT_fin:
        CCSDT_lines = CCSDT_fin.readlines()
    with open(HF_fname, 'r') as HF_fin:
        HF_lines = HF_fin.readlines()

    CCSDT_basis_set_l = _get_basis_sets(CCSDT_lines)
    HF_basis_set_l = _get_basis_sets(HF_lines)

    CCSDT_energy_l = _get_energies(
                            CCSDT_lines,
                            key=CCSDT_KEYWORD,
                            energy_position=CCSDT_E_POSITION,
                            )
    HF_energy_l = _get_energies(
                            HF_lines,
                            key=HF_KEYWORD,
                            energy_position=HF_E_POSITION,
                            )

    assert len(CCSDT_basis_set_l) == len(CCSDT_energy_l)
    assert len(HF_basis_set_l) == len(HF_energy_l)

    assert len(CCSDT_basis_set_l) >= 2
    assert len(HF_basis_set_l) >= 3

    CCSDT_order, CCSDT_zetas = _get_ordered_basis_sets(
                                                CCSDT_basis_set_l,
                                                top_N=2
                                                )
    HF_order, HF_zetas = _get_ordered_basis_sets(
                                                HF_basis_set_l,
                                                top_N=3
                                                )

    CCSDT_sorted_energies = np.array(CCSDT_energy_l)[CCSDT_order]
    HF_sorted_energies = np.array(HF_energy_l)[HF_order]

    assert (CCSDT_zetas == HF_zetas[:2]).all()

    CCSDT_CBS_energy = get_CCSDT_CBS_energy(
                            HF_sorted_energies,
                            CCSDT_sorted_energies,
                            HF_X_l=HF_zetas,
                            CCSDT_X_l=CCSDT_zetas,
                            use_first_two=True
                            )

    return CCSDT_CBS_energy


# -------------------
# -- CLI for the code
# -------------------

def main():

    parser = optparse.OptionParser()

    parser.add_option('-i','--input',
                    dest="input_fname",
                    type="str",
                    help="path to the Gaussian09 CCSD(T) output file")
    parser.add_option('-o','--output',
                    dest="output_fname",
                    type="str",
                    default="if_default",
                    help="output fname with extrapolated CCSD(T) energy")

    options, remainder = parser.parse_args()

    if len(sys.argv[1:])==0:
        parser.print_help()
        sys.exit(0)

    if options.output_fname == "if_default":
        output_fname = options.input_fname + '.energy'
    else:
        output_fname = options.output_fname

    energy = CCSDT_CBS_energy(
            options.input_fname,
            )

    with open(output_fname, 'w') as fout:
        fout.write('# Total energy (Ha) / extrapolated CCSD(T)/CBS energy\n')
        fout.write('%12.7f\n' % energy)

    return


if __name__ == "__main__":
    main()
