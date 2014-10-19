#!/usr/bin/env python

from os import path


# -- Default parameters
DEFAULT_TOTAL_ENERGY_PATTERN = 'total energy ='
DEFAULT_KINETIC_ENERGY_PATTERN = 'kinetic energy ='
DEFAULT_ELECTROSTATIC_ENERGY_PATTERN = 'electrostatic energy ='
DEFAULT_XC_ENERGY_PATTERN = 'exchange-correlation energy ='
DEFAULT_ENERGY_PATTERNS = {
        'total_energy': DEFAULT_TOTAL_ENERGY_PATTERN,
        'kinetic_energy': DEFAULT_KINETIC_ENERGY_PATTERN,
        'electrostatic_energy': DEFAULT_ELECTROSTATIC_ENERGY_PATTERN,
        'xc_energy': DEFAULT_XC_ENERGY_PATTERN,
        }
DEFAULT_MAJ_SPIN_PATTERNS = ['Eigenvalues (eV),', ', spin =  1']
DEFAULT_MIN_SPIN_PATTERNS = ['Eigenvalues (eV),', ', spin =  2']
DEFAULT_HOMO_PATTERN = 'HOMO Eigenvalue (eV)'
DEFAULT_LUMO_PATTERN = 'LUMO Eigenvalue (eV)'
DEFAULT_LAMBDA_PATTERNS = ['lambda', 'spin =']
DEFAULT_CENTER_PATTERNS = ['Charge', 'Centers xyz (Bohr)', 'Spreads', 'spin =']


# -----------------------------
# -- Internal utility functions
# -----------------------------
def get_energy_from_line(line):
    """Given a line string, the code will extract the first object from that
    line that is understandable as a floating point number.
    """
    energy = 0.
    for item in line.split():
        try:
            energy = float(item)
            break
        except:
            pass
    return energy


def find_matching_lines_idx_for_patterns(lines, patterns):
    """Given a list of line strings, patterns are searched within these lines.
    For all the lines containing *all* of the patterns, we keep track of their
    index (i.e. the index of that line in the list).
    """
    found_indices = []
    for idx, line in enumerate(lines):
        select_line = True
        for pattern in patterns:
            if pattern not in line:
                select_line = False
                break
        if select_line:
            found_indices += [idx]
    return found_indices


def get_orbital_energies_for_idx(lines, idx):
    """Extracts all the orbital energies after the ``idx``-th line of ``lines``
    and stops extraction as soon as a line does not contain any floating point
    numbers.

    Warning
    =======

    This code assumes that the orbital energies only start at the ``idx+2``-th
    line, as expected from a CP output file.
    """

    assert len(lines) > idx

    eigenvalues = []
    can_go_to_next_line = True
    counter = idx + 1
    while can_go_to_next_line:
        try:
            objects = lines[counter+1].split()
            eigenvalues += [float(eig) for eig in objects]
            counter += 1
        except:
            can_go_to_next_line = False
    return eigenvalues


def get_homo_level(lines, homo_pattern=DEFAULT_HOMO_PATTERN):
    """Simply extracts the HOMO level from a CP output file
    """
    homo = 0.
    new_lines = lines[::-1]
    for idx, line in enumerate(new_lines):
        if homo_pattern in line:
            homo = float(new_lines[idx-2].split()[0])
            break
    return homo


def get_lumo_level(lines, lumo_pattern=DEFAULT_LUMO_PATTERN):
    """Simply extracts the LUMO level from a CP output file
    """
    lumo = 0.
    new_lines = lines[::-1]
    for idx, line in enumerate(new_lines):
        if lumo_pattern in line:
            lumo = float(new_lines[idx-2].split()[0])
            break
    return lumo


def get_orbital_energies(lines):
    """this code extracts all the orbital energies for all the spin states
    (minority spins, majority spins) at all k-points.
    """

    assert len(lines) > 0

    maj_spin_patterns = DEFAULT_MAJ_SPIN_PATTERNS
    min_spin_patterns = DEFAULT_MIN_SPIN_PATTERNS
    maj_spin_line_idx = find_matching_lines_idx_for_patterns(
            lines,
            maj_spin_patterns
            )
    min_spin_line_idx = find_matching_lines_idx_for_patterns(
            lines,
            min_spin_patterns
            )

    maj_spin_eigenvals = []
    for idx in maj_spin_line_idx:
        maj_spin_eigenvals += [get_orbital_energies_for_idx(lines, idx)]
    min_spin_eigenvals = []
    for idx in min_spin_line_idx:
        min_spin_eigenvals += [get_orbital_energies_for_idx(lines, idx)]

    to_return = {
            'majority_spin_orbital_energies': maj_spin_eigenvals,
            'minority_spin_orbital_energies': min_spin_eigenvals,
            'HOMO': get_homo_level(lines),
            'LUMO': get_lumo_level(lines),
            'orbital_energy_unit': 'electron-volt',
            }

    return to_return


def extract_lambda_matrix_elements(lines, units='eV'):

    assert units in ['Ha', 'eV']

    Ha_to_eV = 27.211396132

    matrix = []
    for line in lines:
        try:
            if units == 'eV':
                matrix_elements = [
                        Ha_to_eV * float(element) for element in line.split()
                        ]
                assert len(matrix_elements) > 0
                matrix += [matrix_elements]
            else:
                matrix_elements = [float(element) for element in line.split()]
                assert len(matrix_elements) > 0
                matrix += [matrix_elements]
        except:
            break

    return matrix


def get_lambda_matrices(lines):

    assert len(lines) > 0

    lambda_patterns = DEFAULT_LAMBDA_PATTERNS
    next_line_pattern = "print only first"
    lambda_matrices_lines = find_matching_lines_idx_for_patterns(
            lines,
            lambda_patterns
            )

    to_return = {}
    for idx in lambda_matrices_lines:
        if lines[idx].split()[-1] == '1':
            spin_component = "up"
        else:
            spin_component = "dw"
        assert next_line_pattern in lines[idx+1]
        matrix_name = "lambda_matrix_%s" % spin_component
        to_return[matrix_name] = extract_lambda_matrix_elements(lines[idx+2:])

    return to_return


def extract_centers_matrix_elements(lines):

    separator = '---'

    lines_containing_centers = []
    for line in lines:
        if separator in line:
            lines_containing_centers += [line]
        else:
            break

    centers = []
    for line in lines_containing_centers:
        coordinates = [float(item) for item in line.split(separator)[1].split()]
        centers += [coordinates]

    return centers


def get_centers(lines):

    assert len(lines) > 0

    center_patterns = DEFAULT_CENTER_PATTERNS
    center_lines = find_matching_lines_idx_for_patterns(
            lines,
            center_patterns
            )

    to_return = {}
    for idx in center_lines:
        if lines[idx].split()[-1] == '1':
            spin_component = "up"
        else:
            spin_component = "dw"
        matrix_name = "centers_matrix_%s" % spin_component
        to_return[matrix_name] = extract_centers_matrix_elements(lines[idx+2:])

    return to_return


# --------------------
# -- Parsing utilities
# --------------------
def parse_cp_output(filename, e_patterns_dict=DEFAULT_ENERGY_PATTERNS):

    assert path.exists(filename)
    assert path.isfile(filename)

    gathered_data = {}
    with open(filename, 'r') as fin:
        raw_data = fin.readlines()
        gathered_data['energy_unit'] = 'Hartree'
        for line in raw_data:
            # -- extraction of the energies
            for name, pattern in e_patterns_dict.items():
                if pattern in line:
                    gathered_data[name] = get_energy_from_line(line)
        # -- extraction of the orbital energies
        orb_energies_dict = get_orbital_energies(raw_data)
        gathered_data["orbital_energies"] = orb_energies_dict
        # -- extracting lambda matrices
        lambda_matrices = get_lambda_matrices(raw_data)
        gathered_data["lambda_matrices"] = lambda_matrices
        # -- extracting minimizing orbital centers
        centers = get_centers(raw_data)
        gathered_data["minimizing_orbital_centers"] = centers
        return gathered_data


def get_alpha_value(
        lumo_reference,
        homo_reference,
        lumo_current,
        homo_current,
        previous_alpha=None):
    """This script will extract HOMO and LUMO levels from a set of CP output
    files and compute the "alpha" parameter used in Non-Koopmans type of
    calculations.
    If a previous value for alpha is given, on top of the other required inputs
    to this function, then the second-method is used to compute the updated
    value for alpha.

    Explanations
    ============

    ``lumo_reference`` and ``homo_reference`` should be the output files for
    respectively the cation and neutral systems. These are the "references"
    because they correspond to the calculation with alpha = 1.0
    ``lumo_current`` and ``homo_current`` correspond to the most recent value of
    alpha. In particular, if ``previous_alpha`` is None, then these should
    correspond respectively to the cation and neutral systems computed in LDA.
    If ``previous_alpha`` is *not* None, these should correspond to respectively
    the cation and neutral systems of the latest a-NK calculations.

    References
    ==========

    The calculation of alpha from scratch (default behavior of this script) is
    possible through the use of formula (43) in the `reference paper`_
    The calculation of an updated value for alpha is made possible through the
    use of formula (44) from the same paper.

    .. _reference paper: http://prb.aps.org/abstract/PRB/v82/i11/e115121
    """
    lumo_nk_reference = \
    parse_cp_output(lumo_reference)["orbital_energies"]["HOMO"]
    homo_nk_reference = \
    parse_cp_output(homo_reference)["orbital_energies"]["HOMO"]
    lumo_current = parse_cp_output(lumo_current)["orbital_energies"]["HOMO"]
    homo_current = parse_cp_output(homo_current)["orbital_energies"]["HOMO"]

    if previous_alpha is None:
        alpha = (lumo_current - homo_current) / ((lumo_current - homo_current) -
                (lumo_nk_reference - homo_nk_reference))
    else:
        assert 0. <= previous_alpha <= 1.
        alpha = previous_alpha + (1. - previous_alpha) * ((lumo_current -
            homo_current) / ((lumo_current - homo_current) - (lumo_nk_reference
                - homo_nk_reference)))
    assert 0. <= alpha <= 1.
    return alpha
