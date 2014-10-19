#!/usr/bin/env python


from elements import ELEMENTS
from os import path
from string import join
import hashlib


# -- getting the dictionnaries that map element symbols to atomic numbers
Z_l = [element.number for element in ELEMENTS]
Symbol_l = [element.symbol for element in ELEMENTS]

symbol_to_atomic_number = {}
atomic_number_to_symbol = {}
for symbol, Z in zip(Symbol_l, Z_l):
    symbol_to_atomic_number[symbol] = int(Z)
    atomic_number_to_symbol[str(Z)] = symbol


def join_strings(list_of_lists):

    final_string = ''
    for a_list in list_of_lists:
        new_list = [str(item) for item in a_list]
        a_list_string = join(new_list)
        final_string += a_list_string

    return final_string


# -- Utilities
def get_valence_electron_number(symbols):

    n_elec = 0
    for S in symbols:
        if S == 'H' or S == 'Li' or S == 'Na' or S == 'K':
            n_elec += 1
        elif S == 'B' or S == 'Al' or S == 'Ga':
            n_elec += 3
        elif S == 'C' or S == 'Si' or S == 'Ge':
            n_elec += 4
        elif S == 'N' or S == 'P' or S == 'As':
            n_elec += 5
        elif S == 'O' or S == 'S' or S == 'Se':
            n_elec += 6
        elif S == 'F' or S == 'Cl' or S == 'Br' or S == 'I':
            n_elec += 7
        else:
            raise ValueError('unrecognized element: %s' % S)

    return n_elec


def normalize_xyz(xyz_line_elements):

    xyz_normalized = []
    for line_elements in xyz_line_elements:
        symbol, x, y, z = line_elements
        try:
            symbol_int = int(symbol)
            xyz_normalized += [
                    [
                        atomic_number_to_symbol[str(symbol_int)],
                        float(x), float(y), float(z)
                    ]
                ]
        except:
            xyz_normalized += [
                    [
                        str(symbol),
                        float(x), float(y), float(z)
                    ]
                ]

    return xyz_normalized


def robust_xyz_extract(fname):

    assert path.isfile(fname)

    with open(fname, 'r') as fin:
        lines = fin.readlines()

    xyz_lines = []
    for line in lines:
        line_items = line.split()
        if len(line_items) == 4:
            try:
                s, x, y, z = line_items
                x, y, z = float(x), float(y), float(z)
                xyz_lines += [[s, x, y, z]]
            except:
                continue

    return normalize_xyz(xyz_lines)


def generate_single_point_energy(
        xyz_lines,
        output_fname,
        scratch_dir='/path/to/gaussian09/scratch/dir/',
        basis_sets=['aug-cc-pVDZ', 'aug-cc-pVTZ'],
        chemistry='PBE1PBE',
        memory_footprint='2GB',
        charge=0,
        ):

    normalized_xyz = normalize_xyz(xyz_lines)
    symbol_l = [item[0] for item in normalized_xyz]

    n_electrons = get_valence_electron_number(symbol_l)
    n_electrons -= charge
    multiplicity = n_electrons%2 + 1

    xyz_file_string = join_strings(normalized_xyz)
    scratch_basename = hashlib.md5(xyz_file_string).hexdigest()
    with open(output_fname, 'w') as fout:
        scratch_path = path.join(scratch_dir, scratch_basename)
        fout.write('%Chk=' + scratch_path + '\n')
        fout.write('%Mem=' + memory_footprint + '\n')
        fout.write('# ' + chemistry + ' / ' + basis_sets[0] + '\n')
        fout.write('\n')
        fout.write('Single Point Energy Calculation\n')
        fout.write('\n')
        fout.write('%i %i\n' % (charge, multiplicity))
        for xyz_line in normalized_xyz:
            s, x, y, z = xyz_line
            fout.write('%s  %10.6f  %10.6f  %10.6f\n' % (
                    s, x, y, z
                ))
        fout.write('\n')
        if len(basis_sets) > 1:
            for basis in basis_sets[1:]:
                fout.write('\n')
                fout.write('--Link1--\n')
                fout.write('%Chk=' + scratch_path + '\n')
                fout.write('%Mem=' + memory_footprint + '\n')
                fout.write('# ' + chemistry + ' / ' + basis + \
                           ' Geom=check Guess=Read\n')
                fout.write('\n')
                fout.write('new basis %s\n' % basis)
                fout.write('\n')
                fout.write('%i %i\n' % (charge, multiplicity))
                fout.write('\n')

    return


def generate_geometry_optimization(
        xyz_lines,
        output_fname,
        scratch_dir='/path/to/gaussian09/scratch/dir',
        relaxation_basis_set='6-31+G*',
        single_point_basis_set='aug-cc-pVTZ',
        chemistry='PBE1PBE',
        memory_footprint='2GB',
        charge=0,
        ):

    normalized_xyz = normalize_xyz(xyz_lines)
    symbol_l = [item[0] for item in normalized_xyz]

    n_electrons = get_valence_electron_number(symbol_l)
    n_electrons -= charge
    multiplicity = n_electrons%2 + 1

    xyz_file_string = join_strings(normalized_xyz)
    scratch_basename = hashlib.md5(xyz_file_string).hexdigest()
    with open(output_fname, 'w') as fout:
        scratch_path = path.join(scratch_dir, scratch_basename)
        fout.write('%Chk=' + scratch_path + '\n')
        fout.write('%Mem=' + memory_footprint + '\n')
        fout.write('# PM6 / Opt=(Tight,CalcAll,maxcycles=100)\n')
        fout.write('\n')
        fout.write('Initial Geometry Optimization\n')
        fout.write('\n')
        fout.write('%i %i\n' % (charge, multiplicity))
        for xyz_line in normalized_xyz:
            s, x, y, z = xyz_line
            fout.write('%s  %10.6f  %10.6f  %10.6f\n' % (
                    s, x, y, z
                ))
        fout.write('\n')
        fout.write('--Link1--\n')
        fout.write('%Chk=' + scratch_path + '\n')
        fout.write('%Mem=' + memory_footprint + '\n')
        fout.write('# ' + chemistry + ' / ' + relaxation_basis_set + \
                   ' Opt=(Tight,CalcFC,maxcycles=100) Geom=AllCheck\n')
        fout.write('\n')
        fout.write('Actual Geometry Optimization\n')
        fout.write('\n')
        fout.write('%i %i\n' % (charge, multiplicity))
        fout.write('\n')
        fout.write('\n')
        fout.write('--Link1--\n')
        fout.write('%Chk=' + scratch_path + '\n')
        fout.write('%Mem=' + memory_footprint + '\n')
        fout.write('# ' + chemistry + ' / ' + single_point_basis_set + \
                   ' Geom=AllCheck Guess=Read\n')
        fout.write('\n')
        fout.write('Single Point Energy\n')
        fout.write('\n')
        fout.write('%i %i\n' % (charge, multiplicity))
        fout.write('\n')

    return
