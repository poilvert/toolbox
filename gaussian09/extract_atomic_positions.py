#!/usr/bin/env python

import sys
import optparse
from os import path
import numpy as np


DEFAULT_ATO_POS_BLK_STR = 'Standard orientation:'
DEFAULT_LINE_SEPARATOR_STR = '------------------'


def extract_Z_and_R(
        fname,
        ato_pos_blk_str=DEFAULT_ATO_POS_BLK_STR,
        line_separator_str=DEFAULT_LINE_SEPARATOR_STR,
        reorder_axes=False,
        ):
    """
    Finds atomic numbers and positions in the Standard Orientation of the last
    molecule from a Gaussian09 output file and dump these metadata into a file

    Re-order axes means that the z axis will be the one for which the standard
    deviation about the mean coordinate is the smallest.
    In the particular case that the molecule is 2D, then the z axis will be
    orthogonal to the molecular plane.
    """

    assert path.exists(fname)
    assert path.isfile(fname)

    with open(fname, 'r') as fin:
        lines = fin.readlines()

    n_lines = len(lines)
    assert n_lines > 1

    start_line_idx = 0
    for idx, line in enumerate(reversed(lines)):
        if ato_pos_blk_str in line:
            start_line_idx = n_lines-1-idx
            break

    separator_line_indices = []
    while ( (len(separator_line_indices) <= 2) and \
            (start_line_idx < n_lines) ):
        line = lines[start_line_idx]
        if line_separator_str in line:
            separator_line_indices += [start_line_idx]
        start_line_idx += 1

    assert len(separator_line_indices) == 3

    atomic_block_start_line = separator_line_indices[1] + 1
    atomic_block_stop_line = separator_line_indices[-1]

    atomic_data = []
    for idx in range(atomic_block_start_line, atomic_block_stop_line):
        line_items = [item for item in lines[idx].split()]
        assert len(line_items) == 6
        Z = int(line_items[1])
        X1, X2, X3 = float(line_items[3]), float(line_items[4]), \
                     float(line_items[5])
        atomic_data += [[Z, X1, X2, X3]]

    assert len(atomic_data) > 0

    atomic_data = np.array(atomic_data)

    if reorder_axes:
        positions = atomic_data[:, 1:].copy()
        assert positions.ndim == 2
        assert positions.shape[1] == 3
        std_devs = positions.std(0)
        assert std_devs.shape == (3,)
        min_idx = std_devs.argmin()
        if min_idx == 0:
            positions = np.roll(positions, 2, axis=1)
        elif min_idx == 1:
            positions = np.roll(positions, 1, axis=1)
        atomic_data[:, 1:] = positions

    return atomic_data


def dump_positions(atomic_data, out_fname):

    with open(out_fname, 'w') as fout:
        fout.write('%i\n' % (len(atomic_data),))
        fout.write("# Zi, xi, yi, zi / Zi = atomic number /" + \
                   " (xi, yi, zi) cartesian coordinates in Angstroms\n")
        for atom_data in atomic_data:
            Z, X1, X2, X3 = atom_data
            fout.write('%2i   %12.7f   %12.7f   %12.7f\n' % (
                int(np.around(Z)), X1, X2, X3))

    return


# -------------------
# -- CLI for the code
# -------------------

def main():

    parser = optparse.OptionParser()
    parser.add_option('-i','--input',
                    dest="input_fname",
                    type="str",
                    help="path to the Gaussian09 output file")
    parser.add_option('-o','--output',
                    dest="output_fname",
                    type="str",
                    default="if_default",
                    help="name of the output position file")
    parser.add_option('-r','--reorder',
                    action="store_true",
                    dest="reorder",
                    default=False,
                    help="set z perpendicular to the molecular plane")
    options, remainder = parser.parse_args()
    if len(sys.argv[1:])==0:
        parser.print_help()
        sys.exit(0)

    if options.output_fname == "if_default":
        output_fname = options.input_fname + '.pos'
    else:
        output_fname = options.output_fname

    if options.reorder:
        atomic_data = extract_Z_and_R(
                options.input_fname,
                output_fname,
                reorder_axes=True,
                )
    else:
        atomic_data = extract_Z_and_R(
                options.input_fname,
                output_fname,
                )

    dump_positions(atomic_data, output_fname)

    return


if __name__ == "__main__":
    main()
