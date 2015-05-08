#!/usr/bin/env python


import os
import sys
import optparse
from os import path


# -- Internal Utilities
def _listfiles(root_folder_path, filename_extension='.dat'):
    for root, folders, files in os.walk(root_folder_path):
        for filename in folders + files:
            if filename.endswith(filename_extension):
                yield path.join(root, filename)

def _call_capture_output(cmdline, cwd=None, error_on_nonzero=True):
    """
    Returns
    =======

    a tuple (return code, stdout_data, stderr_data)

    Reference
    =========

    from Andreas Klochner
    https://github.com/inducer/pytools/blob/master/pytools/prefork.py#L39

    modified by Nicolas Poilvert, 2011
    """

    from subprocess import Popen, PIPE
    import shlex

    try:
        args = shlex.split(cmdline)
        popen = Popen(args, cwd=cwd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout_data, stderr_data = popen.communicate()
        if error_on_nonzero and popen.returncode:
            raise ExecError("status %d invoking '%s': %s"
                    % (popen.returncode, "".join(cmdline), stderr_data))
        return popen.returncode, stdout_data, stderr_data
    except OSError, e:
        raise ExecError("error invoking '%s': %s"
                % ("".join(cmdline), e))

class ExecError(OSError):
    pass


# -- Main conversion routines
def convert_outdir(outdir_path, iotk_exe_path):
    """
    Front-end routine that will convert every binary data file inside the outdir
    to XML files
    """

    assert path.isdir(outdir_path)

    dat_paths = [
            path.abspath(file_path) for file_path in _listfiles(outdir_path)
            ]

    for dat_path in dat_paths:
        dat2xml(iotk_exe_path, dat_path)

    return


def dat2xml(iotk_exe_path, dat_path):
    """
    This function uses the IOTK converter in Quantum Espresso to convert a
    binary file (from IOTK, usually ending in ``.dat``) into a universally
    readable XML file.

    Parameters
    ==========

    ``iotk_exe_path``: str
        path to the ``iotk`` conversion utility

    ``dat_path``: str
        path to the binary file to convert into XML

    Returns
    =======

    nothing. The code will just produce an XML file alongside the original
    binary file. If an error occurs the code will relay the error message
    """

    assert path.isfile(iotk_exe_path)
    assert path.isfile(dat_path)

    # -- transform data file path into absolute path
    abs_dat_path = path.abspath(dat_path)

    # -- building the data XML absolute path
    dirname = path.dirname(abs_dat_path)
    basename = path.basename(abs_dat_path)
    assert basename.endswith('.dat')
    abs_xml_path = path.join(dirname, basename.replace('.dat', '.xml'))

    # -- absolute path to the IOTK executable
    abs_iotk_exe_path = path.abspath(iotk_exe_path)
    assert (
            path.isfile(abs_iotk_exe_path) and
            os.access(abs_iotk_exe_path, os.X_OK)
           )

    cmdline = ("%s %s %s %s" % (
        abs_iotk_exe_path, 'convert',
        abs_dat_path, abs_xml_path)
        )

    return_code, stdout, stderr = _call_capture_output(cmdline)

    assert 'iotk:' in stderr

    if return_code == 0:
        print 'converted "%s" into "%s"' % (
                path.basename(abs_dat_path),
                path.basename(abs_xml_path)
                )
    else:
        print 'An error occured while executing command :'
        print '%s' % cmdline
        raise ExecError("return code %s" % return_code)

    return


# -- CLI
def main():
    """
    Command Line Interface (CLI) to the conversion routine
    """

    parser = optparse.OptionParser()
    parser.add_option('-o', '--outdir',
                dest="outdir_path",
                help="path to Quantum Espresso outdir/restart directory")
    parser.add_option('-e', '--exe',
                dest="iotk_exe_path",
                help="path to the IOTK executable (iotk not iotk.x)")

    options, remainder = parser.parse_args()

    # -- if the program is called without options it prints the help menu
    if len(sys.argv[1:])==0:
        parser.print_help()
        sys.exit(0)

    convert_outdir(options.outdir_path,
                   options.iotk_exe_path)

    return


if __name__ == '__main__':
    main()
