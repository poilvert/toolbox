#!/usr/bin/env python


from os import path
import numpy as np
from xml2dict_library import eigxml_to_dict
from xml2dict_library import evcxml_to_dict
from xml2dict_library import gkvectorxml_to_dict
from xml2dict_library import xml_to_dict


def get_eigenvalues_and_eigenvectors(outdir_path):
    """
    This high-level routine will extract all the eigenvalues and eigenfunctions
    for all the k points in a Quantum Espresso outdir/restart directory.

    Warning
    =======
    This function assumes that all the binary files corresponding to the
    eigenvalues, gk-vectors and eigenfunctions have been converted to XML
    format. If it is not the case, please use the ``binary_to_xml.py`` utility
    first.
    """

    # -- first we build the path to the master data file
    datafile_path = path.join(outdir_path, 'data-file.xml')

    # -- reading all the information
    data_dict = xml_to_dict(datafile_path)
    key_l = data_dict.keys()

    # -- some checks
    keys_to_check = ['EIGENVALUES', 'EIGENVECTORS']
    for key in keys_to_check:
        assert key in key_l

    eigval_keys = data_dict['EIGENVALUES'].keys()
    eigvec_keys = data_dict['EIGENVECTORS'].keys()

    eigval_kpoint_keys = set(
            [key for key in eigval_keys if key.startswith('K-POINT')]
            )
    eigvec_kpoint_keys = set(
            [key for key in eigvec_keys if key.startswith('K-POINT')]
            )
    assert (eigval_kpoint_keys == eigvec_kpoint_keys)
    assert len(eigval_kpoint_keys) >= 1
    kpoint_keys = eigval_kpoint_keys

    # -- it is now time to extract eigenvalues and eigenvectors at each k point
    kpoint_dict = {}
    kpoint_dict['max_number_of_gk-vectors'] = int(
                data_dict['EIGENVECTORS']['MAX_NUMBER_OF_GK-VECTORS']['value']
            )
    for key in kpoint_keys:

        keigval_info = data_dict['EIGENVALUES'][key]
        keigvec_info = data_dict['EIGENVECTORS'][key]
        ik = int(key.split('.')[-1])

        # -- k vector coordinates
        k_vec = np.array(keigval_info['K-POINT_COORDS']['value'])
        # -- weight of the k vector for Brillouin Zone summations
        k_weight = float(keigval_info['WEIGHT']['value'])

        # -- we need the Gk vectors
        gkvectors = gkvectorxml_to_dict(path.join(
                outdir_path,
                keigvec_info['GK-VECTORS']['iotk_link'].replace('dat', 'xml')
                ))
        gk_num = int(keigvec_info['NUMBER_OF_GK-VECTORS']['value'])

        # -- finding the keys to the eigval and eigvec data files
        datafile_keys = [mykey for mykey in keigval_info.keys()
                               if mykey.startswith('DATAFILE')]
        wfcfile_keys = [mykey for mykey in keigvec_info.keys()
                              if mykey.startswith('WFC')]
        assert len(datafile_keys) == 1 or len(datafile_keys) == 2

        # -- paths to the XML files of the eigenvalues and eigenvectors
        eigval_xml_path_l = [path.join(
                outdir_path,
                keigval_info[datafile_key]['iotk_link'].replace('dat', 'xml')
                ) for datafile_key in datafile_keys
                ]
        eigvec_xml_path_l = [path.join(
                outdir_path,
                keigvec_info[wfcfile_key]['iotk_link'].replace('dat', 'xml')
                ) for wfcfile_key in wfcfile_keys
                ]

        # -- eigenvalue dictionnary
        keigval_dict = {}
        for eigval_xml_path in eigval_xml_path_l:
            eigval_key = path.basename(eigval_xml_path)
            keigval_dict[eigval_key] = eigxml_to_dict(eigval_xml_path)
            assert keigval_dict[eigval_key]['ik'] == ik

        # -- eigenvector dictionnary
        keigvec_dict = {}
        for eigvec_xml_path in eigvec_xml_path_l:
            eigvec_key = path.basename(eigvec_xml_path)
            keigvec_dict[eigvec_key] = evcxml_to_dict(eigvec_xml_path)
            #assert keigvec_dict[eigvec_key]['ik'] == ik

        # -- we now add the k point data to the global dictionnary
        kpoint_dict[key] = {
                    'coordinates': k_vec,
                    'weight': k_weight,
                    'eigenvalues': keigval_dict,
                    'eigenvectors': keigvec_dict,
                    'gkvectors': gkvectors,
                    'number_of_gk-vectors': gk_num,
                }

    return kpoint_dict
