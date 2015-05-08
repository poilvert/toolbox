#!/usr/bin/env python


from xml.etree import ElementTree
import collections
import numpy as np
from copy import deepcopy
import os


__all__ = [
    'densityxml_to_dict',
    'evcxml_to_dict',
    'gkvectorxml_to_dict',
    'xml_to_dict',
    'binary_file_paths'
]


# -- Utility functions
def leaf_paths_generator(tree, cur=()):

    if not isinstance(tree, dict):
        yield cur
    else:
        for n, s in tree.items():
            for path in leaf_paths_generator(s, cur+(n,)):
                yield path

def leaf_paths_in_dict(tree):
    return list(leaf_paths_generator(tree))

def get_dict_leaf_val_from_path(my_dict, leaf_path):
    return reduce(dict.get, leaf_path, my_dict)

def update_item_dict(item_dict):

    updated_item_dict = deepcopy(item_dict)

    text = updated_item_dict['text']
    item_type = updated_item_dict['type']
    item_size = int(updated_item_dict['size'])
    item_keys = updated_item_dict.keys()

    if item_type == 'integer':
        text_items = [
                item.strip()
                for item in text.strip().split('\n')
                if item.strip() != ''
                     ]

        if 'columns' in item_keys:
            updated_item_dict['columns'] = int(updated_item_dict['columns'])
            text_items = [
                    [float(number) for number in line.split()]
                    for line in text_items
                    ]
            assert sum([len(sub_l) for sub_l in text_items]) == item_size
        else:
            text_items = [
                    int(number) for number in text_items
                    ]
            assert len(text_items) == item_size

        if 'UP' in item_keys:
            updated_item_dict['UP'] = int(updated_item_dict['UP'])
        if 'DW' in item_keys:
            updated_item_dict['DW'] = int(updated_item_dict['DW'])

    elif item_type == 'logical':
        text_items = [
                item.strip()
                for item in text.strip().split('\n')
                if item.strip() != ''
                ]

        assert len(text_items) == item_size

        text_items = [
                True if item.upper() == 'T' else False
                for item in text_items
                ]

    elif item_type == 'character':
        text_items = [
                item.strip()
                for item in text.strip().split('\n')
                if item.strip() != ''
                ]

        assert len(text_items) == item_size

    elif item_type == 'real':
        text_items = [
                item.strip()
                for item in text.strip().split('\n')
                if item.strip() != ''
                ]

        if 'UP' in item_keys:
            updated_item_dict['UP'] = float(updated_item_dict['UP'])
        if 'DW' in item_keys:
            updated_item_dict['DW'] = float(updated_item_dict['DW'])
        if 'columns' in item_keys:
            updated_item_dict['columns'] = int(updated_item_dict['columns'])
            text_items = [
                    [float(number) for number in line.split()]
                    for line in text_items
                    ]
            assert sum([len(sub_l) for sub_l in text_items]) == item_size
        else:
            text_items = [
                    float(number) for number in text_items
                    ]
            assert len(text_items) == item_size

    elif item_type == 'complex':
        text_items = [
                map(float, item.strip().split(','))
                for item in text.strip().split('\n')
                if item.strip() != ''
                ]
        assert len(text_items) == item_size

    else:
        raise ValueError('unknown type: "%s"' % item_type)

    updated_item_dict['size'] = item_size
    if len(text_items) == 1:
        updated_item_dict['value'] = text_items[0]
    else:
        updated_item_dict['value'] = text_items
    updated_item_dict.pop('text')

    return updated_item_dict

def item_list_to_dict(item_l):

    keys = [key for key, val in item_l]

    if 'type' not in keys:
        return dict(item_l)
    else:
        item_dict = dict(item_l)
        item_dict = update_item_dict(
                    item_dict
                    )
        return item_dict

def flatten_to_tuples_list(l):

    for el in l:
        if (
            isinstance(el, collections.Iterable) and not
            isinstance(el, tuple)
            ):
            for sub in flatten_to_tuples_list(el):
                yield sub
        else:
            yield el

def get_children_and_tags(element):

    children = element.getchildren()
    tags = [child.tag for child in children]
    return children, tags

def recursive_element_extraction(element, full_item_list):

    children, tags = get_children_and_tags(element)

    if len(children) == 0:
        item_l = element.items()
        if element.text is not None:
            item_l += [('text', element.text.strip())]
        else:
            item_l += [('text', None)]
        full_item_list += [item_l]
        item_dict = item_list_to_dict(item_l)
        return item_dict
    else:
        dico = {}
        for child, tag in zip(children, tags):
            dico[tag] = recursive_element_extraction(
                            child,
                            full_item_list
                            )
        return dico


# -- Main extraction routine
def xml_to_dict(xml_path, with_item_l=False):
    """
    Reads an XML formatted file produced by IOTK into a Python dictionnary
    """

    assert os.path.isfile(xml_path)

    tree = ElementTree.parse(xml_path)
    root = tree.getroot()

    full_item_list = []
    data_dict = recursive_element_extraction(root, full_item_list)

    if with_item_l:
        return data_dict, full_item_list
    else:
        return data_dict


def densityxml_to_dict(xml_path):
    """
    Extracts the total charge density or spin polarization

    Conventions
    ===========
    In order to extract the density as a 3D tensor I followed the conventions
    for storing the charge density into the ``charge_density.xml`` or
    ``spin_polarization.xml`` files. The details are given in Quantum Espresso's
    Developer Manual found here::

        http://www.quantum-espresso.org/wp-content/uploads/Doc/developer_man.pdf
    """

    raw_dict = xml_to_dict(xml_path)
    key_l = raw_dict.keys()

    assert 'CHARGE-DENSITY' in key_l

    raw_dict = raw_dict['CHARGE-DENSITY']
    key_l = raw_dict.keys()

    assert 'INFO' in key_l

    info_dict = {}
    for key, val in raw_dict['INFO'].items():
        if key == 'nr1' or key == 'nr2' or key == 'nr3':
            info_dict[key] = int(val)

    # -- now the Fourier coefficients of the charge density
    nr1 = info_dict['nr1']
    nr2 = info_dict['nr2']
    nr3 = info_dict['nr3']
    assert len([key for key in key_l if key != 'INFO']) == nr3
    density_arr = np.empty((nr1, nr2, nr3), dtype=np.float64)
    for k in range(nr3):
        key = 'z.%i' % (k+1,)
        assert key in key_l
        size = int(raw_dict[key]['size'])
        value = np.array(raw_dict[key]['value'], dtype=np.float64)
        assert value.size == size
        assert size == nr1 * nr2
        value = value.reshape(nr2, nr1)
        density_arr[..., k] = value.T
    info_dict['density'] = density_arr

    return info_dict


def gvectorxml_to_dict(xml_path):
    """
    Extracts the global G vector grid information
    """

    raw_dict = xml_to_dict(xml_path)
    key_l = raw_dict.keys()

    assert 'INFO' in key_l
    assert 'g' in key_l

    info_dict = {}
    for key, val in raw_dict['INFO'].items():
        if key == 'gamma_only':
            info_dict['gamma_only'] = False if val == "F" else True
        elif key == 'do_wf_cmplx':
            info_dict['do_wf_cmplx'] = False if val == "F" else True
        elif key == 'gvect_number':
            info_dict['gvect_number'] = int(val)
        elif key == 'nr1s':
            info_dict['nr1s'] = int(val)
        elif key == 'nr2s':
            info_dict['nr2s'] = int(val)
        elif key == 'nr3s':
            info_dict['nr3s'] = int(val)
        else:
            info_dict[key] = val

    size = int(raw_dict['g']['size'])
    cols = int(raw_dict['g']['columns'])
    arr = np.array(raw_dict['g']['value'], dtype=np.int)
    assert arr.size == size
    assert arr.shape == (size / cols, cols)
    info_dict['grid'] = arr

    return info_dict


def gkvectorxml_to_dict(xml_path):
    """
    Extracts the G vector grid information for the given k point
    """

    raw_dict = xml_to_dict(xml_path)
    key_l = raw_dict.keys()

    info_dict = {}
    for key in key_l:
        if key == 'GAMMA_ONLY':
            info_dict['gamma_only'] = raw_dict[key]['value']
        elif key == 'DO_WF_CMPLX':
            info_dict['do_wf_cmplx'] = raw_dict[key]['value']
        elif key == 'MAX_NUMBER_OF_GK-VECTORS':
            info_dict['max_number_of_gk_vectors'] = int(raw_dict[key]['value'])
        elif key == 'NUMBER_OF_GK-VECTORS':
            info_dict['number_of_gk_vectors'] = int(raw_dict[key]['value'])
        elif key == 'INDEX':
            size = int(raw_dict[key]['size'])
            arr = np.array(raw_dict[key]['value'], dtype=np.float64)
            assert arr.ndim == 1
            assert arr.size == size
            info_dict['index'] = np.rint(arr).astype(np.int)
        elif key == 'GRID':
            size = int(raw_dict[key]['size'])
            cols = int(raw_dict[key]['columns'])
            arr = np.array(raw_dict[key]['value'], dtype=np.int)
            assert arr.size == size
            assert arr.shape == (size / cols, cols)
            info_dict['grid'] = arr
        elif key == 'K-POINT_COORDS':
            info_dict['k-point_coords'] = raw_dict[key]

    return info_dict


def eigxml_to_dict(xml_path):
    """
    Extracts the eigenvalue information for a given k point along with the
    occupations for each eigenvalue
    """

    raw_dict = xml_to_dict(xml_path)

    info_dict = {}
    info_dict['eigenvalues'] = np.array(
            raw_dict['EIGENVALUES']['value'],
            dtype=np.float64
            )
    info_dict['ik'] = int(raw_dict['INFO']['ik'])
    info_dict['nbnd'] = int(raw_dict['INFO']['nbnd'])
    info_dict['occupations'] = np.array(
            raw_dict['OCCUPATIONS']['value'],
            dtype=np.float64
            )
    info_dict['units_for_energies'] = raw_dict['UNITS_FOR_ENERGIES']['UNITS']

    return info_dict


def evcxml_to_dict(xml_path):
    """
    Extracts the eigenvectors information for a given k point
    """

    raw_dict = xml_to_dict(xml_path)
    key_l = raw_dict.keys()

    assert 'INFO' in key_l

    # -- global informations on the eigenvectors
    info_dict = {}
    info_key_l = raw_dict['INFO'].keys()
    for key in info_key_l:
        val = raw_dict['INFO'][key]
        if key == 'ngw':
            info_dict['ngw'] = int(val)
        elif key == 'igwx':
            info_dict['igwx'] = int(val)
        elif key == 'do_wf_cmplx':
            info_dict['do_wf_cmplx'] = False if val == "F" else True
        elif key == 'gamma_only':
            info_dict['gamma_only'] = False if val == "F" else True
        elif key == 'nbnd':
            info_dict['nbnd'] = int(val)
        elif key == 'ik':
            info_dict['ik'] = int(val)
        elif key == 'nk':
            info_dict['nk'] = int(val)
        elif key == 'ispin':
            info_dict['ispin'] = int(val)
        elif key == 'nspin':
            info_dict['nspin'] = int(val)
        elif key == 'scale_factor':
            info_dict['scale_factor'] = float(val)
        else:
            info_dict[key] = val

    # -- extraction of the eigenvector Fourier components
    for key in key_l:
        if key != 'INFO':
            assert 'evc' in key
            evc_dict = raw_dict[key]
            evc_dict['size'] = int(evc_dict['size'])
            arr = np.array(evc_dict['value'], dtype=np.float64)
            assert arr.ndim == 2
            size, dim = arr.shape
            assert dim == 2
            assert size == evc_dict['size']
            evc_arr = np.empty((size,), dtype=np.complex128)
            evc_arr.real = arr[:, 0]
            evc_arr.imag = arr[:, 1]
            evc_dict['value'] = evc_arr
            info_dict[key] = evc_dict

    return info_dict


def binary_file_paths(
        outdir_path,
        data_file_name='data-file.xml'
        ):
    """
    This routine returns all the binary file paths found inside file
    ``data-file.xml`` of the Quantum Espresso outdir/restart directory.
    """

    assert os.path.isdir(outdir_path)

    data_file_path = os.path.join(outdir_path, data_file_name)
    assert os.path.isfile(data_file_path)

    data_dict = xml_to_dict(data_file_path)

    # -- all ``iotk_link`` paths
    all_leaf_paths = leaf_paths_in_dict(data_dict)
    iotk_leaf_paths = [
            leaf_path for leaf_path in all_leaf_paths
            if leaf_path[-1] == 'iotk_link'
            ]
    iotk_paths = [
            get_dict_leaf_val_from_path(data_dict, leaf_path)
            for leaf_path in iotk_leaf_paths
            ]

    return iotk_paths
