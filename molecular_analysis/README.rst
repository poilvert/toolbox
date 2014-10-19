Extracting molecular structures from SMILES strings
===================================================

The point here is to leverage the power of OpenBabel_ that can read molecular descriptors like SMILES and convert these descriptors into 3D molecular structures and dump the resulting structure in any format like XYZ for example.

The script `extract_xyz.sh` does just this.

The script `convert_and_describe.sh` reads an XYZ formatted file to introspect a molecular structure and determine all the different distinct organic groups present. These groups can be alkyles, nitriles, etc... for over 200 distincts groups. This has been made possible with the use of the excellent tool designed by **Norbert Haider**. See checkmol_


.. _OpenBabel: http://openbabel.org/wiki/Main_Page
.. _checkmol: http://merian.pch.univie.ac.at/~nhaider/cheminf/cmmm.html
