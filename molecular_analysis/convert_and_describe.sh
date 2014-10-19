#!/bin/bash

# -- Finding the executable for OpenBabel
BAB=`which babel`

# -- Finding the executable for "checkmol/matchmol" (you can find this code
# here: http://merian.pch.univie.ac.at/~nhaider/cheminf/cmmm.html)
CHK="./checkmol"

INPUT=$*

# -- first we convert the input file in XYZ format into a .mol format
$BAB -ixyz $INPUT -omol $INPUT.mol

# -- then we get the functional group information
$CHK $INPUT.mol > $INPUT.mol.info
