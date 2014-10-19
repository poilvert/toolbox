#!/bin/bash

# -- finding the executable of OpenBabel
BABEL=`which obabel`

# -- Input filename containing the molecules in SMILES format
SMIINPUT="<change_me>.smi"

# -- Ouput filename that will contain the same molecules but in XYZ format
XYZOUTPUT="<change_me>.xyz"

# -- to extract the XYZ positions from the raw SMILES strings
$BABEL -ismi $SMIINPUT -oxyz -O$XYZOUTPUT --gen3d -c
