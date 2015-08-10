#!/usr/bin/env bash

## This script should be copied into the directory where all of the
## input files are going to be generated!

# What kind of calculation are we doing? One of {`freq`, `eda`,
# `freq_noCT`}.
CALCTYPE=eda

# The *root* directory for the MBE code.
MBE_ROOT_DIR=${HOME}/development/mbe
DROPLET_DIR=${MBE_ROOT_DIR}/examples/droplet

# Copy over the $rem section and point charges that will be used for
# all the input files we're generating.

cp ${DROPLET_DIR}/droplet_qchem_rem_section_${CALCTYPE} .
cp ${DROPLET_DIR}/point_charges_anion.txt .
cp ${DROPLET_DIR}/point_charges_cation.txt .

# Copy over the driver script that will generate all the inputs.

cp ${DROPLET_DIR}/droplet_make_inputs_${CALCTYPE}.sh .

# We call a shell script that calls a shell script that calls a Python
# script...

# The driver that handles number of QM pairs, number of MM pairs, and
# so on only requires the path to the XYZ files.
./droplet_make_inputs_${CALCTYPE}.sh ../xyz
