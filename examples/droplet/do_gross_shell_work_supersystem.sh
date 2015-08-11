#!/usr/bin/env bash

## This script should be copied into the directory where all of the
## input files are going to be generated!

# What kind of $molecule section do we need?
#  joptype = freq          => supersystem
#  jobtype = freq (w/o CT) => fragments
#  jobtype = eda           => fragments
MOLECULE_SECTION_TYPE=supersystem

# The *root* directory for the MBE code.
MBE_ROOT_DIR="${HOME}/development/mbe"
DROPLET_DIR="${MBE_ROOT_DIR}/examples/droplet"

# Copy over the point charges that will be used for all the input
# files we're generating.

cp "${DROPLET_DIR}/point_charges_anion.txt" .
cp "${DROPLET_DIR}/point_charges_cation.txt" .

# Copy over the driver script that will generate all the inputs.

cp "${DROPLET_DIR}/droplet_make_input_sections_${MOLECULE_SECTION_TYPE}.sh" .

# We call a shell script that calls a shell script that calls a Python
# script...

# The driver that handles number of QM pairs, number of MM pairs, and
# so on only requires the path to the XYZ files.
./droplet_make_input_sections_${MOLECULE_SECTION_TYPE}.sh ../xyz
