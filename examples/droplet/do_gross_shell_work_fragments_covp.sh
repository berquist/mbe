#!/usr/bin/env bash

## This script should be copied into the directory where all of the
## input files are going to be generated!

# What kind of $molecule section do we need?
#  joptype = freq          => supersystem
#  jobtype = freq (w/o CT) => fragments
#  jobtype = eda           => fragments
#  jobtype = eda (w/ COVP) => fragments_covp
MOLECULE_SECTION_TYPE=fragments_covp

# The *root* directory for the MBE code.
MBE_ROOT_DIR="${HOME}/development/mbe"
DROPLET_DIR="${MBE_ROOT_DIR}/examples/droplet"

# Copy over the point charges that will be used for all the input
# files we're generating.

cp "${DROPLET_DIR}/point_charges_anion.txt" .
cp "${DROPLET_DIR}/point_charges_cation.txt" .

# Copy over the driver script that will generate all the inputs.

DRIVER_NAME="droplet_make_qchem_input_sections_${MOLECULE_SECTION_TYPE}.sh"
cp "${DROPLET_DIR}/${DRIVER_NAME}" .

# We call a shell script that calls a shell script that calls a Python
# script...

# The driver that handles number of QM pairs, number of MM pairs, and
# so on only requires the path to the XYZ files.
./${DRIVER_NAME} ../xyz
