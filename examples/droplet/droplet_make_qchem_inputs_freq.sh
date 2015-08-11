#!/usr/bin/env bash

## This script should be copied into the directory where all of the
## input files are going to be generated!

# What kind of $molecule section do we need?
#  joptype = freq          => supersystem
#  jobtype = freq (w/o CT) => fragments
#  jobtype = eda           => fragments
MOLECULE_SECTION_TYPE=supersystem
CALC_TYPE=freq

# The *root* directory for the MBE code.
MBE_ROOT_DIR="${HOME}/development/mbe"
DROPLET_DIR="${MBE_ROOT_DIR}/examples/droplet"

# Copy over the $rem section that will be used for all the input files
# we're generating.

remfile="droplet_qchem_rem_section_${CALC_TYPE}"
cp "${DROPLET_DIR}/${remfile}" .

# Create input files by appending the $rem section to files that
# already contain $molecule/$external_charges sections.

genfiles="$(find ${PWD}/../qchem_molecule_external_charges_stubs_${MOLECULE_SECTION_TYPE} -type f -name "drop_*qm_*mm")"
for genfile in ${genfiles[@]}; do
    genfilebase=$(basename "${genfile}")
    inputfile="${genfilebase//mm/mm_${CALC_TYPE}.in}"
    cp "${genfile}" "${inputfile}"
    cat ${remfile} >> "${inputfile}"
done
