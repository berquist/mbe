#!/usr/bin/env bash

# droplet_make_inputs.sh: Make inputs for droplet (CO2 in ionic
# liquid) calculations, varying the number of ionic liquid pairs that
# are treated using QM and MM.

# The main directory for the MBE code might need to be changed for
# your particular machine.
MBE_ROOT_DIR="${HOME}/development/mbe"
DROPLET_DIR="${MBE_ROOT_DIR}/examples/droplet"

# These files contain the point charges that will be used for the
# cation and anion; the format is such that they can be copied
# directly from a Q-Chem output (Mulliken, Hirshfeld, ChElPG,
# Merz-Kollman...).
pc_output_anion="point_charges_anion.txt"
pc_output_cation="point_charges_cation.txt"

# Vary the number of QM fragments and MM fragments for every droplet.

n_mm_1=$(seq 0 2 16)
n_mm_2=(32 64 128)
n_mm_arr=( ${n_mm_1[@]} ${n_mm_2[@]} )

for n_qm in $(seq 0 1 3); do
    for n_mm in ${n_mm_arr[@]}; do
        str="python ${DROPLET_DIR}/droplet.py --write-input-sections-qchem --num-closest-pairs-qm=${n_qm} --num-closest-pairs-mm=${n_mm} --point-charge-output-cation=${pc_output_cation} --point-charge-output-anion=${pc_output_anion} --make-supersystem --path=${1}"
        echo ${str}
        eval ${str}
    done
    str="python ${DROPLET_DIR}/droplet.py --write-input-sections-qchem --num-closest-pairs-qm=${n_qm} --all-other-pairs-mm --point-charge-output-cation=${pc_output_cation} --point-charge-output-anion=${pc_output_anion} --make-supersystem --path=${1}"
    echo ${str}
    eval ${str}
done
