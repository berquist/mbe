#!/usr/bin/env bash

# droplet_make_inputs.sh: Make inputs for droplet (CO2 in ionic
# liquid) calculations, varying the number of ionic liquid pairs that
# are treated using QM and MM.

# This file holds the Q-Chem $rem section that will be applied to
# every input.
remfile="droplet_qchem_rem_section"

# These are the Q-Chem outputs Mulliken and CHELPG charges that will
# be used as point charges.
pc_output_anion="droplet_anion_chelpg.txt"
pc_output_cation="droplet_cation_chelpg.txt"

# Find all the droplet XYZ files. Pass the path as a command-line
# argument.
xyzfiles=$(find $1 -type f -name "drop_*.xyz")

# A big loop that varies the number of QM fragments and MM fragments
# for every droplet.
for n_qm in $(seq 0 2); do
    for n_mm in $(seq 0 4); do
        for xyzfile in ${xyzfiles[@]}; do
            xyzfilebase=$(basename ${xyzfile})
            inputfile="${xyzfilebase//.xyz/}_${n_qm}qm_${n_mm}mm.in"
            cat ${remfile} > ${inputfile}
            str="python ${HOME}/development/mbe/examples/droplet.py --print-input-sections-qchem --num-closest-pairs-qm=${n_qm} --num-closest-pairs-mm=${n_mm} --point-charge-output-cation=${pc_output_cation} --point-charge-output-anion=${pc_output_anion} --make-supersystem --single=${xyzfile}"
            echo ${str}
            eval ${str} >> ${inputfile}
        done
    done
done
