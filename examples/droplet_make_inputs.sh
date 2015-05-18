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

# Vary the number of QM fragments and MM fragments for every droplet.
for n_qm in $(seq 0 1 2); do
    for n_mm in $(seq 0 1 4); do
        str="python ${HOME}/development/mbe/examples/droplet.py --write-input-sections-qchem --num-closest-pairs-qm=${n_qm} --num-closest-pairs-mm=${n_mm} --point-charge-output-cation=${pc_output_cation} --point-charge-output-anion=${pc_output_anion} --make-supersystem --path=${1}"
        echo ${str}
        eval ${str}
    done
    str="python ${HOME}/development/mbe/examples/droplet.py --write-input-sections-qchem --num-closest-pairs-qm=${n_qm} --all-other-pairs-mm --point-charge-output-cation=${pc_output_cation} --point-charge-output-anion=${pc_output_anion} --make-supersystem --path=${1}"
    echo ${str}
    eval ${str}
done

# Append the $rem section to every input and change the file ending to
# '.in'.
genfiles=$(find ${PWD} -type f -name "drop_*qm_*mm")
for genfile in ${genfiles[@]}; do
    genfilebase=$(basename ${genfile})
    cat ${remfile} >> ${genfile}
    mv ${genfile} ${genfile//mm/mm.in}
done
