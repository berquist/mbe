#!/usr/bin/env bash

# This file holds the Q-Chem $rem section that will be applied to
# every input.
remfile="droplet_qchem_rem_section"

# These are the Q-Chem outputs Mulliken and CHELPG charges that will
# be used as point charges.
pc_output_anion="${HOME}/calc.sgr/droplets/old/drop_318_qchem/anion_chelpg.out"
pc_output_cation="${HOME}/calc.sgr/droplets/old/drop_318_qchem/cation_chelpg.out"

# Find all the droplet XYZ files.
xyzfiles=$(find ${HOME}/calc.sgr/droplets/new/xyz -type f -name "drop_*.xyz")

# A big loop that varies the number of QM fragments and MM fragments
# for every droplet.
for n_qm in $(seq 0 2); do
    for n_mm in $(seq 0 4); do
        for xyzfile in ${xyzfiles[@]}; do
            xyzfilebase=$(basename ${xyzfile})
            inputfile="${xyzfilebase//.xyz/}_${n_qm}qm_${n_mm}mm.in"
            cat ${remfile} > ${inputfile}
            str="python ${HOME}/development/mbe/examples/droplet.py --print-input-sections-qchem --single=${xyzfile} --num-closest-pairs-qm=${n_qm} --num-closest-pairs-mm=${n_mm} --point-charge-output-cation=${pc_output_cation} --point-charge-output-anion=${pc_output_anion} --point-charge-type=chelpg"
            echo ${str}
            eval ${str} >> ${inputfile}
        done
    done
done
