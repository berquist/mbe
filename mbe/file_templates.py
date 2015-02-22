"""mbe.file_templates: (Default) file templates for dealing with PBS,
ORCA, and Q-Chem."""

def eprfile(charge, multiplicity, xyzfile):
    """
    A default template for an EPR input file.
    """
    return """! uks pbe0 def2-sv(p) def2-svp/jk ri rijk pmodel somf(1x) noautostart tightscf grid5

%pal
 nprocs 1
 end

* xyzfile {charge} {multiplicity} {xyzfile}.xyz *

%eprnmr
 tol 1e-10
 gtensor 1
 ori centerofelcharge
 printlevel 5
 end

""".format(charge=charge, multiplicity=multiplicity, xyzfile=xyzfile)


def pbsfile(xyzfile):
    """
    A default template for a PBS job file.
    """
    return """#!/bin/bash

#PBS -N {xyzfile}
#PBS -q shared
#PBS -l nodes=1:ppn=1
#PBS -l walltime=144:00:00
#PBS -j oe
#PBS -l qos=low
#PBS -m abe
#PBS -M erb74@pitt.edu

module purge
module load intel/2013.0
module load openmpi/1.6.5-intel12
module load orca/3.0.1

cp $PBS_O_WORKDIR/{xyzfile}.in $LOCAL
cp $PBS_O_WORKDIR/{xyzfile}.xyz $LOCAL
cd $LOCAL

run_on_exit() {{
    set -v
    cp $LOCAL/* $PBS_O_WORKDIR
}}

trap run_on_exit EXIT

$(which orca) {xyzfile}.in >& $PBS_O_WORKDIR/{xyzfile}.out
""".format(xyzfile=xyzfile)
