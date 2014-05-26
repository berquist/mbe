"""mbe.file_templates: (Default) file templates for dealing with PBS,
ORCA, and Q-Chem."""

def eprfile(charge, multiplicity, xyzfile):
    """
    A default template for an EPR input file.
    """
    return """! uks pbe0 def2-tzvpp def2-tzvpp/jk ri rijk pmodel somf(1x) noautostart tightscf grid5

%pal
 nprocs 1
 end

* xyzfile {0} {1} {2}.xyz *

%eprnmr
 tol 1e-10
 gtensor 1
 ori centerofelcharge
 printlevel 5
 end

""".format(charge, multiplicity, xyzfile)


def pbsfile(xyzfile):
    """
    A default template for a PBS job file.
    """
    return """#!/bin/bash

#PBS -N {0}
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

cp $PBS_O_WORKDIR/{0}.inp $LOCAL
cp $PBS_O_WORKDIR/{0}.xyz $LOCAL
cd $LOCAL

run_on_exit() {{
    set -v
    cp $LOCAL/* $PBS_O_WORKDIR
}}

trap run_on_exit EXIT

$(which orca) {0}.inp >& $PBS_O_WORKDIR/{0}.out
""".format(xyzfile)
