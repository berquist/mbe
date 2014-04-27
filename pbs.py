#!/usr/bin/env python

import os
import subprocess

import MBE
import orca_parser


def submit_fragment_job(fragment):
    """
    Run an ORCA calculation on the given fragment through PBS.
    """
    fragment.write(fragment.fxyz)
    stub = os.path.splitext(os.path.basename(fragment.fxyz))[0]
    eprname = os.path.splitext(fragment.fxyz)[0] + ".inp"
    pbsname = os.path.splitext(fragment.fxyz)[0] + ".pbs"
    eprhandle = open(eprname, "w")
    pbshandle = open(pbsname, "w")
    eprhandle.write(MBE.file_templates.eprfile(fragment.charge,
                                               fragment.multiplicity,
                                               stub))
    pbshandle.write(MBE.file_templates.pbsfile(stub))
    eprhandle.close()
    pbshandle.close()
    os.chdir(os.path.dirname(pbsname))
    subprocess.call(["qsub", pbsname])
    eproutname = os.path.splitext(eprname)[0] + ".out"
    orcafile = orca_parser.ORCAOutputParser(eproutname)
    return orcafile

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('nodefilename')
    parser.add_argument('--nodes', default=1)
    parser.add_argument('--ppn', default=1)
    args = parser.parse_args()
    nodefilename = args.nodefilename

    with open(nodefilename, 'r') as handle:
        nodefile_raw = handle.readlines()

    nodefile = []
    for line in nodefile_raw:
        nodefile.append(line.strip())

    
