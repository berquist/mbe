"""mbe.pbs: Functions for submitting calculations to the Portable
Batch Scheduler job scheduling system."""

import os
import subprocess

import mbe
import orca_parser


def submit_fragment_job(fragment):
    """Run an ORCA calculation on the given fragment through PBS.
    """
    # The fragment coordinates need to be written to disk so the PBS
    # job submission doesn't fail (the XYZ file gets copied to the
    # node's local scratch).
    fragment.write(fragment.fxyz)
    stub = os.path.splitext(os.path.basename(fragment.fxyz))[0]
    eprname = os.path.splitext(fragment.fxyz)[0] + '.in'
    pbsname = os.path.splitext(fragment.fxyz)[0] + '.pbs'
    eprhandle = open(eprname, 'w')
    pbshandle = open(pbsname, 'w')
    eprhandle.write(mbe.file_templates.eprfile(fragment.charge,
                                               fragment.multiplicity,
                                               stub))
    pbshandle.write(mbe.file_templates.pbsfile(stub))
    eprhandle.close()
    pbshandle.close()

    os.chdir(os.path.dirname(pbsname))
    # subprocess.call(['qsub', pbsname])
    eproutname = os.path.splitext(eprname)[0] + '.out'
    # This call is here only for debugging purposes.
    subprocess.check_call(' '.join(['orca', eprname, '|', 'tee', eproutname]), shell=True)
    orcafile = orca_parser.ORCAOutputParser(eproutname)
    return orcafile


def main():
    """If this file is called as a script, ... (does nothing right now)"""

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('nodefilename', help='path to the nodefile generated by PBS')
    parser.add_argument('--nodes', default=1, help='# of nodes to run on')
    parser.add_argument('--ppn', default=1, help='# of cores to use per node')
    args = parser.parse_args()
    nodefilename = args.nodefilename

    with open(nodefilename) as handle:
        nodefile_raw = handle.readlines()

    nodefile = []
    for line in nodefile_raw:
        nodefile.append(line.strip())
    # Eventually, logic for determining which nodes/cores to run on
    # goes here.


if __name__ == '__main__':
    main()
