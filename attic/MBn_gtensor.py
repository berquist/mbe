#!/usr/bin/env python

import os.path
import sympy

from mbe.fragment import Fragment

"""MBn_gtensor.py: Given a set of XYZ files that are individual
fragments of a system, do a many-body expansion of the g-tensor
and compare to the full system."""

def is_active(f):
    """
    Determine of the combination of fragments in the given iterable
    are EPR-active (doublet, not singlet).
    """
    active = False
    for fragment in f:
        if fragment.is_active():
            active = not active
    return active


def combine_write(x, y, fxyz):
    """
    Combine two fragments into one, returning a new fragment
    and writing it to disk.
    """
    new = combine(x, y)
    new.write(fxyz)
    return new

def generate_monomers(fxyzlist):
    """
    Generate all the starting fragments (monomers) from coordinate files.
    """
    monomers = []
    for fxyz in fxyzlist:
        monomer = Fragment()
        monomer.read(fxyz)
        monomers.append(monomer)
    return monomers

def generate_dimers(m):
    """
    Generate all possible dimers from the list of monomers.
    """
    dimers = []
    n = len(m)
    for i in range(n):
        for j in range(i+1, n):
            dimers.append(combine([m[i], m[j]]))
    return dimers

def generate_trimers(m):
    """
    Generate all possible trimers from the list of monomers.
    """
    trimers = []
    n = len(m)
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                trimers.append(combine([m[i], m[j], m[k]]))
    return trimers

def generate_symbol_frgm_map(m):
    """
    Given a list of monomers, create a map from the fragment symbols to the
    Fragment objects they represent.
    """
    newdict = dict()
    for fragment in m:
        newdict[fragment.symbol_repr[0]] = fragment
    return newdict

if __name__ == '__main__':
    import argparse

    # parser = argparse.ArgumentParser()
    # parser.add_argument('fxyz', nargs='+')
    # args = parser.parse_args()
    # args = parser.parse_args(sp.check_output(["ls", "~/calc.epr/mbe/frag\*.xyz"]))
    # fxyzlist = args.fxyz

    # generate the jobs for each of the individual fragments
    # for fxyz in fxyzlist:
    #     dir = os.path.dirname(fxyz)
    #     stub = os.path.splitext(os.path.basename(fxyz))[0]
    #     orcaname = stub + ".inp"
    #     pbsname = stub + ".pbs"
    #     orcahandle = open(os.path.join(dir, orcaname), "w")
    #     pbshandle = open(os.path.join(dir, pbsname), "w")
    #     print orcaname
    #     print >> orcahandle, eprfile(charge, multiplicity, stub)
    #     print pbsname
    #     print >> pbshandle, pbsfile(stub)
    #     orcahandle.close()
    #     pbshandle.close()

    monomers = generate_monomers(fxyzlist)
    dimers = generate_dimers(monomers)
    trimers = generate_trimers(monomers)
