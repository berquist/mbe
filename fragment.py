#!/usr/bin/env python

import os
import sympy

from MBE.xyz_operations import read_fragment_xyz


class Fragment(object):
    """
    A fragment is a molecule that can be combined with another
    fragment.
    """
    def __init__(self):
        self.fxyz = ""
        self.name = ""
        self.comment = ""
        self.charge = 0
        self.multiplicity = 1
        self.atoms = []
        self.coords = []
        self.nfragments = 0
        self.active = False
        self.symbol_repr = []

    def is_active(self):
        """
        Is this fragment going to be EPR-active? Determine based
        on the multiplicity.
        """
        if self.multiplicity == 1:
            return False
        elif self.multiplicity == 2:
            return True
        else:
            raise Exception('Only singlets and doublets are allowed.')

    def read(self, fxyz):
        """
        Load the charge, multiplicity, and coordinates into the fragment.
        """
        self.fxyz = fxyz
        self.name = os.path.splitext(os.path.basename(fxyz))[0]
        xyzhandle = open(fxyz, "r")
        xyzhandle.readline()
        self.charge, self.multiplicity \
            = list(map(int, xyzhandle.readline().split()))
        raw_coords = xyzhandle.readlines()
        xyzhandle.close()
        for line in raw_coords:
            if line != '\n':
                self.atoms.append(line.split()[0])
                self.coords.append(list(map(float, line.split()[1:])))
        self.nfragments += 1
        self.active = self.is_active()
        self.symbol_repr.append(sympy.S(self.name))

    def write(self, fxyz):
        """
        Write an XYZ file to disk, with the charge and multiplicity stored
        in the comment line.
        """
        stemp = "{:2s} {:f} {:f} {:f}\n"
        with open(fxyz, "w") as xyzhandle:
            for atom, coords in zip(self.atoms, self.coords):
                xyzhandle.write(stemp.format(atom, *coords))

    def __str__(self):
        return self.name

    def __add__(self, other):
        return self.combine(other)

    def combine(self, other):
        """
        Combine this fragment with another fragment, returning a new fragment.
        """
        new = Fragment()
        new.name = other.name + self.name
        new.charge = other.charge + self.charge
        if other.multiplicity == self.multiplicity:
            new.multiplicity = 1
        else:
            new.multiplicity = 2
        new.atoms = other.atoms + self.atoms
        new.coords = other.coords + self.coords
        new.nfragments = other.nfragments + self.nfragments
        new.active = other.active ^ self.active
        new.symbol_repr = other.symbol_repr + self.symbol_repr
        return new

    def combine_write(self, other, fxyz):
        """
        Combine this fragment with another fragment, returning a new fragment
        and writing it to disk.
        """
        new = self.combine(other)
        new.write(fxyz)
        return new


def generate_fragment_objects(filename):
    """
    Generate a list of Fragment() objects given a fragment-style
    XYZ file (Q-Chem, Psi4).
    """
    sys_charge, sys_multiplicity, frag_charges, frag_multiplicities, \
        atoms, coords, comments, atom_count = read_fragment_xyz(filename)

    sntemp = os.path.splitext(os.path.basename(filename))[0] + '_F{}'
    sftemp = os.path.join(os.path.dirname(filename), sntemp + '.xyz')

    fragments = []
    n = len(frag_charges)
    for fid in range(n):
        fragment = Fragment()
        fragment.fxyz = sftemp.format(fid)
        fragment.name = sntemp.format(fid)
        fragment.comment = comments[fid]
        fragment.charge = frag_charges[fid]
        fragment.multiplicities = frag_multiplicities[fid]
        fragment.atoms = atoms[fid]
        fragment.coords = coords[fid]
        fragment.nfragments += 1
        fragment.active = fragment.is_active()
        fragment.symbol_repr.append(sympy.S(fragment.name))
        fragments.append(fragment)

    return fragments
