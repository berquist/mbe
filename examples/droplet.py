#!/usr/bin/env python

"""droplet.py: Routines for working with the CO2-IL droplet XYZ
files."""

from __future__ import print_function

from glob import glob
from itertools import count
from itertools import cycle

import os
import os.path
import numpy as np

try:
    import pybel as pb
    import openbabel as ob
except ImportError:
    pass

import mbe
from mbe.utils import pad_left_zeros

import periodic_table as pt


def rename_droplet_files():
    """Rename all the original droplet files from 'XXX_drop.xyz' to
    'drop_XXX.xyz', where the number is left-padded with zeros.
    """

    filenames = glob('*_drop.xyz')
    print(filenames)

    # Store the maximum length of the internal number.
    maxlen = 0

    # Unfortunately, two full traversals need to be performed.
    # 1. Find the maximum length of the internal number and store it.
    for oldfilename in filenames:
        newlen = len(oldfilename.split('_')[0])
        if newlen > maxlen:
            maxlen = newlen

    newfilenames = []

    # 2. Go through each file again and rename.
    for oldfilename in filenames:
        splitname = oldfilename.split('_')
        filenumlen = len(splitname[0])
        if filenumlen < maxlen:
            splitname[0] = pad_left_zeros(splitname[0], maxlen)
        newfilename = 'drop_{}.xyz'.format(splitname[0])
        os.rename(oldfilename, newfilename)
        print(oldfilename + ' -> ' + newfilename)
        newfilenames.append(newfilename)

    return newfilenames


def get_bond_graph(obmol):
    """Test a few possible ways for getting the bond_connectivities using
    Open Babel. Not used for the droplets.
    """

    bond_pair_indices1 = []
    bond_pair_indices2 = []

    for obbond in ob.OBMolBondIter(obmol):
        bond_pair_indices1.append((obbond.GetBeginAtomIdx(),
                                   obbond.GetEndAtomIdx(),
                                   obbond.GetBondOrder()))

    for obbond in ob.OBMolBondIter(obmol):
        bond_pair_indices2.append((obbond.GetBeginAtom().GetIndex(),
                                   obbond.GetEndAtom().GetIndex(),
                                   obbond.GetBondOrder()))

    bonds = [{'atoms': [bond.GetBeginAtom().GetIndex(),
                        bond.GetEndAtom().GetIndex()],
              'order': bond.GetBondOrder(),
              'symbols': [pt.Element[bond.GetBeginAtom().GetAtomicNum()],
                          pt.Element[bond.GetEndAtom().GetAtomicNum()]]}
             for bond in ob.OBMolBondIter(obmol)]

    print('bond pair indices (1)')
    print(bond_pair_indices1)
    print('bond pair indices (2)')
    print(bond_pair_indices2)
    print('bond pair indices (3)')
    for bond in bonds:
        print(bond)


def make_obmol_from_file(filename, obconv):
    """Make an Open Babel molecule from a file, using an OBConv instance
    with the correct conversion type already set.
    """

    print(filename)
    obmol = ob.OBMol()
    obconv.ReadFile(obmol, filename)
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()
    return obmol


def make_ob_fragments_from_obmol(obmol):
    """Make a list of Open Babel molecules (fragments) from an OBMol
    containing multiple (non-bonded) molecules.
    """

    return obmol.Separate()


def determine_fragment_grouping(atoms):
    """From the main list of atoms (probably right from an input file),
    return a list of lists, where each list contains the indices for atoms
    on distinct molecules.
    """

    grouping = []

    # We have some a priori knowledge of the ordering of atoms in each
    # of the molecules that can be used to our advantage.
    ordering_anion = ['F', 'F', 'F', 'F', 'F', 'F', 'P']
    ordering_cation = ['N', 'C', 'N', 'C', 'C', 'H', 'H', 'H', 'C', 'H', 'H',
                       'H', 'C', 'H', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'C',
                       'H', 'H', 'H']
    ordering_CO2 = ['C', 'O', 'O']

    # An iterator over each ordering that will never run out.
    possible_orderings = cycle((ordering_anion, ordering_cation, ordering_CO2))

    start = 0
    while start < len(atoms):
        # For each possible ordering (representing a single molecule),
        # make a "window" over the atom ordering from the total input,
        # seeing if they match.
        for possible_ordering in possible_orderings:
            # The end of the window always needs to be updated.
            end = start + len(possible_ordering)
            if atoms[start:end] == []:
                break
            # If we have a match, add its indices to a list as a
            # "group" representing a molecule. Only update the start
            # of the window when we have a match.
            if possible_ordering == atoms[start:end]:
                grouping.append(list(range(start, end)))
                start = end

    # A list of lists, where each list contains the indices for atoms
    # on distinct molecules.
    return grouping


if __name__ == '__main__':
    import argparse
    import sys
    import logging

    parser = argparse.ArgumentParser()

    parser.add_argument('--verbose', action='store_true', help="""Print more information to stdout.""")
    parser.add_argument('--debug', action='store_true', help="""Print debug-level information to stdout.""")
    parser.add_argument('--rename', action='store_true', help="""whether or not to run the rename method (only for the original files)""")
    parser.add_argument('--single', help="""pass a single droplet XYZ file to only operate on that file""")
    parser.add_argument('--path', default='.', help="""the path all of the droplet coordinate files are contained (can be relative or absolute""")
    parser.add_argument('--write-fragments', action='store_true', help="""Write the new fragments from each droplet to disk as individual XYZ files.""")
    parser.add_argument('--write-fragment-input-qchem', action='store_true', help="""Write the new fragments from each droplet to disk as a single fragment input file (Q-Chem style).""")
    parser.add_argument('--print-fragment-input-qchem', action='store_true', help="""Print the fragment input (Q-Chem style).""")
    parser.add_argument('--write-fragment-input-psi', action='store_true', help="""Write the new fragments from each droplet to disk as a single fragment input file (Psi style).""")
    parser.add_argument('--print-fragment-input-psi', action='store_true', help="""Print the fragment input (Psi style).""")

    args = parser.parse_args()

    # logger = logging.

    # If we want to rename all the files, be safe and don't try any
    # other operations.
    if args.rename:
        rename_droplet_files()
        sys.exit(0)

    if args.single:
        filenames = [args.single]
    else:
        filenames = [name for name in glob(os.path.join(args.path, '*.xyz'))]

    if args.debug:
        print(filenames)

    obconv = ob.OBConversion()
    obconv.SetInAndOutFormats('xyz', 'xyz')

    # These are the total charges for the each of the ionic liquid
    # components.
    map_charges = {
        'C8H15N2': +1,
        'F6P': -1,
        'CO2': 0,
    }

    for filename in filenames:

        # Here is how it might be done with Open Babel as a more general case..
        # obmol_droplet = make_obmol_from_file(filename, obconv)
        # obmol_fragments = make_ob_fragments_from_obmol(obmol_droplet)

        # fragments = []

        # for i, obmol_fragment in zip(count(start=1), obmol_fragments):

        #     pbmol_fragment = pb.Molecule(obmol_fragment)
        #     atoms = [pt.Element[atom.atomicnum] for atom in pbmol_fragment.atoms]
        #     coords = [atom.coords for atom in pbmol_fragment.atoms]
        #     formula = pbmol_fragment.formula
        #     if formula not in map_charges.keys():
        #         print(formula)
        #     if formula in map_charges.keys():
        #         charge = map_charges[formula]
        #     else:
        #         charge = -999

        #     fragment = mbe.fragment.Fragment()
        #     fragment.atoms = atoms
        #     fragment.coords = coords
        #     fragment.charge = charge
        #     fragment.nfragments = 1
        #     fragment.comment = ' '.join([formula, '({})'.format(i)])
        #     fragment.name = fragment.comment
        #     print(fragment)
        #     fragments.append(fragment)

        basename = os.path.splitext(os.path.basename(filename))[0]

        # Read in the entire XYZ file before breaking it apart into fragments.
        natoms, comment, atoms, coords = mbe.xyz_operations.read_xyz(filename)

        print(basename, natoms)

        fragments = []

        grouping = determine_fragment_grouping(atoms)

        for i, group in zip(count(start=1), grouping):
            start, end = group[0], group[-1] + 1
            fragment = mbe.fragment.Fragment()
            fragment.atoms = atoms[start:end]
            fragment.formula_string = fragment._make_canonical_formula_string()
            fragment.coords = coords[start:end]
            fragment.charge = map_charges[fragment.formula_string]
            fragment.nfragments = 1
            fragment.comment = ' '.join([fragment.formula_string, '({})'.format(i)])
            fragment.name = fragment.comment
            print(fragment)
            fragments.append(fragment)

        if args.write_fragment_input_qchem:
            mbe.xyz_operations.write_fragment_section_qchem(fragments, filename='frag_{}'.format(os.path.basename(filename)))
        if args.write_fragment_input_psi:
            mbe.xyz_operations.write_fragment_section_psi(fragments, filename='frag_{}'.format(os.path.basename(filename)))
        if args.print_fragment_input_qchem:
            mbe.xyz_operations.write_fragment_section_qchem(fragments)
        if args.print_fragment_input_psi:
            mbe.xyz_operations.write_fragment_section_psi(fragments)
