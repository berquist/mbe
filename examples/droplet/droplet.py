#!/usr/bin/env python

"""droplet.py: Routines for working with the CO2-IL droplet XYZ
files."""

from __future__ import print_function

from glob import glob
from itertools import count
from itertools import cycle
from math import sqrt

import os
import numpy as np
from sympy import S

import mbe
from mbe.utils import pad_left_zeros
import periodic_table as pt


def rename_old_droplet_files():
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


def rename_new_droplet_files():
    """Rename all the original droplet files from 'clusterXXX.xyz' to
    'drop_XXX.xyz', where the number is left-padded with zeros.
    """

    filenames = glob('cluster*.xyz')
    print(filenames)

    # Store the maximum length of the internal number.
    maxlen = 0

    # Unfortunately, two full traversals need to be performed.
    # 1. Find the maximum length of the internal number and store it.
    for oldfilename in filenames:
        numstr = oldfilename[7:-4]
        newlen = len(numstr)
        if newlen > maxlen:
            maxlen = newlen

    newfilenames = []

    # 2. Go through each file again and rename.
    for oldfilename in filenames:
        numstr = oldfilename[7:-4]
        filenumlen = len(numstr)
        if filenumlen < maxlen:
            numstr = pad_left_zeros(numstr, maxlen)
        newfilename = 'drop_{}.xyz'.format(numstr)
        os.rename(oldfilename, newfilename)
        print(oldfilename + ' -> ' + newfilename)
        newfilenames.append(newfilename)

    return newfilenames


# def get_bond_graph(obmol):
#     """Test a few possible ways for getting the bond_connectivities using
#     Open Babel. Not used for the droplets.
#     """

#     bond_pair_indices1 = []
#     bond_pair_indices2 = []

#     for obbond in ob.OBMolBondIter(obmol):
#         bond_pair_indices1.append((obbond.GetBeginAtomIdx(),
#                                    obbond.GetEndAtomIdx(),
#                                    obbond.GetBondOrder()))

#     for obbond in ob.OBMolBondIter(obmol):
#         bond_pair_indices2.append((obbond.GetBeginAtom().GetIndex(),
#                                    obbond.GetEndAtom().GetIndex(),
#                                    obbond.GetBondOrder()))

#     bonds = [{'atoms': [bond.GetBeginAtom().GetIndex(),
#                         bond.GetEndAtom().GetIndex()],
#               'order': bond.GetBondOrder(),
#               'symbols': [pt.Element[bond.GetBeginAtom().GetAtomicNum()],
#                           pt.Element[bond.GetEndAtom().GetAtomicNum()]]}
#              for bond in ob.OBMolBondIter(obmol)]

#     print('bond pair indices (1)')
#     print(bond_pair_indices1)
#     print('bond pair indices (2)')
#     print(bond_pair_indices2)
#     print('bond pair indices (3)')
#     for bond in bonds:
#         print(bond)


# def make_obmol_from_file(filename, obconv):
#     """Make an Open Babel molecule from a file, using an OBConv instance
#     with the correct conversion type already set.
#     """

#     print(filename)
#     obmol = ob.OBMol()
#     obconv.ReadFile(obmol, filename)
#     obmol.ConnectTheDots()
#     obmol.PerceiveBondOrders()
#     return obmol


# def make_ob_fragments_from_obmol(obmol):
#     """Make a list of Open Babel molecules (fragments) from an OBMol
#     containing multiple (non-bonded) molecules.
#     """

#     return obmol.Separate()


def determine_fragment_grouping(atoms):
    """From the main list of atoms (probably right from an input file),
    return a list of lists, where each list contains the indices for atoms
    on distinct molecules.
    """

    grouping_anions = []
    grouping_cations = []
    grouping_CO2 = []

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
                grouping = list(range(start, end))
                if possible_ordering == ordering_anion:
                    grouping_anions.append(grouping)
                elif possible_ordering == ordering_cation:
                    grouping_cations.append(grouping)
                else:
                    grouping_CO2.append(grouping)
                start = end

    # A list of lists, where each list contains the indices for atoms
    # on distinct molecules.
    return grouping_anions, grouping_cations, grouping_CO2


def make_fragments_from_grouping(atoms, coords, grouping, start=1, pc=None, pc_addition='every'):
    """Given a set of atoms (list of atomic symbols), coordinates (list of
    length 3 lists), and a grouping of atoms into molecules, make a list
    of Fragment objects corresponding to the grouping.

    Optionally, deal with point charge business.
    """

    fragments = []

    for i, group in zip(count(start=start), grouping):
        start, end = group[0], group[-1] + 1
        fragment = mbe.fragment.Fragment()
        fragment.atoms = atoms[start:end]
        fragment.formula_string = fragment._make_canonical_formula_string()
        fragment.coords = coords[start:end]
        fragment.charge = map_charges[fragment.formula_string]
        fragment.nfragments = 1
        fragment.comment = ' '.join([fragment.formula_string, '({})'.format(i)])
        fragment.name = 'f{}'.format(i)
        fragment.symbol_repr = {S(fragment.name)}
        if pc is not None:
            if pc_addition == 'every':
                fragment.pointcharges = pc
            elif pc_addition == 'unique':
                fragment.pointcharges = pc[start:end]
            else:
                pass
        fragments.append(fragment)

    return fragments


def combine_fragments_for_covp(fragments):
    """COVP analysis with Q-Chem can only be performed with two
    fragments. Given a list of N fragments, combine the first N-1 into
    one, leaving the last one uncombined.
    """

    return [combine_fragment_sequence(fragments[:-1]), fragments[-1]]


def get_point_charges_qmout(pc_file_path, pc_type):
    """Use cclib to parse a QM output file for point charges (Mulliken,
    Lowdin, CHELPG, ...).
    """

    from cclib.parser import ccopen

    job = ccopen(pc_file_path)
    data = job.parse()

    return data.atomcharges[pc_type]


def get_point_charges_txt(pc_file_path):
    """Extract point charges from a plain text file, where the charges
    themselves are in a column (assumed to be the last).
    """

    pointcharges = []

    with open(pc_file_path) as pc_file:
        for line in pc_file:
            pointcharges.append(float(line.split()[-1]))

    return pointcharges


def distance_twopoint(l1, l2):
    """Return the distance between two Cartesian points (given as
    lists).
    """

    return sqrt(((l1[0] - l2[0])**2) + ((l1[1] - l2[1])**2) + ((l1[2] - l2[2])**2))


def distance_atomic(fragment1, fragment2):
    """Return the distance between the atom in fragment 1 and the atom in
    fragment 2 that are closest to each other.
    """

    shortest_distance = 99999.9

    for c1 in fragment1.coords:
        for c2 in fragment2.coords:
            current_distance = distance_twopoint(c1, c2)
            if current_distance < shortest_distance:
                shortest_distance = current_distance

    return shortest_distance


def distance_centerofmass(fragment1, fragment2):
    """Return the distance between the center of mass of fragment 1 and
    the center of mass of fragment 2."""

    # implement me!
    pass


def get_n_closest_fragments(n, target_fragment, other_fragments, method='atomic'):
    """Return the n closest fragments to the target fragment, based on
    distance between centers of mass or pairwise between atoms.
    """

    closest_fragments = []
    distances = []

    # Calculate the distance between the target fragment and the other
    # fragments.
    for fragidx, other_fragment in enumerate(other_fragments):
        if method == 'atomic':
            distance = distance_atomic(target_fragment, other_fragment)
        elif method == 'centerofmass':
            # implement me!
            pass
        else:
            raise NotImplementedError
        distances.append((fragidx, distance))

    # Now that all the distances are calculated, sort them in
    # increasing order.
    distances = sorted(distances, key=lambda x: x[1])

    # Return the n closest fragments by:
    # 1. Get the index of the ith fragment.
    # 2. Look into all other fragments for it.
    # 3. Append the Fragment object to the list.
    for i in range(n):
        closest_fragments.append(other_fragments[distances[i][0]])

    return closest_fragments


def get_n_closest_pairs(n, fragments_anions, fragments_cations, fragment_CO2, method='atomic'):
    """Return the n closest ionic liquid pairs to the CO2, given lists of
    Fragment objects of each.
    """

    closest_fragments_anions = get_n_closest_fragments(n, fragment_CO2, fragments_anions, method)
    closest_fragments_cations = get_n_closest_fragments(n, fragment_CO2, fragments_cations, method)

    return closest_fragments_anions, closest_fragments_cations


if __name__ == '__main__':
    import argparse
    import sys
    # import logging

    parser = argparse.ArgumentParser()

    parser.add_argument('--verbose',
                        action='store_true',
                        help="""Print more information to stdout.""")
    parser.add_argument('--debug',
                        action='store_true',
                        help="""Print debug-level information to stdout.""")

    parser.add_argument('--rename-old',
                        action='store_true',
                        help="""whether or not to run the rename method on the \
                        old set of droplet XYZ files (only for the original files!)""")
    parser.add_argument('--rename-new',
                        action='store_true',
                        help="""whether or not to run the rename method on the \
                        new set of droplet XYZ files (only for the original files!)""")

    parser.add_argument('--single',
                        help="""pass a single droplet XYZ file to only \
                        operate on that file""")
    parser.add_argument('--path',
                        default='.',
                        help="""the path all of the droplet coordinate files \
                        are contained (can be relative or absolute)""")

    parser.add_argument('--write-fragment-input-qchem',
                        action='store_true',
                        help="""Write the new fragments from each droplet to \
                        disk as a single fragment input file (Q-Chem style).""")
    parser.add_argument('--print-fragment-input-qchem',
                        action='store_true',
                        help="""Print the fragment input (Q-Chem style).""")
    parser.add_argument('--write-fragment-input-psi',
                        action='store_true',
                        help="""Write the new fragments from each droplet to \
                        disk as a single fragment input file (Psi style).""")
    parser.add_argument('--print-fragment-input-psi',
                        action='store_true',
                        help="""Print the fragment input (Psi style).""")

    parser.add_argument('--write-input-sections-qchem',
                        action='store_true',
                        help="""Same as --write-fragment-input-xxx, but with \
                        possible point charges (Q-Chem style).""")
    parser.add_argument('--print-input-sections-qchem',
                        action='store_true',
                        help="""Same as --print-fragment-input-xxx, but with \
                        possible point charges (Q-Chem style).""")

    parser.add_argument('--mbe-order',
                        type=int,
                        default=0,
                        help="""Order of the many-body expansion to go up to.""")

    parser.add_argument('--point-charge-type',
                        choices=('mulliken', 'lowdin', 'chelpg', 'hirshfeld'),
                        default='mulliken',
                        help="""The type of point charges to extract from \
                        an output file.""")
    parser.add_argument('--point-charge-file-type',
                        choices=('qmout', 'txt'),
                        default='txt',
                        help="""If 'qmout', the files containing point charges \
                        are QM outputs to be parsed by cclib. If 'txt', they \
                        are two columns, the first containing atomic symbols, \
                        the second containing the charge; --point-charge-type \
                        doesn't matter.""")
    parser.add_argument('--point-charge-output-cation',
                        help="""An output file containing point charges \
                        that will be applied to every cation.""")
    parser.add_argument('--point-charge-output-anion',
                        help="""An output file containing point charges \
                        that will be applied to every anion.""")
    parser.add_argument('--point-charge-output-unique',
                        help="""An output file containing unique point charges \
                        for every atom.""")

    parser.add_argument('--num-closest-pairs-qm',
                        type=int,
                        default=0,
                        help="""The number of closest ionic liquid pairs to the \
                        CO2 that will be treated quantum mechanically.""")
    parser.add_argument('--num-closest-pairs-mm',
                        type=int,
                        default=0,
                        help="""The number of closest ionic liquid pairs to the \
                        CO2 that will be treated as point charges. If any pairs \
                        are treated using QM, these are the closest outside the \
                        QM region.""")
    parser.add_argument('--distance-metric',
                        choices=('centerofmass', 'atomic'),
                        default='atomic',
                        help="""The metric for calculating distances between \
                        fragments. 'atomic' is ...""")
    parser.add_argument('--all-pairs-qm',
                        action='store_true',
                        help="""Treat all ionic liquid pairs quantum mechanically.""")
    parser.add_argument('--all-other-pairs-mm',
                        action='store_true',
                        help="""Treat all possible ionic liquid pairs as point \
                        charges. If any pairs are treated using QM, treat all \
                        others as point charges.""")

    parser.add_argument('--make-supersystem',
                        action='store_true',
                        help="""When writing or printing fragment inputs, \
                        combine all fragments into a single block, making a \
                        'traditional' molecule input. For Q-Chem, this is \
                        required when only working with a single fragment.""")
    parser.add_argument('--qchem-covp',
                        action='store_true',
                        help="""If writing a fragment (non-supersystem) input for Q-Chem, combine
                        all the fragments except the last, so COVP analysis can be
                        performed.""")

    args = parser.parse_args()

    if args.debug:
        args.verbose = True

    # if args.verbose:
    #     level = logging.INFO
    # if args.debug:
    #     level = logging.DEBUG
    # else:
    #     level = logging.WARNING
    # logging.basicConfig(level=level)

    # If we want to rename all the files, be safe and don't try any
    # other operations.
    if args.rename_old:
        rename_old_droplet_files()
        sys.exit(0)
    if args.rename_new:
        rename_new_droplet_files()
        sys.exit(0)

    if args.single:
        filenames = [args.single]
    else:
        filenames = [name for name in glob(os.path.join(args.path, 'drop_*.xyz'))]

    if args.debug:
        print(filenames)

    # obconv = ob.OBConversion()
    # obconv.SetInAndOutFormats('xyz', 'xyz')

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

        if args.verbose:
            print(basename, natoms)

        fragments = []

        # Determine which atoms from the XYZ file belong to individual
        # anions, cations, and the CO2, returning a list of lists,
        # each containing the indices of an atom.
        grouping_anions, grouping_cations, grouping_CO2 = determine_fragment_grouping(atoms)
        grouping = grouping_anions + grouping_cations + grouping_CO2
        if args.debug:
            print('grouping of anions:')
            print(grouping_anions)
            print('grouping of cations:')
            print(grouping_cations)
            print('grouping of CO2:')
            print(grouping_CO2)
            print('all groupings:')
            print(grouping)

        # Figure out whether or not we need to add point charges to
        # fragments based on if certain command-line arguments were
        # passed.
        if args.point_charge_file_type == 'qmout':
            if args.point_charge_output_unique:
                # There is a unique point charge for every atom in the
                # input XYZ file.
                pointcharges_unique = get_point_charges_qmout(args.point_charge_output_unique, args.point_charge_type)
            elif args.point_charge_output_anion and args.point_charge_output_cation:
                # There are common point charges for each anion and cation
                # in the XYZ file.
                pointcharges_anion = get_point_charges_qmout(args.point_charge_output_anion, args.point_charge_type)
                pointcharges_cation = get_point_charges_qmout(args.point_charge_output_cation, args.point_charge_type)
            else:
                # Don't worry about adding point charges.
                pass
        else:
            # Must be the simple format.
            if args.point_charge_output_unique:
                pointcharges_unique = get_point_charges_txt(args.point_charge_output_unique)
            elif args.point_charge_output_anion and args.point_charge_output_cation:
                pointcharges_anion = get_point_charges_txt(args.point_charge_output_anion)
                pointcharges_cation = get_point_charges_txt(args.point_charge_output_cation)


        # From the indices, ...
        if args.point_charge_output_unique:
            start = 1
            fragments_anions = make_fragments_from_grouping(atoms, coords, grouping_anions, start=start, pc=pointcharges_unique, pc_addition='unique')
            start = len(fragments_anions)
            fragments_cations = make_fragments_from_grouping(atoms, coords, grouping_cations, start=start, pc=pointcharges_unique, pc_addition='unique')
            start = len(fragments_anions) + len(fragments_cations)
            fragment_CO2 = make_fragments_from_grouping(atoms, coords, grouping_CO2, start=start)[0]
        elif args.point_charge_output_anion and args.point_charge_output_cation:
            start = 1
            fragments_anions = make_fragments_from_grouping(atoms, coords, grouping_anions, start=start, pc=pointcharges_anion, pc_addition='every')
            start = len(fragments_anions)
            fragments_cations = make_fragments_from_grouping(atoms, coords, grouping_cations, start=start, pc=pointcharges_cation, pc_addition='every')
            start = len(fragments_anions) + len(fragments_cations)
            fragment_CO2 = make_fragments_from_grouping(atoms, coords, grouping_CO2, start=start)[0]
        else:
            start = 1
            fragments_anions = make_fragments_from_grouping(atoms, coords, grouping_anions, start=start)
            start = len(fragments_anions)
            fragments_cations = make_fragments_from_grouping(atoms, coords, grouping_cations, start=start)
            start = len(fragments_anions) + len(fragments_cations)
            fragment_CO2 = make_fragments_from_grouping(atoms, coords, grouping_CO2, start=start)[0]

        if args.debug:
            if args.point_charge_output_unique:
                for f in fragments_anions:
                    print(f.pointcharges)
                for f in fragments_cations:
                    print(f.pointcharges)
            elif args.point_charge_output_cation and args.point_charge_output_anion:
                for f in fragments_anions:
                    print(f.pointcharges)
                for f in fragments_cations:
                    print(f.pointcharges)
            else:
                pass

        fragments = fragments_anions + fragments_cations + [fragment_CO2]

        disk_fragments_qm = []
        disk_fragments_mm = []

        # Control flow:
        # 1. all fragments are QM
        # 2. some fragments are QM
        #  remainder: all are MM
        #  remainder: some are MM
        #  remainder: ignore
        # 3. all other than CO2 are MM
        # 4. some fragments are MM
        # 5. CO2 alone
        if args.all_pairs_qm:
            disk_fragments_qm = fragments
        elif args.num_closest_pairs_qm:
            n_qm = args.num_closest_pairs_qm
            closest_pairs_qm = get_n_closest_pairs(n_qm, fragments_anions, fragments_cations, fragment_CO2, method=args.distance_metric)
            closest_anions_qm = closest_pairs_qm[0]
            closest_cations_qm = closest_pairs_qm[1]
            all_other_anions = [f for f in fragments_anions if f not in closest_anions_qm]
            all_other_cations = [f for f in fragments_cations if f not in closest_cations_qm]
            if args.all_other_pairs_mm:
                disk_fragments_mm = all_other_anions + all_other_cations
            elif args.num_closest_pairs_mm:
                n_mm = args.num_closest_pairs_mm
                closest_pairs_mm = get_n_closest_pairs(n_mm, all_other_anions, all_other_cations, fragment_CO2, method=args.distance_metric)
                closest_anions_mm = closest_pairs_mm[0]
                closest_cations_mm = closest_pairs_mm[1]
                disk_fragments_mm = closest_anions_mm + closest_cations_mm
            else:
                pass
            disk_fragments_qm = closest_anions_qm + closest_cations_qm + [fragment_CO2]
        elif args.all_other_pairs_mm:
            disk_fragments_qm = [fragment_CO2]
            disk_fragments_mm = fragments_anions + fragments_cations
        elif args.num_closest_pairs_mm:
            n_mm = args.num_closest_pairs_mm
            closest_pairs_mm = get_n_closest_pairs(n_mm, fragments_anions, fragments_cations, fragment_CO2, method=args.distance_metric)
            closest_anions_mm = closest_pairs_mm[0]
            closest_cations_mm = closest_pairs_mm[1]
            disk_fragments_mm = closest_anions_mm + closest_cations_mm
            disk_fragments_qm = [fragment_CO2]
        else:
            disk_fragments_qm = [fragment_CO2]

        # Do any of the fragments need to be recombined for any reason
        # (such as COVP analysis)?
        if args.qchem_covp:
            disk_fragments_qm = combine_fragments_for_covp(disk_fragments_qm)

        # For now, we do nothing with the many-body expansion other
        # than generate and print it out if requested.
        if args.mbe_order > 0:
            monomer_symbols = []
            for monomer in fragments:
                for symbol in monomer.symbol_repr:
                    monomer_symbols.append(symbol)
            if args.debug:
                print('Symbolic representation of monomers:')
                print(monomer_symbols)
            mbe_expression = mbe.expressions.MBEn(monomer_symbols, args.mbe_order)
            if args.debug:
                print('Full MBE{} expression:'.format(args.mbe_order))
                print(mbe_expression)

        # Write or print full Q-Chem $molecule/$external_charges sections.
        if args.write_input_sections_qchem:
            n_qm = (len(disk_fragments_qm) - 1) // 2
            n_mm = len(disk_fragments_mm) // 2
            filename = '{}_{}qm_{}mm'.format(os.path.splitext(os.path.basename(filename))[0], n_qm, n_mm)
            mbe.xyz_operations.write_input_sections_qchem(disk_fragments_qm, disk_fragments_mm, supersystem=args.make_supersystem, filename=filename, stdout=False)
        if args.print_input_sections_qchem:
            mbe.xyz_operations.write_input_sections_qchem(disk_fragments_qm, disk_fragments_mm, supersystem=args.make_supersystem)

        # Write fragments to disk or print them to stdout.
        if args.write_fragment_input_qchem:
            filename = 'frag_{}'.format(os.path.basename(filename))
            mbe.xyz_operations.write_fragment_section_qchem(disk_fragments_qm, filename=filename)
        if args.write_fragment_input_psi:
            filename = 'frag_{}'.format(os.path.basename(filename))
            mbe.xyz_operations.write_fragment_section_psi(disk_fragments_qm, filename=filename)
        if args.print_fragment_input_qchem:
            mbe.xyz_operations.write_fragment_section_qchem(disk_fragments_qm)
        if args.print_fragment_input_psi:
            mbe.xyz_operations.write_fragment_section_psi(disk_fragments_qm)
