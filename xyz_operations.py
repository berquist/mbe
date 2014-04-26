#!/usr/bin/env python

def read_xyz(filename):
    """
    A simple example for reading from an XYZ file line by line, not all at once
    using readlines().
    """
    atoms = []
    coords = []
    with open(filename, 'r') as xyzfile:
        natoms = int(xyzfile.readline().strip())
        comment = xyzfile.readline().strip()
        for line in xyzfile:
            if len(line) > 1:
                atoms.append(line.split()[0])
                coords.append(map(float, line.split()[1:]))

    return (natoms, comment, atoms, coords)

def read_fragment_xyz(filename):
    """
    A more complicated example, reading an XYZ file divided into fragments,
    again doing so incrementally.
    """
    comments = []
    frag_charges = []
    frag_multiplicities = []
    atoms = []
    coords = []
    frag_atoms = []
    frag_coords = []
    with open(filename, 'r') as fragfile:
        # the very first line is the system's total charge and multiplicity
        sys_charge, sys_multiplicity = map(int, fragfile.readline().split())
        for line in fragfile:
            # each fragment section is separated by '--' as a delimiter
            # if we find it, we're at the start of a new block
            if line[0:2] == '--':
                comment = " ".join(line.split()[1:])
                comments.append(comment)
                # the fragment's own charge and multiplicity must come immediately afterwards
                frag_charge, frag_multiplicity = map(int, next(fragfile).split())
                frag_charges.append(frag_charge)
                frag_multiplicities.append(frag_multiplicity)
                # reset the fragment atoms and coordinates after append
                if len(frag_atoms) > 0:
                    atoms.append(frag_atoms)
                if len(frag_coords) > 0:
                    coords.append(frag_coords)
                frag_atoms = []
                frag_coords = []
            # if we aren't dealing with the two fragment header lines, it
            # must be a fragment's atom entry
            else:
                atom = line.split()[0]
                atom_coords = map(float, line.split()[1:])
                frag_atoms.append(atom)
                frag_coords.append(atom_coords)
        # these need to be here otherwise the last fragment will never be appended
        atoms.append(frag_atoms)
        coords.append(frag_coords)

    return (sys_charge, sys_multiplicity, frag_charges, frag_multiplicities, atoms, coords, comments)

def write_individual_fragments(filename):
    """
    Read a combined fragment XYZ file, split it up into individual fragment XYZ
    files, and write them to disk.
    """
    sys_charge, sys_multiplicity, frag_charges, frag_multiplicities, atoms, coords, comments = read_fragment_xyz(filename)
    # the number of fragments is the length of any of the returned lists
    nfragments = len(frag_charges)
    # template for individual fragment filename(s)
    s = os.path.splitext(filename)[0] + '.F{}.xyz'
    # generate a single XYZ file for each fragment
    for f in range(nfragments):
        with open(s.format(f), 'w') as xyzfile:
            xyzfile.write('{}\n'.format(len(atoms[f])))
            xyzfile.write('{} {} {}\n'.format(frag_charges[f], frag_multiplicities[f], comments[f]))
            for atom in range(len(atoms[f])):
                xyzfile.write('{:3} {:12.7f} {:12.7f} {:12.7f}\n'.format(atoms[f][atom], *coords[f][atom]))

if __name__ == '__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('fraginpfile')
    args = parser.parse_args()
    fraginpfile = args.fraginpfile

    write_individual_fragments(fraginpfile)
