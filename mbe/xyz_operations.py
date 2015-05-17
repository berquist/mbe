"""mbe.xyz_operations: Functions for dealing with XYZ (*.xyz) files.
(reading, writing, etc.)"""

from __future__ import print_function

import os.path

import mbe


def read_xyz(filename):
    """A simple example for reading from an XYZ file line by line, not
    all at once using readlines().
    """
    atoms = []
    coords = []
    with open(filename) as xyzfile:
        natoms = int(xyzfile.readline().strip())
        comment = xyzfile.readline().strip()
        for line in xyzfile:
            if len(line) > 1:
                atoms.append(line.split()[0])
                coords.append(list(map(float, line.split()[1:])))

    return (natoms, comment, atoms, coords)


def read_fragment_xyz(filename):
    """A more complicated example, reading an XYZ file divided into
    fragments, again doing so incrementally.
    """
    comments = []
    frag_charges = []
    frag_multiplicities = []
    atoms = []
    coords = []
    frag_atoms = []
    frag_coords = []
    atom_count = 0
    with open(filename) as fragfile:
        # the very first line is the system's total charge and multiplicity
        sys_charge, sys_multiplicity \
            = list(map(int, fragfile.readline().split()))
        for line in fragfile:
            # each fragment section is separated by '--' as a delimiter
            # if we find it, we're at the start of a new block
            if line[0:2] == '--':
                comment = " ".join(line.split()[1:])
                comments.append(comment)
                # the fragment's own charge and multiplicity
                # must come immediately afterwards
                frag_charge, frag_multiplicity \
                    = list(map(int, next(fragfile).split()))
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
                atom_coords = list(map(float, line.split()[1:]))
                frag_atoms.append(atom)
                frag_coords.append(atom_coords)
                atom_count += 1
        # these need to be here otherwise the last fragment
        # will never be appended
        atoms.append(frag_atoms)
        coords.append(frag_coords)

    return (sys_charge, sys_multiplicity, frag_charges, frag_multiplicities,
            atoms, coords, comments, atom_count)


def write_fragment_section_qchem(fragments, filename=None, stdout=True, bookends=False):
    """From an iterable of Fragment()s, write a fragment file to disk that
    can either be read by these scripts or used as a Q-Chem input.
    """

    # Make a supersystem fragment so that we can get the total charge
    # and multiplicity.
    superfrag = mbe.fragment.combine_fragment_sequence(fragments)

    blocks = []

    if bookends:
        blocks.append('$molecule')

    blocks.append('{} {}'.format(superfrag.charge, superfrag.multiplicity))

    satemp = '{:3} {:15.10f} {:15.10f} {:15.10f}'

    for fragment in fragments:
        blocks.append('-- {}'.format(fragment.comment))
        blocks.append('{} {}'.format(fragment.charge, fragment.multiplicity))
        for atomsym, atomcoords in zip(fragment.atoms, fragment.coords):
            blocks.append(satemp.format(atomsym, *atomcoords))

    if bookends:
        blocks.append('$end')

    section = '\n'.join(blocks)

    if not filename:
        if stdout:
            print(section)
    else:
        with open(filename, 'w') as outfile:
            print(section, file=outfile)

    return section


def write_fragment_section_psi(fragments, filename=None, stdout=True):
    """From an iterable of Fragment()s, write a fragment file to disk that
    can be used as a Psi input.

    The primary difference between Psi and Q-Chem fragment inputs is
    that Psi doesn't explicitly require the charge and multiplicity of
    the supersystem.
    """

    # ...so we don't need to build the supersystem.

    blocks = []

    satemp = '{:3} {:15.10f} {:15.10f} {:15.10f}'

    for fragment in fragments:
        blocks.append('--')
        blocks.append('{} {}'.format(fragment.charge, fragment.multiplicity))
        for atomsym, atomcoords in zip(fragment.atoms, fragment.coords):
            blocks.append(satemp.format(atomsym, *atomcoords))
    del blocks[0]

    if not filename:
        if stdout:
            print(section)
    else:
        with open(filename, 'w') as outfile:
            print(section, file=outfile)

    return section


def make_pointcharge_section_qchem(fragments, bookends=False):
    """If any fragments have point charges attached to them, placed at the
    atomic centers, form the point charge contribution section for a
    Q-Chem input.
    """

    lines = []

    if bookends:
        lines.append('$external_charges')

    t = ' {:f} {:f} {:f} {:f}'.format

    for fragment in fragments:
        if hasattr(fragment, 'pointcharges'):
            for (x, y, z), c in zip(fragment.coords, fragment.pointcharges):
                line = t(x, y, z, c)
                lines.append(line)

    if bookends:
        lines.append('$end')

    return '\n'.join(lines)


def write_input_sections_qchem(fragments_qm, fragments_mm, filename=None, stdout=True):
    """Given two sets of fragments, one meant to be treated quantum
    mechanically and other meant to be represented as point charges, write
    a combined $molecule/$external_charges part of a Q-Chem input.
    """

    section_molecule = write_fragment_section_qchem(fragments_qm, stdout=False, bookends=True)
    section_external_charges = make_pointcharge_section_qchem(fragments_mm, bookends=True)

    section = '\n'.join([
        section_molecule,
        '',
        section_external_charges
    ])

    if not filename:
        if stdout:
            print(section)
    else:
        with open(filename, 'w') as outfile:
            print(section, file=outfile)

    return section


def write_individual_fragments(filename):
    """Read a combined fragment XYZ file, split it up into individual
    fragment XYZ files, and write them to disk.
    """
    sys_charge, sys_multiplicity, frag_charges, frag_multiplicities, \
        atoms, coords, comments, atom_count = read_fragment_xyz(filename)
    # the number of fragments is the length of any of the returned lists
    nfragments = len(frag_charges)
    # string templates
    sftemp = os.path.splitext(filename)[0] + '_F{}.xyz'
    sctemp = '{} {} {}\n'
    satemp = '{:3} {:15.10f} {:15.10f} {:15.10f}\n'
    # generate a single XYZ file for each fragment
    for fid in range(nfragments):
        with open(sftemp.format(fid), 'w') as xyzfile:
            xyzfile.write('{}\n'.format(len(atoms[fid])))
            xyzfile.write(sctemp.format(frag_charges[fid],
                                        frag_multiplicities[fid],
                                        comments[fid]))
            for atom in range(len(atoms[fid])):
                xyzfile.write(satemp.format(atoms[fid][atom],
                                            *coords[fid][atom]))


def write_full_system(filename):
    """Read a combined fragment XYZ file and write a normal XYZ file of
    the full system.
    """
    sys_charge, sys_multiplicity, frag_charges, frag_multiplicities, \
        atoms, coords, comments, atom_count = read_fragment_xyz(filename)
    nfragments = len(frag_charges)
    sftemp = os.path.splitext(filename)[0] + '_FULL.xyz'
    satemp = '{:3} {:15.10f} {:15.10f} {:15.10f}\n'
    with open(sftemp, 'w') as xyzfile:
        xyzfile.write('{}\n\n'.format(atom_count))
        for fid in range(nfragments):
            for atom in range(len(atoms[fid])):
                xyzfile.write(satemp.format(atoms[fid][atom],
                                            *coords[fid][atom]))

def main():
    """If this file is called as a script, read in a fragment-style input
    file, and write the monomers and supersystem to disk as XYZ files.
    """

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('fraginpfile', help='name of the fragment-style input file')
    args = parser.parse_args()
    fraginpfile = args.fraginpfile

    write_individual_fragments(fraginpfile)
    write_full_system(fraginpfile)


if __name__ == '__main__':
    main()
