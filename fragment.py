"""mbe.fragment: a class for dealing with Fragments and other standalone
functions."""

import os
import sympy

import mbe
from mbe.xyz_operations import read_fragment_xyz


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
        self.symbol_repr = set()

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
        self.symbol_repr = set(sympy.S(self.name))

    def xyz(self):
        """Form the XYZ file representation, without the number of atoms or
        the comment.
        """
        xyz_repr = []
        satemp = "{:3s} {:12.7f} {:12.7f} {:12.7f}\n"
        for atom, coords in zip(self.atoms, self.coords):
            xyz_repr.append(satemp.format(atom, *coords))
        return "".join(xyz_repr).rstrip('\n')

    def write(self, fxyz):
        """
        Write an XYZ file to disk, with the charge and multiplicity stored
        in the comment line.
        """
        with open(fxyz, "w") as xyzhandle:
            xyzhandle.write(str(len(self.atoms)) + "\n")
            xyzhandle.write(self.comment + "\n")
            xyzhandle.write(self.xyz())

    def __str__(self):
        return self.name

    def __repr__(self):
        return 'Fragment({})'.format(self.name)

    def __add__(self, other):
        return self.combine(other)

    def combine(self, other):
        """
        Combine this fragment with another fragment, returning a new fragment.
        """
        new = Fragment()
        new.name = self.name + other.name
        new.fxyz = os.path.join(os.path.dirname(self.fxyz), new.name + ".xyz")
        new.comment = self.comment + " :: " + other.comment
        new.charge = self.charge + other.charge
        if other.multiplicity == self.multiplicity:
            new.multiplicity = 1
        else:
            new.multiplicity = 2
        new.atoms = self.atoms + other.atoms
        new.coords = self.coords + other.coords
        new.nfragments = self.nfragments + other.nfragments
        new.active = self.active ^ other.active
        new.symbol_repr = set.union(self.symbol_repr, other.symbol_repr)
        return new

    def combine_write(self, other, fxyz):
        """
        Combine this fragment with another fragment, returning a new fragment
        and writing it to disk.
        """
        new = self.combine(other)
        new.write(fxyz)
        return new


def combine_two_fragments(x, y):
    """
    Combine two fragments into one, returning a new fragment.
    """
    new = Fragment()
    new.name = x.name + y.name
    new.fxyz = os.path.join(os.path.dirname(x.fxyz), new.name + ".xyz")
    new.comment = x.comment + " :: " + y.comment
    new.charge = x.charge + y.charge
    if x.multiplicity == y.multiplicity:
        new.multiplicity = 1
    else:
        new.multiplicity = 2
    new.atoms = x.atoms + y.atoms
    new.coords = x.coords + y.coords
    new.nfragments = x.nfragments + y.nfragments
    new.active = x.active ^ y.active
    new.symbol_repr = set.union(x.symbol_repr, y.symbol_repr)
    return new


def combine_fragment_sequence(f):
    """
    Combine multiple fragments into one, returning a new fragment.
    """
    new = Fragment()
    for fragment in f:
        new.name += fragment.name
        new.comment += "{} :: ".format(fragment.comment)
        new.charge += fragment.charge
        if new.multiplicity == fragment.multiplicity:
            new.multiplicity = 1
        else:
            new.multiplicity = 2
        new.atoms += fragment.atoms
        new.coords += fragment.coords
        new.nfragments += fragment.nfragments
        new.symbol_repr = new.symbol_repr.union(fragment.symbol_repr)
    # assume that all fragments we're combining live in the same directory
    new.fxyz = os.path.join(os.path.dirname(f[0].fxyz), new.name + ".xyz")
    return new


def _normalize_name(name):
    """Replace all non-alphanumeric symbols by underscores."""
    newname = []
    for c in name:
        cl = c.lower()
        if (cl >= 'a' and cl <= 'z') or (cl >= '0' and cl <= '9'):
            newname.append(c)
        else:
            newname.append('_')
    return ''.join(newname)


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
        fragment.fxyz = _normalize_name(sftemp.format(fid))
        fragment.name = _normalize_name(sntemp.format(fid))
        fragment.comment = comments[fid]
        fragment.charge = frag_charges[fid]
        fragment.multiplicity = frag_multiplicities[fid]
        fragment.atoms = atoms[fid]
        fragment.coords = coords[fid]
        fragment.nfragments += 1
        fragment.active = fragment.is_active()
        fragment.symbol_repr.add(sympy.S(fragment.name))
        fragments.append(fragment)

    return fragments


def generate_fragment_from_term(term, fragmap):
    """
    Given a SymPy term (part of an expression) and a map from symbols
    to fragments, generate the fragment corresponding to the term.
    """
    atoms = term.atoms()
    union = list()
    for atom in atoms:
        union.append(fragmap[atom])
    fragment = combine_fragment_sequence(union)
    return fragment


def generate_term_from_fragment(fragment):
    """
    Given a fragment object, generate a SymPy term corresponding to it
    (either a Symbol or a Mul).
    """
    symbol_repr = list(fragment.symbol_repr)
    return mbe.utils.gen_term_from_symlist(symbol_repr)


def generate_fragments_from_expr(expr, monomers):
    """
    Given a full SymPy expression (multiple terms) and a sequence of monomers,
    generate a list of the proper fragments (without coefficients).
    """
    # The term dict maps the symbolic representation of fragments (with their
    # union shown as multiplication) to their coefficient in the final expression.
    fragmap = mbe.utils.gen_dict_symbol2fragment(monomers)
    termdict = dict(expr.as_coefficients_dict())
    fragments = []
    for term in termdict:
        fragment = generate_fragment_from_term(term, fragmap)
        fragments.append(fragment)
    return fragments


def generate_expr_from_fragments(fragments):
    pass
