#!/usr/bin/env python

"""mbe.utils: Utility functions that don't fit anywhere else."""

def term2str(term):
    """
    Convert a term (sympy.core.mul.Mul) into its string representation
    by splitting and joining at each multiplication.
    """
    s = str(term).replace('*', '')
    return s


def gen_dict_term2str(termdict):
    """
    Given a dictionary where all the keys are terms, return a new dictionary
    that maps each term to its custom string form.
    """
    newdict = dict()
    for term in termdict:
        newdict[term] = term2str(term)
    return newdict


def invert_map(m):
    """
    Return the map m, inverted.
    """
    im = {v: k for k, v in m.items()}
    return im


def gen_dict_symbol2monomer(monomers):
    """
    Given a sequence of monomers, create a map of monomer symbols
    to monomer objects.
    """
    newdict = dict()
    for monomer in monomers:
        newdict[next(iter(monomer.symbol_repr))] = monomer
    return newdict


