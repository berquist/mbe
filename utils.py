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


def gen_dict_symbol2fragment(fragments):
    """
    Given a sequence of fragments, create a map of fragment symbols
    to fragment objects.
    """
    newdict = dict()
    for fragment in fragments:
        newdict[next(iter(fragment.symbol_repr))] = fragment
    return newdict


def gen_dict_fragment2symbol(fragments):
    """
    Given a sequence of fragments, create a map of fragment objects
    to fragment symbols.
    """
    newdict = gen_dict_symbol2fragment(fragments)
    inewdict = invert_map(newdict)
    return inewdict


def gen_term_from_symlist(seq):
    """
    Given a sequence of symbols, generate a SymPy Symbol/Mul from them.
    """
    # We can't index into sets, so...
    seq = list(seq)
    if len(seq) == 1:
        return seq[0]
    else:
        return seq[0] * gen_term_from_symlist(seq[1:])
