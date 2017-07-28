"""mbe.epr: Functions useful specifically for EPR calculations.
"""

import numpy as np

import mbe
from .fragment import generate_fragment_from_term


def remove_inactive_terms(expr, fragments):
    """
    Given a SymPy expression that represents combinations of individual
    fragments or the union of fragments, remove those terms that would be
    EPR silent (multiplicity == 1), and return the new expression.
    """
    terms = expr.as_coefficients_dict()
    termdict = mbe.utils.gen_dict_symbol2fragment(fragments)
    fragments = []
    # Convert every term in the expression into a Fragment().
    # Perform the check.
    for term in terms:
        fragment = generate_fragment_from_term(term, termdict)
        if fragment.is_active():
            fragments.append(fragment)
    # Reconstruct the expression from the fragments with the
    # proper coefficients.
    new_expr = 0
    for fragment in fragments:
        symbol_repr = fragment.symbol_repr
        term = mbe.utils.gen_term_from_symlist(symbol_repr)
        coeff = terms[term]
        new_expr += (coeff * term)
    return new_expr

def calculate_gtensor(expr, map_frag_results):
    """
    ...
    """
    tot_gmat = 0
    map_symbols_frag = mbe.utils.gen_dict_symbol2fragment(map_frag_results.keys())
    exact_term_tuples = expr.as_terms()[0]
    for exact_term in exact_term_tuples:
        term = exact_term[0].as_coeff_Mul()
        coeff, symbol = term
        gmat = map_frag_results[map_symbols_frag[symbol]].molecule.gtensor.gmatrix
        tot_gmat += (coeff * gmat)
    tot_gmat = np.sqrt(np.linalg.eigvals(np.dot(tot_gmat.T, tot_gmat)).real)
    return tot_gmat
