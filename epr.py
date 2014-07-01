"""mbe.epr: Functions useful specifically for EPR calculations.
"""

import mbe
from mbe.fragment import generate_fragment_from_term


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
