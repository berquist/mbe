#!/usr/bin/env python


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


# def gen_dict_str2term(strdict):
#     """
#     Given a dictionary that maps terms to strings, return a new dictionary
#     of the inverse map.
#     """
#     ...
#     return newdict


def T1(f):
    """
    Generate one-body terms.
    """
    n = len(f)
    terms = 0
    for i in range(n):
        terms += f[i]
    return terms


def T2(f):
    """
    Generate two-body terms.
    """
    n = len(f)
    terms = 0
    for i in range(n):
        for j in range(i+1, n):
            terms += (f[i]*f[j]) - f[i] - f[j]
    return terms


def T3(f):
    """
    Generate three-body terms.
    """
    n = len(f)
    terms = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                terms += ((f[i]*f[j]*f[k] - f[i] - f[j] - f[k]) -
                          (f[i]*f[j] - f[i] - f[j]) -
                          (f[j]*f[k] - f[j] - f[k]) -
                          (f[k]*f[i] - f[k] - f[i]))
    return terms

if __name__ == '__main__':
    from sympy import var, I, symbols
    x, y = var('x, y')
    g = 3*I+5*x**2 + 7*y + 10*x*y

    f2 = list(symbols('frag0:2'))
    f3 = list(symbols('frag0:3'))
    f4 = list(symbols('frag0:4'))
    f5 = list(symbols('frag0:5'))
    f6 = list(symbols('frag0:6'))
    f7 = list(symbols('frag0:7'))
