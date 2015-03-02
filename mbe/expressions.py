"""mbe.expressions: Functions that generate SymPy expressions
corresponding to arbitrary-order many-body expansions."""

from sympy import prod
from itertools import combinations


def dEn(f, n):
    """Generate a single n-body term."""
    if n < 1:
        return 0
    elif n == 1:
        return sum(f)
    else:
        acc = prod(f)
        for i in range(1, n):
            acc -= sum(dEn(g, i) for g in combinations(f, i))
        return acc


def En(f, n):
    """Generate all possible n-body terms."""
    return sum(dEn(g, n) for g in combinations(f, n))


def MBEn(f, n):
    """Generate the full nth-order many-body expansion expression."""
    if n < 1:
        return 0
    elif n == 1:
        return En(f, n)
    else:
        return En(f, n) + MBEn(f, n - 1)


def main():
    """If this file is run as a script, ... (does nothing right now)"""

    from sympy import symbols

    f2 = symbols('a, b')
    f3 = symbols('a, b, c')
    f4 = symbols('a, b, c, d')
    f5 = symbols('a, b, c, d, e')
    f6 = symbols('a, b, c, d, e, f')
    f7 = symbols('a, b, c, d, e, f, g')


if __name__ == '__main__':
    main()
