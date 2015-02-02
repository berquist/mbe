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

    f2 = list(symbols('f0:2'))
    f3 = list(symbols('f0:3'))
    f4 = list(symbols('f0:4'))
    f5 = list(symbols('f0:5'))
    f6 = list(symbols('f0:6'))
    f7 = list(symbols('f0:7'))


if __name__ == '__main__':
    main()
