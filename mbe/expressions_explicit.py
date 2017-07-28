"""mbe.expressions_explicit: Explicit formulation of formulas in
mbe.expressions for some low orders. Useful for testing."""

from itertools import combinations


def dE1(i):
    """Generate a single one-body term."""
    return i
def dE2(i, j):
    """Generate a single two-body term."""
    return i*j - i - j
def dE3(i, j, k):
    """Generate a single three-body term."""
    return i*j*k - dE2(i, j) - dE2(i, k) - dE2(j, k) - i - j - k
def dE4(i, j, k, l):
    """Generate a single four-body term."""
    return i*j*k*l \
        - dE3(i, j, k) - dE3(j, k, l) - dE3(i, j, l) - dE3(i, k, l) \
        - dE2(i, j) - dE2(j, k) - dE2(k, l) \
        - dE2(i, k) - dE2(i, l) - dE2(j, l) \
        - i - j - k - l
def dE2gen(i, j):
    """Generate a single two-body term."""
    return i*j \
        - sum(dE1(*g) for g in combinations([i, j], 1))
def dE3gen(i, j, k):
    """Generate a single three-body term."""
    return i*j*k \
        - sum(dE2(*g) for g in combinations([i, j, k], 2)) \
        - sum(dE1(*g) for g in combinations([i, j, k], 1))
def dE4gen(i, j, k, l):
    """Generate a single four-body term."""
    return i*j*k*l \
        - sum(dE3(*g) for g in combinations([i, j, k, l], 3)) \
        - sum(dE2(*g) for g in combinations([i, j, k, l], 2)) \
        - sum(dE1(*g) for g in combinations([i, j, k, l], 1))
def dE5gen(i, j, k, l, m):
    """Generate a single five-body term."""
    return i*j*k*l*m \
        - sum(dE4(*g) for g in combinations([i, j, k, l, m], 4)) \
        - sum(dE3(*g) for g in combinations([i, j, k, l, m], 3)) \
        - sum(dE2(*g) for g in combinations([i, j, k, l, m], 2)) \
        - sum(dE1(*g) for g in combinations([i, j, k, l, m], 1))
def dE6gen(i, j, k, l, m, n):
    """Generate a single six-body term."""
    return i*j*k*l*m*n \
        - sum(dE5gen(*g) for g in combinations([i, j, k, l, m, n], 5)) \
        - sum(dE4gen(*g) for g in combinations([i, j, k, l, m, n], 4)) \
        - sum(dE3gen(*g) for g in combinations([i, j, k, l, m, n], 3)) \
        - sum(dE2gen(*g) for g in combinations([i, j, k, l, m, n], 2)) \
        - sum(dE1(*g) for g in combinations([i, j, k, l, m, n], 1))
def dE7gen(i, j, k, l, m, n, o):
    """Generate a single seven-body term."""
    return i*j*k*l*m*n*o \
        - sum(dE6gen(*g) for g in combinations([i, j, k, l, m, n, o], 6)) \
        - sum(dE5gen(*g) for g in combinations([i, j, k, l, m, n, o], 5)) \
        - sum(dE4gen(*g) for g in combinations([i, j, k, l, m, n, o], 4)) \
        - sum(dE3gen(*g) for g in combinations([i, j, k, l, m, n, o], 3)) \
        - sum(dE2gen(*g) for g in combinations([i, j, k, l, m, n, o], 2)) \
        - sum(dE1(*g) for g in combinations([i, j, k, l, m, n, o], 1))

def E1(f):
    """Generate all possible one-body terms."""
    return sum(dE1(*g) for g in combinations(f, 1))
def E2(f):
    """Generate all possible two-body terms."""
    return sum(dE2(*g) for g in combinations(f, 2))
def E3(f):
    """Generate all possible three-body terms."""
    return sum(dE3(*g) for g in combinations(f, 3))
def E4(f):
    """Generate all possible four-body terms."""
    return sum(dE4(*g) for g in combinations(f, 4))
def E2gen(f):
    """Generate all possible two-body terms."""
    return sum(dE2gen(*g) for g in combinations(f, 2))
def E3gen(f):
    """Generate all possible three-body terms."""
    return sum(dE3gen(*g) for g in combinations(f, 3))
def E4gen(f):
    """Generate all possible four-body terms."""
    return sum(dE4gen(*g) for g in combinations(f, 4))
def E5gen(f):
    """Generate all possible five-body terms."""
    return sum(dE5gen(*g) for g in combinations(f, 5))
def E6gen(f):
    """Generate all possible six-body terms."""
    return sum(dE6gen(*g) for g in combinations(f, 6))
def E7gen(f):
    """Generate all possible seven-body terms."""
    return sum(dE7gen(*g) for g in combinations(f, 7))

def T1(f):
    """Generate all possible one-body terms, using explicit loops."""
    n = len(f)
    terms = 0
    for i in range(n):
        terms += f[i]
    return terms


def T2(f):
    """Generate all possible two-body terms, using explicit loops."""
    n = len(f)
    terms = 0
    for i in range(n):
        for j in range(i+1, n):
            terms += (f[i]*f[j]) - f[i] - f[j]
    return terms


def T3(f):
    """Generate all possible three-body terms, using explicit loops."""
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


def MBE1_explicit(f):
    """Generate the full 1st-order many-body expansion expression, using
    functions containing explicit loops."""
    return T1(f)
def MBE2_explicit(f):
    """Generate the full 2nd-order many-body expansion expression, using
    functions containing explicit loops."""
    return MBE1_explicit(f) + T2(f)
def MBE3_explicit(f):
    """Generate the full 3rd-order many-body expansion expression, using
    functions containing explicit loops."""
    return MBE2_explicit(f) + T3(f)

def MBE1_itertools(f):
    return E1(f)
def MBE2_itertools(f):
    return MBE1_itertools(f) + E2gen(f)
def MBE3_itertools(f):
    return MBE2_itertools(f) + E3gen(f)
def MBE4_itertools(f):
    return MBE3_itertools(f) + E4gen(f)
def MBE5_itertools(f):
    return MBE4_itertools(f) + E5gen(f)
def MBE6_itertools(f):
    return MBE5_itertools(f) + E6gen(f)
def MBE7_itertools(f):
    return MBE6_itertools(f) + E7gen(f)
