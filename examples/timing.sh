#!/usr/bin/env bash

# SETUP="
# from sympy import symbols
# from mbe.expressions_explicit import E1, E2, E3
# from mbe.expressions_explicit import T1, T2, T3
# f7 = list(symbols('a, b, c, d, e, f, g'))
# "
# echo ${SETUP}

# python -m timeit -r 5 -s "${SETUP}" "E1(f7)"
# python -m timeit -r 5 -s "${SETUP}" "T1(f7)"
# python -m timeit -r 5 -s "${SETUP}" "E2(f7)"
# python -m timeit -r 5 -s "${SETUP}" "T2(f7)"
# python -m timeit -r 5 -s "${SETUP}" "E3(f7)"
# python -m timeit -r 5 -s "${SETUP}" "T3(f7)"

SETUP="
from sympy import symbols
from mbe.expressions_explicit import dE2, dE3, dE4
from mbe.expressions_explicit import dE2gen, dE3gen, dE4gen
f7 = list(symbols('a, b, c, d, e, f, g'))
f2 = f7[:2]
f3 = f7[:3]
f4 = f7[:4]

f5 = f7[:5]
f6 = f7[:6]
"
echo "${SETUP}"

python -m timeit -r 5 -s "${SETUP}" "dE2(*f2)"
python -m timeit -r 5 -s "${SETUP}" "dE2gen(*f2)"
python -m timeit -r 5 -s "${SETUP}" "dE3(*f3)"
python -m timeit -r 5 -s "${SETUP}" "dE3gen(*f3)"
python -m timeit -r 5 -s "${SETUP}" "dE4(*f4)"
python -m timeit -r 5 -s "${SETUP}" "dE4gen(*f4)"

# python -m timeit -r 5 -s "${SETUP}" "dE5gen(*f5)"
# python -m timeit -r 5 -s "${SETUP}" "dE6gen(*f6)"

SETUP="
from sympy import symbols
from mbe.expressions_explicit import MBE1_explicit, MBE2_explicit, MBE3_explicit
from mbe.expressions import MBEn
f7 = list(symbols('a, b, c, d, e, f, g'))
f2 = f7[:2]
f3 = f7[:3]
f4 = f7[:4]

f5 = f7[:5]
f6 = f7[:6]
"
echo "${SETUP}"

echo "1st-order explicit vs. recursive"
python -m timeit -r 5 -s "${SETUP}" "MBE1_explicit(f7)"
python -m timeit -r 5 -s "${SETUP}" "MBEn(f7, 1)"

echo "2nd-order explicit vs. recursive"
python -m timeit -r 5 -s "${SETUP}" "MBE2_explicit(f7)"
python -m timeit -r 5 -s "${SETUP}" "MBEn(f7, 2)"

echo "3rd-order explicit vs. recursive"
python -m timeit -r 5 -s "${SETUP}" "MBE3_explicit(f7)"
python -m timeit -r 5 -s "${SETUP}" "MBEn(f7, 3)"

SETUP="
from sympy import symbols
from mbe.expressions_explicit import MBE1_explicit, MBE2_explicit, MBE3_explicit
from mbe.expressions import MBEn
alphabet = symbols('a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z')
"
echo "${SETUP}"

echo "1st-order explicit vs. recursive"
python -m timeit -r 5 -s "${SETUP}" "MBE1_explicit(alphabet)"
python -m timeit -r 5 -s "${SETUP}" "MBEn(alphabet, 1)"

echo "2nd-order explicit vs. recursive"
python -m timeit -r 5 -s "${SETUP}" "MBE2_explicit(alphabet)"
python -m timeit -r 5 -s "${SETUP}" "MBEn(alphabet, 2)"

echo "3rd-order explicit vs. recursive"
python -m timeit -r 5 -s "${SETUP}" "MBE3_explicit(alphabet)"
python -m timeit -r 5 -s "${SETUP}" "MBEn(alphabet, 3)"
