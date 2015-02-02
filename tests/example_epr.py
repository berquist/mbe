#!/usr/bin/env python

# pylint: disable=C0103,C0301

"""A practical example ...
"""

from __future__ import print_function

import argparse
import os

import mbe


parser = argparse.ArgumentParser()
parser.add_argument('--submit', action='store_true')
args = parser.parse_args()
submit = args.submit

# load the fragment input file
filename = os.path.join(os.environ['PWD'], 'example.xyz')

# generate the fragment (monomer) objects from the input
monomers = mbe.fragment.generate_fragment_objects(filename)

monomer_symbols = []
for m in monomers:
    for symbol in m.symbol_repr:
        monomer_symbols.append(symbol)

# generate our many-body expansion equations from the monomer symbols
two_body_expression = mbe.expressions_explicit.MBE2(monomer_symbols)
three_body_expression = mbe.expressions_explicit.MBE3(monomer_symbols)
four_body_expression = mbe.expressions.MBEn(monomer_symbols, 4)
five_body_expression = mbe.expressions.MBEn(monomer_symbols, 5)
six_body_expression = mbe.expressions.MBEn(monomer_symbols, 6)

print('Symbolic representation of monomers:\n')
print(monomer_symbols)
print('\nFull MBE2 expression:\n')
print(two_body_expression)
print('\nFull MBE3 expression:\n')
print(three_body_expression)
print('\nFull MBE4 expression:\n')
print(four_body_expression)
print('\nFull MBE5 expression:\n')
print(five_body_expression)
print('\nFull MBE6 expression:\n')
print(six_body_expression)
print('\n')

print('Removing the EPR silent fragments...\n')
two_body_expression = mbe.epr.remove_inactive_terms(two_body_expression, monomers)
three_body_expression = mbe.epr.remove_inactive_terms(three_body_expression, monomers)
four_body_expression = mbe.epr.remove_inactive_terms(four_body_expression, monomers)
five_body_expression = mbe.epr.remove_inactive_terms(five_body_expression, monomers)
six_body_expression = mbe.epr.remove_inactive_terms(six_body_expression, monomers)

print('EPR MBE2 expression:\n')
print(two_body_expression)
print('\nEPR MBE3 expression:\n')
print(three_body_expression)
print('\nEPR MBE4 expression:\n')
print(four_body_expression)
print('\nEPR MBE5 expression:\n')
print(five_body_expression)
print('\nEPR MBE6 expression:\n')
print(six_body_expression)
print('\n')

print('Forming fragments for MBE2...')
two_body_fragments = mbe.fragment.generate_fragments_from_expr(two_body_expression, monomers)
print('Forming fragments for MBE3...')
three_body_fragments = mbe.fragment.generate_fragments_from_expr(three_body_expression, monomers)
print('Forming fragments for MBE4...')
four_body_fragments = mbe.fragment.generate_fragments_from_expr(four_body_expression, monomers)
print('Forming fragments for MBE5...')
five_body_fragments = mbe.fragment.generate_fragments_from_expr(five_body_expression, monomers)
print('Forming fragments for MBE6...')
six_body_fragments = mbe.fragment.generate_fragments_from_expr(six_body_expression, monomers)

print('Generated all fragment combinations.')

if submit:
    results = []
    for fragment in two_body_fragments:
        # submit a PBS job for every fragment
        print('Submitting calculation for fragment {}:'.format(fragment.name))
        result = mbe.pbs.submit_fragment_job(fragment)
        results.append(result)

    # now that we have our results, do stuff with the g-tensor
    fragments_to_results = {k: v for k, v in zip(two_body_fragments, results)}
    symbols_to_fragments = mbe.utils.gen_dict_symbol2fragment(two_body_fragments)
