#!/usr/bin/env python

import os

import mbe

# load the fragment input file
filename = os.path.join(os.environ['PWD'], 'example.xyz')

# generate the fragment (monomer) objects from the input
monomers = mbe.fragment.generate_fragment_objects(filename)

# create a map to go from a fragment's symbol representation to the full object
frag_map = mbe.utils.gen_dict_symbol2fragment(monomers)

monomer_symbols = []
for m in monomers:
    for symbol in m.symbol_repr:
        monomer_symbols.append(symbol)

two_body_expression = mbe.expressions_explicit.MBE2(monomer_symbols)
three_body_expression = mbe.expressions_explicit.MBE3(monomer_symbols)

print('Symbolic representation of monomers:')
print(monomer_symbols)
print('\n')
print('Full MBE2 expression:')
print(two_body_expression)
print('Full MBE3 expression:')
print(three_body_expression)
print('\n')
print('Forming dimers...')
# print('Forming dimers and trimers...')

# The term dict maps the symbolic representation of fragments (with their
# union shown as multiplication) to their coefficient in the final expression.
term_dict = dict(two_body_expression.as_coefficients_dict())
# term_dict = dict(three_body_expression.as_coefficients_dict())
fragments = []
for term in term_dict:
    # A term represents a fragment or a union of fragments.
    fragment_monomer_symbols = term.atoms()
    if len(fragment_monomer_symbols) > 1:
        fragment_monomers = []
        for symbol in fragment_monomer_symbols:
            fragment_monomers.append(frag_map[symbol])
        fragment = mbe.fragment.combine_fragment_sequence(fragment_monomers)
    else:
        fragment = frag_map[next(iter(fragment_monomer_symbols))]
    fragments.append(fragment)

print('Generated all fragment combinations.')

results = []
for fragment in fragments:
    # submit a PBS job for every fragment
    print('Submitting calculation for fragment {}:'.format(fragment.name))
    result = mbe.pbs.submit_fragment_job(fragment)
    results.append(result)
