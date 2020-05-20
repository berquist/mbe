import unittest

from sympy import symbols

import mbe.expressions as e
import mbe.expressions_explicit as ee


class TestExpressions(unittest.TestCase):

    def setUp(self):
        self.ten_symbols = symbols('a, b, c, d, e, f, g, h, i, j')

    def test_MBE1_equivalence(self):
        self.assertEqual(ee.MBE1_explicit(self.ten_symbols),
                         e.MBEn(self.ten_symbols, 1))

    def test_MBE2_equivalence(self):
        self.assertEqual(ee.MBE2_explicit(self.ten_symbols),
                         e.MBEn(self.ten_symbols, 2))

    def test_MBE3_equivalence(self):
        self.assertEqual(ee.MBE3_explicit(self.ten_symbols),
                         e.MBEn(self.ten_symbols, 3))

    def test_MBE4_equivalence(self):
        self.assertEqual(ee.MBE4_itertools(self.ten_symbols),
                         e.MBEn(self.ten_symbols, 4))

    def test_MBE5_equivalence(self):
        self.assertEqual(ee.MBE5_itertools(self.ten_symbols),
                         e.MBEn(self.ten_symbols, 5))

    @unittest.skip("takes too long to run")
    def test_MBE6_equivalence(self):
        self.assertEqual(ee.MBE6_itertools(self.ten_symbols),
                         e.MBEn(self.ten_symbols, 6))

    @unittest.skip("takes too long to run")
    def test_MBE7_equivalence(self):
        self.assertEqual(ee.MBE7_itertools(self.ten_symbols),
                         e.MBEn(self.ten_symbols, 7))


class TestExpressionsExplicitSingleTerms(unittest.TestCase):

    def setUp(self):
        self.ten_symbols = symbols('a, b, c, d, e, f, g, h, i, j')

    def test_single_one_body_term_explicit(self):
        self.assertEqual(self.ten_symbols[0],
                         ee.dE1(self.ten_symbols[0]))

    def test_single_two_body_term_equivalence(self):
        s = self.ten_symbols[:2]
        self.assertEqual(ee.dE2(*s),
                         ee.dE2gen(*s))

    def test_single_three_body_term_equivalence(self):
        s = self.ten_symbols[:3]
        self.assertEqual(ee.dE3(*s),
                         ee.dE3gen(*s))

    def test_single_four_body_term_equivalence(self):
        s = self.ten_symbols[:4]
        self.assertEqual(ee.dE4(*s),
                         ee.dE4gen(*s))


class TestExpressionsExplicitExpansions(unittest.TestCase):

    def setUp(self):
        self.ten_symbols = symbols('a, b, c, d, e, f, g, h, i, j')

    def test_one_body_terms_equivalence(self):
        self.assertEqual(ee.E1(self.ten_symbols),
                         ee.T1(self.ten_symbols))

    def test_two_body_terms_equivalence(self):
        self.assertEqual(ee.E2(self.ten_symbols),
                         ee.T2(self.ten_symbols))
        self.assertEqual(ee.E2(self.ten_symbols),
                         ee.E2gen(self.ten_symbols))

    def test_three_body_terms_equivalence(self):
        self.assertEqual(ee.E3(self.ten_symbols),
                         ee.T3(self.ten_symbols))
        self.assertEqual(ee.E3(self.ten_symbols),
                         ee.E3gen(self.ten_symbols))

    def test_MBE1_explicit(self):
        self.assertEqual(sum(self.ten_symbols),
                         ee.MBE1_explicit(self.ten_symbols))
        self.assertEqual(ee.MBE1_explicit(self.ten_symbols),
                         ee.MBE1_itertools(self.ten_symbols))

    def test_MBE2_explicit(self):
        self.assertEqual(ee.MBE2_explicit(self.ten_symbols),
                         ee.MBE2_itertools(self.ten_symbols))

    def test_MBE3_explicit(self):
        self.assertEqual(ee.MBE3_explicit(self.ten_symbols),
                         ee.MBE3_itertools(self.ten_symbols))


if __name__ == "__main__":
    unittest.main()
