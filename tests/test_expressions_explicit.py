import unittest

from sympy import symbols

import mbe.expressions_explicit as ee


class TestExpressionsExplicit(unittest.TestCase):

    def setUp(self):
        self.list_ten_bodies = symbols([
            'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
        ])

    def test_one_body_term_explicit(self):
        term_lhs = self.list_ten_bodies[0]
        term_rhs = ee.dE1(self.list_ten_bodies[0])
        self.assertEquals(term_lhs, term_rhs)

    def test_E1(self):
        pass

    def test_T1(self):
        pass

    def test_MBE1(self):
        expression_lhs = sum(self.list_ten_bodies)
        expression_rhs = ee.MBE1(self.list_ten_bodies)
        self.assertEquals(expression_lhs, expression_rhs)

if __name__ == "__main__":
    unittest.main()
