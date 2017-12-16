import sympy

t1 = (1, {'a'})
t2 = (1, {'b'})
t3 = (1, {'c'})
t4 = (1, {'d'})

t12 = (1, {'a', 'b'})
t13 = (1, {'a', 'c'})
t14 = (1, {'a', 'd'})
t23 = (1, {'b', 'c'})
t24 = (1, {'b', 'd'})
t34 = (1, {'c', 'd'})

t123 = (1, {'a', 'b', 'c'})
t124 = (1, {'a', 'b', 'd'})
t134 = (1, {'a', 'c', 'd'})
t234 = (1, {'b', 'c', 'd'})

t1234 = (1, {'a', 'b', 'c', 'd'})

class Expression:
    """
    An expression is a set of terms.

    >>> Expression( {(1, {'a'}), (1, {'a'})} )
    Expression( {(2, {'a'}) }
    >>> Expression( {(1, {'a'}), (1, {'a', 'b'})} )
    Expression( {(1, {'a'}), (1, {'a', 'b'})} )
    >>> Expression( {(1, {'a'}), (2, {'a'})} )
    Expression( {(3, {'a'})} )
    """
    def __init__(self, expr = set()):
        self.terms = expr

    def __add__(self, other):
        """
        Add two expressions together. If they have matching terms,
        simplify them.
        """
        new = Expression()
        return new

class Term:
    """
    A term is a combination of one or more fragments with an integer weight
    that might appear in a many-body expansion Expression.
    """
    def __init__(self):
        self.weight = 0
        self.fragments = set()

    def __add__(self, other):
        """
        We only add two terms together if their respective fragments are
        identical.
        """
        if self.fragments == other.fragments:
            new = Term()
            new.weight = self.weight + other.weight
            new.fragments = self.fragments
            return new
        pass

    def __sub__(self, other):
        """
        We only subtract two terms from each other if their respective
        fragments are identical.
        """
        if self.fragments == other.fragments:
            new = Term()
            new.weight = self.weight - other.weight
            new.fragments = self.fragments
            return new
        pass

    def __iadd__(self, other):
        pass

    def __isub__(self, other):
        pass

def create_map(sympy_expression, list_of_strings):
    """
    Return a dictionary 
