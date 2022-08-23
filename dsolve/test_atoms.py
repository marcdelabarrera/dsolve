import pytest

from dsolve.atoms import Variable

def test_variable():
    assert str(Variable('x_{t}'))=='x_{t}'