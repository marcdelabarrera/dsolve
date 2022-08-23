import pytest

from dsolve.expressions import DynamicEcuation, DynamicExpression

def test_fractions():
    assert str(DynamicExpression(r'\frac{a}{b}'))=='a/b'
    assert str(DynamicExpression(r'\frac{x+y}{b}'))=='(x+y)/b'
    assert str(DynamicExpression(r'\frac{a+b}{c+d}'))=='(a+b)/(c+d)'
    assert str(DynamicExpression(r'\frac{\frac{a}{b}}{c+d}'))=='a/(b*(c+d))'

def test_sums():
    assert str(DynamicExpression('\sum_{i=0}^{1}{x_{i}}'))=='x_{0}+x_{1}'
    assert str(DynamicExpression('\sum_{i=0}^{1}{x_{i,t}}'))=='x_{0,t}+x_{1,t}'