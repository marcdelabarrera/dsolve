import pytest

from dsolve.expressions import DynamicExpression

def test_expressions():
    assert str(DynamicExpression('E\pi_{t+1}+1')) == 'E_{t}[\pi_{t+1}]+1'
    assert str(DynamicExpression('E_t\pi_{t+1}+1')) == 'E_{t}[\pi_{t+1}]+1'
    assert str(DynamicExpression('E_t[\pi_{t+1}]+1')) == 'E_{t}[\pi_{t+1}]+1'

def test_fractions():
    assert str(DynamicExpression(r'\frac{a}{b}'))=='a/b'
    assert str(DynamicExpression(r'\frac{x+y}{b}'))=='(x+y)/b'
    assert str(DynamicExpression(r'\frac{a+b}{c+d}'))=='(a+b)/(c+d)'
    assert str(DynamicExpression(r'\frac{\frac{a}{b}}{c+d}'))=='a/(b*(c+d))'

def test_sums():
    assert str(DynamicExpression('\sum_{i=0}^{1}{x_{i}}'))=='x_{0}+x_{1}'
    assert str(DynamicExpression('\sum_{i=0}^{1}{x_{i,t}}'))=='x_{0,t}+x_{1,t}'

def test_special_characters():
    assert r'\rho'=='\\rho'
    assert str(DynamicExpression('\rho+1')) == '\\rho+1'