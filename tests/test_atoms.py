import pytest

from dsolve.atoms import Variable, Parameter, normalize_string

def test_normalize_string():
    assert normalize_string('+')=='+'
    
def test_variable():
    assert str(Variable('x_{t}'))=='x_{t}'
    assert str(Variable('x_t'))=='x_{t}'
    assert str(Variable('x_t^p'))=='x^{p}_{t}'
    assert str(Variable('E_{t}[x_{t+1}]'))=='E_{t}[x_{t+1}]'
    assert str(Variable('\pi^{p}_{t}'))=='\pi^{p}_{t}'
    assert str(Variable('x_{t}')(4))=='x_{4}'
    assert str(Variable('x_{t}').lag(4))=='x_{t-4}'
    assert str(Variable('x_{i,t}').subs({'i':0}))=='x_{0,t}'
    assert str(Variable('\theta_{t}'))==r'\theta_{t}'
    assert str(Variable('\theta_{t}'))=='\\theta_{t}'

def test_parameter():
    assert str(Parameter('\beta'))=='\\beta'
    assert str(Parameter('\beta'))==r'\beta'
    assert str(Parameter('\rho_{i}'))=='\\rho_{i}'
    assert str(Parameter('\sigma'))=='\sigma'
    assert str(Parameter('\theta'))=='\\theta'
    assert str(Parameter('\rho_{\theta}')) == '\\rho_{\\theta}'

