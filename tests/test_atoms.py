import pytest
from dsolve.atoms import Variable, Parameter, E
from sympy import Symbol, sympify

def test_variable():
    assert str(Variable('x_{t}'))=='x_{t}'
    assert str(Variable('x_t'))=='x_{t}'
    assert str(Variable('x_t^p'))=='x^{p}_{t}'
    assert str(Variable('E_{t}[x_{t+1}]'))=='E_{t}[x_{t+1}]'
    assert str(Variable('\pi^{p}_{t}'))=='\pi^{p}_{t}'
    assert str(Variable('x_{t}')(4))=='x_{4}'
    assert str(Variable('x_{t}').lag(4))=='x_{t-4}'
    assert Variable('x_{i,t}').eval(2).value==2.
    assert str(Variable('x_{i,t}').subs({'i':0}))=='x_{0,t}'
    assert str(Variable('x_{i,t+1}').subs({'t':0}))=='x_{i,1}'
    assert str(Variable('\theta_{t}'))==r'\theta_{t}'
    assert str(Variable('\theta_{t}'))=='\\theta_{t}'
    assert Variable('x_{i,t}').indexed
    assert Variable('x_{it}').indexed
    assert Variable('x_{it+1}').indexed




def test_split():
    assert Variable.split('x_{t}')==(None, Symbol('x'),[Symbol('t')])
    assert Variable.split('E_{t}[x_{t+1}]')==(sympify('t'), Symbol('x'),[sympify('t+1')])

def test_Parameter():
    assert str(Parameter('\beta'))=='\\beta'
    assert str(Parameter('\beta'))==r'\beta'
    assert Parameter('\rho_{v}').indexed==False
    assert str(Parameter('\rho_{i}'))=='\\rho_{i}'
    assert str(Parameter('\sigma'))=='\sigma'
    assert str(Parameter('\theta'))=='\\theta'
    assert str(Parameter('\rho_{\theta}')) == '\\rho_{\\theta}'
    assert str(Parameter('\rho_{i}').subs({'i':4})) == '\\rho_{4}'
    assert str(Parameter('\rho').subs({'i':4})) == '\\rho'

def test_E():
    assert str(E(Variable('x_{t+1}'),'t'))=='E_{t}[x_{t+1}]'
    assert str(E(Variable('x_{t+1}')))=='E_{t}[x_{t+1}]'
    assert str(E(Variable('x_{t+1}'),0))=='E_{0}[x_{t+1}]'