import pytest

from dsolve.expressions import DynamicExpression, DynamicEquation, classify_string
import sympy as sym 

def test_expressions():
    assert str(DynamicExpression('E\pi_{t+1}+1')) == 'E_{t}[\pi_{t+1}]+1'
    assert str(DynamicExpression('E_t\pi_{t+1}+1')) == 'E_{t}[\pi_{t+1}]+1'
    assert str(DynamicExpression('E_t[\pi_{t+1}]+1')) == 'E_{t}[\pi_{t+1}]+1'

def test_fractions():
    assert str(DynamicExpression(r'\frac{a}{b}'))=='a/b'
    assert str(DynamicExpression(r'\frac{x+y}{b}'))=='(x+y)/b'
    assert str(DynamicExpression(r'\frac{a+b}{c+d}'))=='(a+b)/(c+d)'
    assert str(DynamicExpression(r'\frac{\frac{a}{b}}{c+d}'))=='a/(b*(c+d))'

def test_scientific_notation():
    assert str(DynamicExpression('10e-4+x_{t}'))=='x_{t}+0.001'
    assert str(DynamicExpression('10e+4+x_{t}'))=='x_{t}+100000.0'
    
def test_sums():
    assert str(DynamicExpression('\sum_{i=0}^{1}{x_{i}}'))=='x_{0}+x_{1}'
    assert str(DynamicExpression('\sum_{i=0}^{1}{x_{i,t}}'))=='x_{0,t}+x_{1,t}'

def test_special_characters():
    assert r'\rho'=='\\rho'
    assert str(DynamicExpression('\rho+1')) == '\\rho+1'

def test_classify():
    assert classify_string('x_t')=='variable'
    assert classify_string('x_{i,t}')=='variable'
    assert classify_string('x_{\theta}')=='parameter'

def test_from_sympy():
    assert str(DynamicExpression(sym.Symbol('x')+1))=='x+1'
    assert str(DynamicEquation.from_sympy(sym.Eq(sym.Symbol('x'),1)))=='1 = x'

def test_subs():
    assert DynamicExpression('i_t-E_t[x_{t+1}]').subs({'i_t':3})
    assert DynamicExpression('i_t-E_t[x_{t+1}]').subs({'i_t':3, 'E_t[x_{t+1}]':2})
    assert float(DynamicExpression('i_t-E_t[x_{t+1}]').subs({'i_t':3, 'E_t[x_{t+1}]':2}))==1.

def test_indices():
    assert DynamicEquation('x_{it}=y_{it}').indices == ['i','t']