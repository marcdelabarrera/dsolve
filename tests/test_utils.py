from dsolve import utils
import pytest

def test_encode():
    assert utils.encode('\beta')=='\\beta'
    assert utils.encode('\alpha')=='\\alpha'
    assert utils.encode('\sigma')=='\sigma'

def test_normalize_string():
    assert utils.normalize_string('Ex_{t+1}')=='E_{t}[x_{t+1}]'
    assert utils.normalize_string('E[x_{t+1}]')=='E_{t}[x_{t+1}]'

def test_normalize_dict():
    assert utils.normalize_dict({'e_t':3})=={'e_{t}':3}