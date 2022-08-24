import pytest

import numpy as np
from dsolve.solvers import Klein

def test_klein():
    eq = ['x_{t}=\rho*x_{t-1}+\sigma*eps_{t}']  
    calibration = {'\rho':0.8,'\sigma':1}      
    system = Klein(eq, x='x_{t-1}', z='eps_{t}', calibration=calibration)
    assert system.system_solution == {'Theta_x': np.array([[0.8]]), 'L': np.array([[1.]])}