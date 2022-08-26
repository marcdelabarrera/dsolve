import pytest

import numpy as np
from dsolve.solvers import Klein, SystemVariables


def test_SystemVariables():
    assert str(SystemVariables(x='x_t').x[0])=='x_{t}'
    assert str(SystemVariables(x='x_{i,t}', indices={'i':(0,1)}).x[1])=='x_{1,t}'


def test_Klein():
    #simple AR(1)
    eq = ['x_{t}=\rho*x_{t-1}+\sigma*eps_{t}']  
    calibration = {'\rho':0.8,'\sigma':1}      
    system = Klein(eq, x='x_{t-1}', z='eps_{t}', calibration=calibration)
    assert system.system_solution == {'Theta_x': np.array([[0.8]]), 'L': np.array([[1.]])}

    # NK Model with static equation
    eq=[
    'x_{t}=b*Ex_{t+1}+kappa*y_{t}',
    'y_{t}=Ey_{t+1}-\frac{1}{\sigma}*(i_t-x_{t}-rho)',
    'i_t=rho+phi*x_{t}+v_t',
    'v_t = rho_v*v_{t-1}+eps_t'
    ]
    calibration = {'b':0.98,'kappa':0.1,'phi':1.1,'rho_v':0.8, 'rho':0.02,'\sigma':1}
    system = Klein(eq, x='v_{t-1}', p='x_t,y_t', z='eps_t', s='i_t', calibration=calibration)
    assert isinstance(system.system_solution, dict)
    assert system.n_s == 1 