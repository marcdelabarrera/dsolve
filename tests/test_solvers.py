import pytest

import numpy as np
from dsolve.solvers import Klein, SystemVariables, SystemEquations
import sympy as sym

def test_SystemVariables():
    assert str(SystemVariables(x=[sym.Symbol('x_{t}')]).x[0])=='x_{t}'
    #assert str(SystemVariables(x='x_{i,t}', indices={'i':(0,1)}).x[1])=='x_{1,t}'

#def test_SystemEquations():
#    assert SystemEquations(equations=['x_{t}+y_{t}=z_{t}']).static_equations==[]
#    s = SystemEquations(equations=['x_t+y_t=z_t', 'i_t=3*x_t'], vars=SystemVariables(s='i_t'))
#    assert s.static_equations==[sym.Eq(sym.Symbol('i_{t}'), 3*sym.Symbol('x_{t}'))]
#    s = SystemEquations(equations=['x_{i,t}=1'], indices={'i':(0,3)})
#    assert s.dynamic_equations[2]==sym.Eq(sym.Symbol('x_{2,t}'),1)




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
    assert np.all(system.normalize_z({'eps_{0}':1},T=4)==np.array([1.,0,0,0]))
    assert np.all(system.normalize_z({'eps_{0}':1, 'eps_{2}':2.},T=4)==np.array([1.,0,2.,0]))

def test_indices():
    eq = ['x_{it}=rho*x_{it-1}+eps_{t}']
    #system = Klein(eq, indices={'i':(0,1)})
    #assert str(system.equations.dynamic.symbolic[0]) == 'Eq(x_{0,t}, eps_{t} + rho*x_{0,t-1})'
    eq = ['x_{ijt}=rho_i*x_{ijt-1}+eps_{t}']
    system = Klein(eq, x = 'x_{ijt-1}', z='eps_{t}', indices={'i':(0,1), 'j':(0,1)})
    assert str(system.equations.dynamic.symbolic[0]) == 'Eq(x_{0,0,t}, eps_{t} + rho_{0}*x_{0,0,t-1})'
