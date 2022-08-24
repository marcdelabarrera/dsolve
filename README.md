# dsolve

A package to solve systems of dynamic equations with Python. It understands
Latex syntax and it requires minimum specifications from the user end.


```
from solvers import Klein

# Your latex equations here as a list of strings
eq=[
    '\pi_{t}=\beta*E\pi_{t+1}+\kappa*y_{t}+u_{t}',
    'y_{t}=Ey_{t+1}+(1-\phi)*E[\pi_{t+1}]+\epsilon_{t}',
    '\epsilon_{t} = \rho_v*\epsilon_{t-1}+v_{t}'
]

# Your calibration here as a dictionary
calibration = {'\beta':0.98,'\kappa':0.1,'\phi':1.1,'\rho_v':0.8}

# Define pre-determined variables, forward looking variables, and shocks as strings separated by commas or a list of strings.

x = '\epsilon_{t-1}'
p = '\pi_t, y_t'
z = 'v_t, u_t'

system = Klein(eq, x=x, p=p, z=z, calibration=calibration)
```
