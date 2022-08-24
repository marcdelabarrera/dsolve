# dsolve

A package to solve systems of dynamic equations with Python. It understands
Latex syntax and it requires minimum specifications from the user end.


```python
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

system = Klein(eq = eq, x=x, p=p, z=z, calibration=calibration)
```

## Flexible input reading

The standarized way to write a variable is `E_{t}[x_{s}]` to represent the expectation of `x_{s}` at time `t`. but `dsolve` understands other formats. `Ex_{s}`, `E[x_s]` and `Ex_s` are quivalents to  `E_{t}[x_{s}]`, and the subscript `t` is assumed. 

Greek symbols can be writen as `\rho` or just `rho`. ´

`dsolve` understands fractions and sums. `\sum_{i=0}^{2}{x_{i,t}}` produces `x_{0,t}+x_{1,t}+x_{2,t}`