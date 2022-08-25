# dsolve

A package to solve systems of dynamic equations with Python. It understands $\LaTeX$ syntax and it requires minimum specifications from the user end. It solves problems of the form:

$$A_0\begin{bmatrix}x_{t+1}\\ E_{t}[p_{t+1}]\end{bmatrix}=A_1\begin{bmatrix}x_{t}\\ p_{t}\end{bmatrix}+\gamma z_t$$

Following Blanchard Kahn notation, $x_{t}$ are state variables (known at time $t$) while $p_{t}$ are forward-looking variables, and $z_t$ are shocks with $E_t[z_{t+1}]=0$.

Returns the matrix solution


$$p_t=\Theta_p x_t+Nz_t$$
$$x_{t+1}=\Theta_x x_t+Lz_t$$

and methods to plot impulse responses given a sequence of $z_t$

The main class of the package is `Klein`, which stores and solves the dynamic system. It takes a list of strings that are written as $\LaTeX$ equations, a dictionary that define the numeric values of the parameters, and the specification of `x`, `p` and `z`, specified as a list of $\LaTeX$ strings or a long string separated by commas.  


```python
from dsolve.solvers import Klein

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

# Simulate the inpulse response of a shock v_{0}=0 for 12 periods when \epsilon_{-1}=0

system.simulate(x0=0, z = {'v_{t}':1}, T=12)
```

## Flexible input reading

The standarized way to write a variable is `E_{t}[x_{s}]` to represent the expectation of `x_{s}` at time `t`. but `dsolve` understands other formats. `Ex_{s}`, `E[x_s]` and `Ex_s` are quivalents to  `E_{t}[x_{s}]`, and the subscript `t` is assumed. 

Greek symbols can be writen as `\rho` or just `rho`. 

`dsolve` understands fractions and sums. `\sum_{i=0}^{2}{x_{i,t}}` produces `x_{0,t}+x_{1,t}+x_{2,t}` and fraction `\frac{a}{b}` produces `(a)/(b)`