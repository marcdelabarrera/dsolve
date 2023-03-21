from jax import Array
import jax
import jax.numpy as jnp
from scipy.optimize import root

def build_F(f: callable, T:int, ss0: Array, ssT= None)->Array:
    '''
    Builds F, which stacks f(x_{t-1},x_{t},x_{t+1},eps_{t}) T times with initial condition x_{-1}=ss0 and
    x_{T+1}=ssT. If ssT is not given, ssT=ss0
    
    Parameters
    ----------
    f: callable
        Function to be stacked. The signature is (x_, x, x1, eps)
    T: int
        Number of time periods to consider
    ss0: Array
        Array with the steady state value
    ssT: Array
        Steady state at the terminal condition. If None, ss0
    
    Returns
    -------
    F: callable
        Equilibrium conditions stacked
    '''
    ssT = ss0 if ssT is None else ssT
    @jax.jit
    def F(X, Eps):
        if len(X)!=T+1 or len(Eps)!=T+1:
            raise ValueError('Incorrect shapes')
        out = jnp.zeros((T+1, len(ss0)))
        out = out.at[0].set(f(x_ = ss0, x = X[0,:], x1 = X[1,:], eps = Eps[0,:]))
        for t in range(1,T-1):
            out = out.at[t].set(f(x_ = X[t-1,:], x= X[t,:], x1 = X[t+1,:], eps = Eps[t,:]))
        out = out.at[T].set(f(x_ = X[T-1,:], x = X[T,:], x1 = ssT, eps = Eps[t,:]))
        return out
    F.T = T
    F.ss0 = ss0
    return F


def solve_model(F:callable, Eps:Array)->Array:
    '''
    Solves the system F(X,E)=0 using scipy root function. 
    Parameters
    ----------
    F: callable
        Callable that takes X, E and returns the error of the dynamic conditions
    Eps: Array
        Vector of shocks
    Returns
    -------
    X: Array
        Returns the solution of the dynamic system.
    '''
    n_x = len(F.ss0)
    X_guess = jnp.tile(F.ss0,(F.T+1,1))
    sol = root(lambda x: F(x.reshape(-1,n_x),Eps).flatten(), x0=X_guess.flatten())
    X = sol.x.reshape(-1,n_x)
    return X