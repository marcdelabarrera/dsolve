from jax import Array
import jax
import jax.numpy as jnp
from scipy.optimize import root
from typing import Callable
from dataclasses import dataclass


equilibrium_conditions = Callable[[Array, Array, Array, Array, Array, Array],Array]

@dataclass
class SequenceSpaceModel:
    '''
    Dataclass that defines a model to be solved using the sequence space method, which involves solving
    a system of equations F(X,Eps)=0.
    '''
    f: equilibrium_conditions
    T: int
    ss0: Array
    ssT: Array = None

    def __post_init__(self):
        self.ssT = self.ss0 if self.ssT is None else self.ssT
        self.F = _build_F(self.f, self.T, self.ss0, self.ssT, jit=True)
        self.n_x = len(self.ss0)


def _build_F(f: equilibrium_conditions, T:int, ss0: Array, ssT:Array= None, jit:bool=True)->Array:
    '''
    Builds F, which stacks f(x_{t-1},x_{t},x_{t+1},eps_{t-1},eps_{t},eps_{t+1}) T times with initial condition x_{-1}=ss0 and
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
        Equilibrium conditions stacked. Shape is (T,n_x).
    '''
    f = lambda **x: f(**x)
    def F(X, Eps):
        # check correct array size for X and Eps
        if len(X)!=T+1 or len(Eps)!=T+1: 
            raise ValueError(f'Incorrect shapes')
        
        out = jnp.zeros((T+1, len(ss0)))
        out = out.at[0].set(f(x_ = ss0, x = X[0,:], x1 = X[1,:], eps_ = Eps[0,:]*0, eps = Eps[0,:], eps1 = Eps[1,:]))
        for t in range(1,T):
            out = out.at[t].set(f(x_ = X[t-1,:], x= X[t,:], x1 = X[t+1,:], eps_ = Eps[t-1,:], eps = Eps[t,:], eps1 = Eps[t+1,:]))
        out = out.at[T].set(f(x_ = X[T-1,:], x = X[T,:], x1 = ssT, eps_ = Eps[T-1,:], eps = Eps[T,:], eps1 = 0*Eps[T,:]))
        return out
    
    F = jax.jit(F) if jit else F
    return F


def solve_impulse_response(model:SequenceSpaceModel, Eps: Array, X_guess: Array = None)->Array:
    '''
    Solves a SequenceSpaceModel for a given path of perfect-foresight shocks Eps
    Parameters
    ----------
    model: SequenceSpaceModel
    Eps: Array
    Returns
    -------
    X: Array
    '''
    if len(Eps)!=model.T+1:
        raise ValueError(f'Wrong dimensions for Eps, len should be {model.T+1}')
    
    n_x = model.n_x
    X_guess = X_guess if X_guess is not None else jnp.tile(model.ss0,(model.T+1,1))

    sol = root(lambda x: model.F(x.reshape(-1,n_x),Eps).flatten(), x0=X_guess.flatten())

    X = sol.x.reshape(-1,n_x)

    if jnp.max(jnp.abs(sol.fun))>1e-5:
        print(sol)
        raise ValueError(f'Solution not achieved')
    return X

