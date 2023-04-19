from jax import Array
import jax
import jax.numpy as jnp
from scipy.optimize import root
from typing import Callable
from dataclasses import dataclass
from functools import partial

equilibrium_conditions = Callable[[Array, Array, Array, Array, Array, Array],Array]

@partial(jax.jit, static_argnames='h')
def lead(X, h=1, fill_value=None):
    fill_value = X[-1:] if fill_value is None else fill_value
    if len(X.shape)==1:
        return jnp.concatenate([X[h:], jnp.tile(fill_value, h)])
    if len(X.shape)==2:
        return jnp.concatenate([X[h:], jnp.tile(fill_value, (h,1))])
    
@partial(jax.jit, static_argnames='h')
def lag(X, h=1, fill_value=None):
    fill_value = X[0] if fill_value is None else fill_value
    if len(X.shape)==1:
        return jnp.concatenate([jnp.tile(fill_value, h), X[:-h], ])
    if len(X.shape)==2:
        return jnp.concatenate([jnp.tile(fill_value, (h,1)), X[:-h]])

def _build_F(f: equilibrium_conditions, initial_ss: Array, final_ss:Array= None, jit:bool=True)->Array:
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
    final_ss = initial_ss if final_ss is None else final_ss
    def F(X, Eps):
        X_, X1 = lag(X, fill_value = initial_ss), lead(X, fill_value = final_ss)
        Eps_, Eps1 = lag(Eps,fill_value=jnp.zeros(Eps.shape[1])), lead(Eps)
        return jax.vmap(f)(X_,X,X1,Eps_, Eps, Eps1)
    
    F = jax.jit(F) if jit else F
    return F

def solve_impulse_response(f: equilibrium_conditions, Eps:Array, initial_ss: Array, final_ss:Array= None):
    '''
    Solves the equilibrium for a given path of shocks Eps.
    Parameters
    ----------
    f: callable
        Residual of the equilibrium conditions f(x_{t-1},x_{t},x_{t+1},eps_{t-1},eps_{t},eps_{t+1})=0
    Eps: Array
        Path of shocks for which to solve the impulse response. It has shape (T x n_shocks)
    initial_ss: Array
        Initial steady state that defines x_{-1}
    final_ss: Array
        Final steady state that defines x_{T+1}. If None, assume the model reverts back to the initial steady state
    '''
    F = _build_F(f, initial_ss, final_ss)
    T = Eps.shape[0]
    n_x = len(initial_ss)
    X_guess = jnp.tile(initial_ss,(T,1))
    sol = root(lambda x: F(x.reshape(-1,n_x), Eps).flatten(), x0=X_guess.flatten())
    X = sol.x.reshape(-1,n_x)
    return X


# @dataclass
# class SequenceSpaceModel:
#     '''
#     Dataclass that defines a model to be solved using the sequence space method, which involves solving
#     a system of equations F(X,Eps)=0.
#     '''
#     f: equilibrium_conditions
#     T: int
#     ss0: Array #TODO: ss0 could be computed from equilibrium_conditions
#     ssT: Array = None
#     jit: bool = True

#     def __post_init__(self):
#         self.ssT = self.ss0 if self.ssT is None else self.ssT
#         self.F = _build_F(self.f, self.T, self.ss0, self.ssT, jit=self.jit)
#         self.n_x = len(self.ss0)


# def solve_impulse_response(model:SequenceSpaceModel, Eps: Array, X_guess: Array = None)->Array:
#     '''
#     Solves a SequenceSpaceModel for a given path of perfect-foresight shocks Eps
#     Parameters
#     ----------
#     model: SequenceSpaceModel
#     Eps: Array
#     Returns
#     -------
#     X: Array
#     '''
#     if len(Eps)!=model.T+1:
#         raise ValueError(f'Wrong dimensions for Eps, len should be {model.T+1}')
    
#     n_x = model.n_x
#     X_guess = X_guess if X_guess is not None else jnp.tile(model.ss0,(model.T+1,1))

#     sol = root(lambda x: model.F(x.reshape(-1,n_x),Eps).flatten(), x0=X_guess.flatten())

#     X = sol.x.reshape(-1,n_x)

#     #if jnp.max(jnp.abs(sol.fun))>1e-5:
#     #    print(sol)
#     #    raise ValueError(f'Solution not achieved')
#     return X


# # def lead(x: Array, xT: float=None, h:int=1)->Array:
#     '''
#     Given vector x_{t}, returns x_{t+h}
#     Example
#     -------
#     x = jnp.array([1,2,3])
#     lead(x)
#     [2,3,3]
#     '''
#     #TODO: implement h
#     xT = xT if xT is not None else x[-1]
#     return jnp.concatenate([x[1:], jnp.atleast_1d(xT)])

# def lag(x: Array, xT: float=None, h:int=1)->Array:
#     '''
#     Given vector x_{t}, returns x_{t-h}
#     Example
#     -------
#     x = jnp.array([1,2,3])
#     lag(x)
#     [1,1,2]
#     '''
#     xT = xT if xT is not None else x[-1]
#     return jnp.concatenate([jnp.atleast_1d(xT), x[1:]])