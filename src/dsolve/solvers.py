from __future__ import annotations
from asyncore import read
import pandas as pd
from dsolve.expressions import DynamicEquation, close_brackets, DynamicExpression
from dsolve.atoms import Variable, E, Parameter
from dsolve.utils import normalize_string, normalize_dict
from scipy.linalg import ordqz, inv
import matplotlib.pyplot as plt
import re
import numpy as np
import sympy as sym
from dataclasses import dataclass, field
from itertools import product

@dataclass
class SystemVariables:
    '''
    Class that reads and contains all information regarding system variables.
    '''
    x : list[sym.Symbol] = field(default_factory=list)
    p : list[sym.Symbol] = field(default_factory=list)
    z : list[sym.Symbol] = field(default_factory=list)
    s : list[sym.Symbol] = field(default_factory=list)
    x1 : list[sym.Symbol] = field(default_factory=list)
    p1 : list[sym.Symbol] = field(default_factory=list)
 
@dataclass
class SystemEquations():
    dynamic: SystemEquationsDynamic
    static: SystemEquationsStatic = None

@dataclass
class SystemEquationsDynamic():
    symbolic: list[sym.Eq]
    calibrated: list[sym.Eq] = None

@dataclass
class SystemEquationsStatic():
    symbolic: list[sym.Eq]
    calibrated: list[sym.Eq] = None

@dataclass
class SystemParameters:
    parameters: list[sym.Symbol]
    calibration: dict[float] = None

    def __call__(self):
        return self.parameters

@dataclass
class SystemIndices:
    indices: list[str]
    start:list[int]
    end: list[int]

class DynamicSystem:
    @property
    def free_symbols(self):
        return set.union(*[eq.free_symbols for eq in self.equations.dynamic.symbolic+self.equations.static.symbolic])

    def __init__(self, 
                    equations:list[str]|str, 
                    x: list[str]|str = None, 
                    p: list[str]|str = None, 
                    z: list[str]|str = None,
                    s: list[str]|str = None, 
                    indices: dict[list[int]]= None,
                    calibration:dict = None):

        self.indexed = indices is not None
        self.indices, self.vars, self.equations, self.parameters = self.read_system(equations, x, p, z, s, indices)
        self.n_eq, self.n_x, self.n_p, self.n_z, self.n_s = len(self.equations.dynamic.symbolic), len(self.vars.x), len(self.vars.p), len(self.vars.z), len(self.vars.s)
        if self.n_eq>(self.n_x+self.n_p):
            raise ValueError(f'More equations ({self.n_eq}) than unknowns ({self.n_x+self.n_p})')
        elif self.n_eq<(self.n_x+self.n_p):
            raise ValueError(f'More unknowns ({self.n_x+self.n_p}) than equations ({self.n_eq}) ')
        self.type = self.classify_system()
        self.system_symbolic = self.get_matrices()
   
        self.system_numeric = None
        self.system_solution = None
        if calibration is not None:
            self.calibrate(calibration)
            self.system_solution = self.solve()
        else:
            self.calibration, self.system_numeric, self.solution = None, None, None
    
    def read_indices(self, indices:dict[list[int]])->SystemIndices:
        if indices is None:
            return None
        indices, start, end  = list(indices.keys()), [indices[i][0] for i in indices], [indices[i][1] for i in indices]
        indices = SystemIndices(indices, start, end)
        return indices

    def read_variables(self, x,p,z,s, indices:SystemIndices=None)->SystemVariables:
        x,p,z,s = [self.split(i) if isinstance(i,str) else i for i in [x,p,z,s]]
        x,p,z,s = [[Variable(j) for j in i] if i is not None else [] for i in [x,p,z,s]]
        if indices is not None:
            x,p,z,s = [self.expand_indices(i, indices) for i in [x,p,z,s]]
        x1 = [ix.lead(1) for ix in x]
        p1 = [E(ip.lead(1),'t') for ip in p]
        x,p,z,s,x1,p1 = [[j.sympy for j in i] for i in [x,p,z,s,x1,p1]]
        return SystemVariables(x,p,z,s,x1,p1)
    
    def read_equations(self, equations, s, indices):
        equations = [DynamicEquation(eq) for eq in equations]
        if indices is not None:
            equations = self.expand_indices(equations, indices)
        if s == []:
            dynamic_equations =SystemEquationsDynamic([eq.sympy for eq in equations])
            static_equations = SystemEquationsStatic([])
        else:
            static_equations = [eq for eq in equations if eq.lhs in s]
            d = {str(eq.lhs):eq.rhs for eq in static_equations}
            dynamic_equations = [eq for eq in equations if eq not in static_equations]
            for i, eq in enumerate(dynamic_equations):
                if len(eq.free_symbols.intersection(s))>0:
                    dynamic_equations[i]=eq.subs(d)
            static_equations =SystemEquationsStatic([eq.sympy for eq in static_equations])
            dynamic_equations =SystemEquationsDynamic([eq.sympy for eq in dynamic_equations])
        return SystemEquations(dynamic_equations, static_equations)

    def read_system(self,
                    equations:list[str]|str, 
                    x: list[str]|str = None, 
                    p: list[str]|str = None, 
                    z: list[str]|str = None,
                    s: list[str]|str = None, 
                    indices: dict[list[int]]= None):
        
        indices = self.read_indices(indices)
        vars = self.read_variables(x,p,z,s,indices)
        equations = self.read_equations(equations, vars.s, indices)
        free_symbols = set.union(*[eq.free_symbols for eq in equations.dynamic.symbolic+equations.static.symbolic])
        parameters = SystemParameters(list(free_symbols.difference(vars.x+vars.p+vars.z+vars.s+vars.x1+vars.p1)))
        return ( indices,vars,equations,parameters)
    
    def expand_index(self, l, i, start, stop):
        out = []
        for el in l:
            if sym.Symbol(i) in el.indices:
                out=out+[el.subs({i:j}) for j in range(start,stop+1)]
            else:
                out.append(el)
        return out

    def expand_indices(self, l:list=None, indices:SystemIndices=None)->list[Variable]:
        if indices is None or l is None:
            return l 
        out = l.copy()
        for i in range(len(indices.indices)):
            out = self.expand_index(out, indices.indices[i],indices.start[i],indices.end[i])
        return out
    
    def expand_calibration_index(self, calibration, index):
        index, start, stop = index
        for k,v in calibration.items(): 
            if Parameter(k).indices is not None and sym.Symbol(index) in Parameter(k).indices:
                calibration=calibration|{str(Parameter(k).subs({index:i})):v[i-start] for i in range(start, stop+1)}
                calibration.pop(k)
        return calibration

    def expand_calibration(self, calibration:dict)->dict:
        indices = self.indices
        for i in range(len(indices.indices)):
            calibration=self.expand_calibration_index(calibration, (indices.indices[i], indices.start[i],indices.end[i]))
        return calibration

    @staticmethod
    def split(string)->list[str]:
        l = re.split('(?<=,)|(?=,)',string)
        l = close_brackets(l)
        l = [i for i in l if i!='' and i!=',']
        return l

    def get_matrices(self)->dict[np.ndarray]:
        '''
        Given the system of equations, write it in the form: 
        A_0y(t+1) = A_1@y(t)+gamma@z(t)
        '''
        A,_ = sym.linear_eq_to_matrix([i for i in self.equations.dynamic.symbolic], self.vars.x1+self.vars.p1)
        B,_ = sym.linear_eq_to_matrix(_,self.vars.x+self.vars.p)
        gamma, C = sym.linear_eq_to_matrix(_, self.vars.z)
        return {'A':A, 'B':B, 'gamma':-gamma, 'C': C}

    def calibrate(self, calibration: dict[float], inplace=False)->dict[np.array]:
        '''
        Substitute numerical variables to 
        '''
        calibration = normalize_dict(calibration)
        if self.indexed:
            calibration = self.expand_calibration(calibration)
        self.equations.dynamic.calibrated =  [eq.subs(calibration) for eq in self.equations.dynamic.symbolic]
        self.equations.static.calibrated  =  [eq.subs(calibration) for eq in self.equations.static.symbolic]
        self.parameters.calibration = calibration

        try:
            self.system_numeric = {k: np.array(v.subs(calibration)).astype(np.float64) for k,v in self.system_symbolic.items()}
        except:
            print({str(i) for i in self.parameters()}.difference(calibration.keys()))
            raise ValueError('Error')

    def solve(self)->dict[np.ndarray]:
        '''
        Solves the system:

        p(t) = Theta_p*x(t)+Nz(t)
        x(t+1) = Theta_x*x(t)+Lz(t)
        '''
        
        if self.system_numeric is None:
            raise ValueError('System is not calibrated.')

        system_numeric = self.system_numeric
        A, B, gamma = system_numeric['A'], system_numeric['B'], system_numeric['gamma']
        S, T, _, _, Q, Z = ordqz(A, B, output='complex',sort=lambda alpha,beta: np.round(np.abs(beta/np.maximum(alpha,1e-15)),6)<=1)
        Q = Q.conjugate().T
        #n_s = len([i for i in np.abs(np.diag(T)/np.diag(S)) if i<1]) #number of stable eigenvalues
        n_s = len([_ for i in range(S.shape[0]) if np.abs(S[i,i])>1e-6 and np.round(np.abs(T[i,i]/S[i,i]),6)<=1])
        #print(f'System with {n_s} stable eigenvalues and {self.n_x} pre-determined variables.')
        
        if n_s>len(self.vars.x):
            print(f'Eigenvalues: {np.diag(np.abs(T/S))}')
            raise ValueError(f'Multiple solutions: {n_s} stable eigenvalues for {len(self.vars.x)} pre-determined variables')

        elif n_s<len(self.vars.x):
            print(f'Eigenvalues: {np.diag(np.abs(T/S))}')
            raise ValueError(f'No solution: {n_s} stable eigenvalues for {len(self.vars.x)} pre-determined variables')

        if self.type == 'forward-looking': 
            return {'N': np.real(Z@(-inv(T)@Q@gamma))}

        elif self.type=='backward-looking':
            Theta_x = Z@inv(S)@T@inv(Z)
            L = Z@inv(S)@Q@gamma
            return {'Theta_x': np.real(Theta_x), 'L':np.real(L), 'Theta_p':None, 'N':None}

        else:
            Theta_p = Z[n_s:,:n_s]@inv(Z[:n_s,:n_s])
            Theta_x = Z[:n_s,:n_s]@inv(S[:n_s,:n_s])@T[:n_s,:n_s]@inv(Z[:n_s,:n_s])
            M = -inv(T[n_s:,n_s:])@Q[n_s:,:]@gamma
            N = (Z[n_s:,n_s:]-Z[n_s:,:n_s]@inv(Z[:n_s,:n_s])@Z[:n_s,n_s:])@M
            L = Z[:n_s,:n_s]@inv(S[:n_s,:n_s])@((-T[:n_s,:n_s]@inv(Z[:n_s,:n_s])@Z[:n_s,n_s:]+T[:n_s,n_s:])@M+Q[:n_s,:]@gamma)
            return {'Theta_x':np.real(Theta_x),'Theta_p':np.real(Theta_p), 'N':np.real(N),'L':np.real(L)}

    def classify_system(self):
        if self.vars.x==[]:
            return 'forward-looking'
        elif self.vars.p==[]:
            return 'backward-looking'
        else:
            return 'mixed'

    def normalize_z(self, z:dict, T=None)->np.array:
        '''
        
        >>> self.normalize_z({'z_{0}':1.},T=4)
        np.array([1.,0.,0.,0.])

        >>> self.normalize_z({'z_{0}':1.,'z_2':2.},T=4)
        np.array([1.,0.,2.,0.])
        '''
        z = normalize_dict(z)
        if T is not None:
            out = {str(k):np.zeros(T+1, dtype=float) for k in self.vars.z}
            for iz in z:
                t = Variable(iz).indices[-1]
                out[str(Variable(iz).reset_t())][t]=z[str(Variable(iz))]
        else:
            T = np.max([len(v) for v in z.values()])
            out = {str(k):np.zeros(T,dtype=float) for k in self.vars.z}
            for k,v in z.items():
                out[k][:len(v)]=v
        return np.array([v for v in out.values()])
    
    def impulse_response(self, z:str):
        if self.type=='mixed':
            return self.simulate_mixed_system(z)
        
        if self.type=='backward-looking':
            return self.simulate_backward_looking_system(z)
        
        if self.type=='forward-looking':
            return self.simulate_forward_looking_system(z)

    def simulate_forward_looking_system(self, z:dict[np.array]):
        raise ValueError('Purely forward looking systems are not implemented')


def simulate_backward_looking_system(z, x0, Theta_x, L):
        T = z.shape[1]
        x = np.zeros((Theta_x.shape[0],T+1))
        x[:,0] = x0
        for t,iz in enumerate(z.T):
            iz = iz.reshape(-1, 1)
            x[:,[t+1]] = Theta_x@x[:,[t]]+L@iz
        return x

def simulate_mixed_system(z, x0, Theta_x, Theta_p, N, L):
     T = z.shape[1]
     x = np.zeros((Theta_x.shape[0],T+1))
     p = np.zeros((Theta_p.shape[0], T))
     x[:,0] = x0
     for t,iz in enumerate(z.T):
         iz = iz.reshape(-1, 1)
         x[:,[t+1]] = Theta_x@x[:,[t]]+L@iz
         p[:,[t]] = Theta_p@x[:,[t]]+N@iz
     return x, p
 
def add_static_variables(system:DynamicSystem, impulse_response:pd.DataFrame)->pd.DataFrame:
    if system.vars.s == []:
        return impulse_response
    for static_equation in system.equations.static.calibrated:
        static_variable = str(static_equation.lhs)
        impulse_response[static_variable] = 0
        for t in impulse_response['t']:
            impulse_response.loc[t, static_variable] = static_equation.rhs.subs(dict(impulse_response.loc[t]))
    return impulse_response 


def simulate(system:DynamicSystem, z, x0:np.array = None, T:int=None)->pd.DataFrame:
    '''
    Implementation of Klein 2000 method to solve dyamic systems of forward equations
    '''
    x0 = np.zeros_like(system.vars.x) if x0 is None else x0
    z = system.normalize_z(z,T)
    T = z.shape[1]-1 if T is None else T
    sol = system.system_solution
    Theta_x, Theta_p, N, L = sol['Theta_x'], sol['Theta_p'], sol['N'], sol['L']
    if system.type=='mixed':
        x, p = simulate_mixed_system(z, x0, Theta_x, Theta_p, N, L)
        impulse_response = pd.DataFrame(np.row_stack((np.arange(T+1), x[:,:-1], x[:,1:], p, z)).T,
                           columns = ['t']+system.vars.x+system.vars.x1+system.vars.p+system.vars.z)

    if system.type == 'backward-looking':
        x = simulate_backward_looking_system(z, x0, Theta_x, L)
        impulse_response = pd.DataFrame(np.row_stack((np.arange(T+1), x[:,:-1], x[:,1:], z)).T,
                           columns = ['t']+system.vars.x+system.vars.x1+system.vars.z)
    impulse_response = impulse_response.astype({'t':int})
    impulse_response = add_static_variables(system, impulse_response)
    return impulse_response




class MITShock:
    def __init__(self, d:dict, model:Klein):
        self.d = d
        self.model = model
    
    def __call__(self, var:str):
        var = str(Variable(var))
        return self.d[var]
        
    def plot(self, ax, vars:str, **kwargs):
        
        vars=[str(Variable(i)) for i in vars.split(',')]
        out=[]
        for ivar in vars:
            label = kwargs['label'] if 'label' in kwargs else rf'${ivar}$' 
            out.append(ax.plot(self.d['t'],self.d[ivar], label=label))
        if 'title' in kwargs:
            out.append(ax.set(title=kwargs['title']))
        if 'legend' in kwargs and kwargs['legend']: 
            out.append(ax.legend())
        
        return out

    def plot_expr(self, ax, eq:str, **kwargs):
        eq = DynamicEquation(eq)
        eq = eq.subs(self.model.parameters.calibration)
        label = kwargs['label'] if 'label' in kwargs else None
        y_t = []
        out = []
        for t in self.d['t']:
            d_t = {k:v[t] for k,v in self.d.items()}
            y_t.append(float(eq.subs(d_t).rhs))
        out.append(ax.plot(self.d['t'], y_t, label=label))
        if 'title' in kwargs:
            out.append(ax.set(title=kwargs['title']))
        if 'legend' in kwargs and kwargs['legend']:
            out.append(ax.legend())
        
        return out
