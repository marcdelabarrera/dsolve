from __future__ import annotations
from asyncore import read
from .expressions import DynamicEquation, close_brackets, DynamicExpression
from .atoms import Variable, E
from .utils import normalize_string, normalize_dict
from scipy.linalg import ordqz, inv
import matplotlib.pyplot as plt
import re
import numpy as np
import sympy as sym
from dataclasses import dataclass

@dataclass
class SystemVariables:
    '''
    Class that reads and contains all information regarding system variables.
    '''
    x : list[sym.Symbol]
    p : list[sym.Symbol]
    z : list[sym.Symbol]
    s : list[sym.Symbol]
    x1 : list[sym.Symbol]
    p1 : list[sym.Symbol]
 
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

class Klein():
    def __init__(self, 
                    equations:list[str]|str, 
                    x: list[str]|str = None, 
                    p: list[str]|str = None, 
                    z: list[str]|str = None,
                    s: list[str]|str = None, 
                    indices: dict[list[int]]= None,
                    calibration:dict = None):

        if indices is not None and len(indices)>1:
            raise ValueError('Systems with more than one index still not implemented')

        self.equations, self.vars, self.parameters = self.read_system(equations, x, p, z, s, indices)
        self.n_eq, self.n_x, self.n_p, self.n_z, self.n_s = len(self.equations.dynamic.symbolic), len(self.vars.x), len(self.vars.p), len(self.vars.z), len(self.vars.s)

        #if self.n_eq>(self.n_x+self.n_p):
        #    raise ValueError(f'More equations ({self.n_eq}) than unknowns ({self.n_x+self.n_p})')
        #elif self.n_eq<(self.n_x+self.n_p):
        #    raise ValueError(f'More unknowns ({self.n_x+self.n_p}) than equations ({self.n_eq}) ')
        self.type = self.classify_system()
        self.system_symbolic = self.get_matrices()
        self.system_numeric = None
        self.system_solution = None
        if calibration is not None:
            self.calibrate(calibration)
            self.system_solution = self.solve()
        else:
            self.calibration, self.system_numeric, self.solution = None, None, None

    def read_system(self,
                    equations:list[str]|str, 
                    x: list[str]|str = None, 
                    p: list[str]|str = None, 
                    z: list[str]|str = None,
                    s: list[str]|str = None, 
                    indices: dict[list[int]]= None):
        
        x,p,z,s = [self.split(i) if isinstance(i,str) else i for i in [x,p,z,s]]
        x,p,z,s = [[Variable(j) for j in i] if i is not None else [] for i in [x,p,z,s]]
        if indices is not None:
            x,p,z,s = [self.expand_indices(i, indices) for i in [x,p,z,s]]
        x1 = [ix.lead(1) for ix in x]
        p1 = [E(ip.lead(1),'t') for ip in p]
        x,p,z,s,x1,p1 = [[j.sympy for j in i] for i in [x,p,z,s,x1,p1]]
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
        free_symbols = set.union(*[eq.free_symbols for eq in dynamic_equations.symbolic+static_equations.symbolic])
        parameters = free_symbols.difference(x+p+z+s+x1+p1)
        return( SystemEquations(dynamic_equations, static_equations),\
                SystemVariables(x,p,z,s,x1,p1),\
                SystemParameters(parameters))
    
    @staticmethod
    def expand_indices(l:list=None, indices:dict[list[int]]=None)->list[Variable]:
        if indices is None or l is None:
            return l
        index = list(indices.keys())[0]
        start, end  = indices[index]
        out = []
        for el in l:
            if el.indexed:
                out = out + [el.subs({index:i}) for i in range(start, end+1)]
            else:
                out.append(el)
        return out
    
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
        A_0,_ = sym.linear_eq_to_matrix([i for i in self.equations.dynamic.symbolic], self.vars.x1+self.vars.p1)
        A_1,_ = sym.linear_eq_to_matrix(_,self.vars.x+self.vars.p)
        gamma, _ = sym.linear_eq_to_matrix(_, self.vars.z)
        return {'A_0':A_0, 'A_1':A_1, 'gamma':-gamma}

    def calibrate(self, calibration: dict[float], inplace=False)->dict[np.array]:
        '''
        Substitute numerical variables to 
        '''
        calibration = normalize_dict(calibration)
        self.equations.dynamic.calibrated =  [eq.subs(calibration) for eq in self.equations.dynamic.symbolic]
        self.equations.static.calibrated  =  [eq.subs(calibration) for eq in self.equations.static.symbolic]
        self.parameters.calibration = calibration
        self.system_numeric = {k: np.array(v.subs(calibration)).astype(np.float64) for k,v in self.system_symbolic.items()}


    def solve(self)->dict[np.ndarray]:
        '''
        Solves the system:

        p(t) = Theta_p*x(t)+Nz(t)
        x(t+1) = Theta_x*x(t)+Lz(t)
        '''
        
        if self.system_numeric is None:
            raise ValueError('System is not calibrated.')

        system_numeric = self.system_numeric
        A_0, A_1, gamma = system_numeric['A_0'], system_numeric['A_1'], system_numeric['gamma']
        S, T, _, _, Q, Z = ordqz(A_0, A_1, output='complex',sort=lambda alpha,beta: np.abs(beta/(alpha+1e-10))<1)
        Q = Q.conjugate().T
        n_s = len([i for i in np.abs(np.diag(T)/np.diag(S)) if i<1]) #number of stable eigenvalues
        
        #print(f'System with {n_s} stable eigenvalues and {self.n_x} pre-determined variables.')
        
        #if n_s>len(self.x):
        #    raise ValueError('Multiple solutions')

        #elif n_s<len(self.x):
        #    raise ValueError('No solution')

        if self.type == 'forward-looking': 
            return {'N': np.real(Z@(-inv(T)@Q@gamma))}

        elif self.type=='backward-looking':
            Theta_x = Z@inv(S)@T@inv(Z)
            L = Z@inv(S)@Q@gamma
            return {'Theta_x': np.real(Theta_x), 'L':np.real(L)}

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
            out = {str(k):np.zeros(T, dtype=float) for k in self.vars.z}
            for iz in z:
                t = Variable(iz).indices[0]
                out[f'{str(Variable(iz).base)}_{{t}}'][t]=z[str(Variable(iz))]
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

    def simulate(self, z:dict[np.array], x0: np.array=None, T:int=None):
        '''
        Simulates for a given path of shocks and initial conditions for predetermined variables.
        x0: initial conditions. Set to 0 if not specified.
        
        '''
        if x0 is None:
            x0 = np.zeros_like(self.vars.x)        
        z = self.normalize_z(z,T)
        
        if self.type=='mixed':
            return self.simulate_mixed_system(z, x0)
        
        if self.type=='backward-looking':
            return self.simulate_backward_looking_system(z, x0)
        
        if self.type=='forward-looking':
            return self.simulate_forward_looking_system(z)

    def simulate_mixed_system(self, z:np.array, x0: np.array)->dict[np.array]:
        sol = self.system_solution
        Theta_x, Theta_p, N, L = sol['Theta_x'], sol['Theta_p'], sol['N'], sol['L']            
        T = z.shape[1]
        x = np.zeros((self.n_x,T+1))
        x[:,0] = x0
        p=np.zeros((self.n_p,T+1))
        for t,iz in enumerate(z.T):
            iz = iz.reshape(self.n_z,-1)
            x[:,[t+1]] = Theta_x@x[:,[t]]+L@iz
            p[:,[t]] = Theta_p@x[:,[t]]+N@iz
        p[:,[t+1]] = Theta_p@x[:,[t+1]]+N@iz
        
        x1 = {str(ix1):ix1_t for ix1, ix1_t in zip(self.vars.x1, x[:,1:])}
        p1 = {str(ip1):ip1_t for ip1, ip1_t in zip(self.vars.p1, p[:,1:])}
        z = {str(iz):iz_t for iz, iz_t in zip(self.vars.z, z)}
        x = {str(ix):ix_t for ix, ix_t in zip(self.vars.x, x[:,:-1])}
        p = {str(ip):ip_t for ip, ip_t in zip(self.vars.p, p[:,:-1])}
        d =  z|x|p|x1|p1|{'t':np.array(range(T))}
        d = self.solve_static(d)
        return d

    def solve_static(self, d:dict[np.array])->dict[np.array]:
        if self.vars.s == []:
            return d
        for s in self.equations.static.calibrated:
            s_t = []
            for t in d['t']:
                d_t = {k:v[t] for k,v in d.items()}
                s_t.append(float(s.rhs.subs(d_t)))
            d[str(s.lhs)]= np.array(s_t)
        return d    

    def simulate_backward_looking_system(self, z:dict[np.array], x0: np.array):
        sol = self.system_solution
        Theta_x, L = sol['Theta_x'], sol['L']
        T = z.shape[1]
        x = np.zeros((self.n_x,T+1))
        x[:,0] = x0
        for t,iz in enumerate(z.T):
            iz = iz.reshape(self.n_z,-1)
            x[:,[t+1]] = Theta_x@x[:,[t]]+L@iz
        z = {str(iz):iz_t for iz, iz_t in zip(self.vars.z, z)}
        x = {str(ix):ix_t for ix, ix_t in zip(self.vars.x, x[:,:-1])}
        return z|x

    def simulate_forward_looking_system(self, z:dict[np.array]):
        raise ValueError('Purely forward looking systems are not implemented')

    def solve_static(self, d:dict[np.array])->dict[np.array]:
        if self.vars.s == []:
            return d
        for s in self.equations.static.calibrated:
            s_t = []
            for t in d['t']:
                d_t = {k:v[t] for k,v in d.items()}
                s_t.append(float(s.rhs.subs(d_t)))
            d[str(s.lhs)]= np.array(s_t)
        return d   

    def plot(self, ax, d, vars:str, param_dict:dict=None):
        vars=[str(Variable(i)) for i in vars.split(',')] 
        out=[]
        for ivar in vars:
            out.append(ax.plot(d['t'],d[ivar], label=rf'${ivar}$'))
        return out

    def plot_expr(self, ax, d, expr:str):
        expr = DynamicExpression(expr)
        y_t = []
        for t in d['t']:
            d_t = {k:v[t] for k,v in d.items()}
            y_t.append(float(expr.subs(d_t)))
        return ax.plot(d['t'],y_t)   
