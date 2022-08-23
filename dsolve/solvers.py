from expressions import DynamicEcuation, close_brackets
from atoms import Variable, E
from scipy.linalg import ordqz, inv
import matplotlib.pyplot as plt
import re
import numpy as np
import sympy as sym

class Klein():



    def __init__(self, 
                    equations:list[str]|str, 
                    x: list[str]|str = None, 
                    p: list[str]|str = None, 
                    z: list[str]|str = None, 
                    indices: dict[list[int]]= None,
                    calibration:dict = None):

        if indices is not None and len(indices)>1:
            raise ValueError('Systems with more than one index still not implemented')

        self.indices = indices
        self.indexed = indices is not None
        self.equations = self.read_equations(equations)
        self.x, self.x1, self.n_x = self.read_x(x)
        self.p, self.p1, self.n_p = self.read_p(p)
        self.z, self.n_z = self.read_z(z)
        self.type = self.classify_system()
        self.parameters = {k:v for d in [e.parameters for e in self.equations] for k,v in d.items()}
        self.system_symbolic = self.get_matrices()
        self.calibration = calibration
        self.system_numeric = None
        self.system_solution = None
        if calibration is not None:
             self.system_numeric = self.calibrate(calibration)
             self.system_solution = self.solve()
        
    def get_matrices(self)->dict[np.ndarray]:
        '''
        Given the system of equations, write it in the form: 
        A_0y(t+1) = A_1@y(t)+gamma@z(t)
        '''
        A_0,_ = sym.linear_eq_to_matrix([i.sympy for i in self.equations], self.x1+self.p1)
        A_1,_ = sym.linear_eq_to_matrix(_,self.x+self.p)
        gamma, _ = sym.linear_eq_to_matrix(_, self.z)
        return {'A_0':A_0, 'A_1':A_1, 'gamma':-gamma}

    @staticmethod
    def split(string)->list[str]:
        l = re.split('(?<=,)|(?=,)',string)
        l = close_brackets(l)
        l = [i for i in l if i!='' and i!=',']
        return l

    def calibrate(self, calibration: dict[float], inplace=False)->dict[np.array]:
        '''
        Substitute numerical variables to 
        '''
        #calibration={str(Parameter(k)):v for k,v in calibration.items()}
        self.calibration = calibration
        if inplace:
            self.system_numeric = {k: np.array(v.subs(calibration)).astype(np.float64) for k,v in self.system_symbolic.items()}
        else:
            return {k: np.array(v.subs(calibration)).astype(np.float64) for k,v in self.system_symbolic.items()}

    def solve(self)->dict[np.ndarray]:
        '''
        Solves the system:

        p(t)=Theta_p*x(t)+Nz(t)
        x(t+1)=Theta_x*x(t)+Lz(t)
        '''
        
        if self.calibration is None:
            raise ValueError('Calbrate the system first')

        system_numeric = self.system_numeric
        A_0, A_1, gamma = system_numeric['A_0'], system_numeric['A_1'], system_numeric['gamma']
        S, T, _, _, Q, Z = ordqz(A_0, A_1, output='complex',sort=lambda alpha,beta: np.abs(beta/(alpha+1e-10))<1)
        Q = Q.conjugate().T
        n_s = len([i for i in np.abs(np.diag(T)/np.diag(S)) if i<1]) #number of stable eigenvalues
        
        print(f'System with {n_s} stable eigenvalues and {self.n_x} pre-determined variables.')
        
        if n_s>len(self.x):
            raise ValueError('Multiple solutions')

        elif n_s<len(self.x):
            raise ValueError('No solution')

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

    def expand_indices(self, l:list)->list:
        index = list(self.indices.keys())[0]
        start, end  = self.indices[index]
        out = []
        for el in l:
            if el.indexed:
                out = out + [el.subs({index:i}) for i in range(start, end+1)]
            else:
                out.append(el)
        return out

    def read_equations(self, equations:list[str])->list[DynamicEcuation]:
        equations = [DynamicEcuation(eq) for eq in equations]
        if self.indexed:
            equations = self.expand_indices(equations)
        return [eq for eq in equations]
            
    def read_x(self, x:list[str]|str)->list[sym.Symbol]:
        if x is None:
            return [],[],0
        if isinstance(x,str):
            x = self.split(x)
        x = [Variable(ix) for ix in x]
        x1 = [ix.lead(1) for ix in x]
        if self.indexed:
            x = self.expand_indices(x)
            x1 = self.expand_indices(x1)
        x = [ix.sympy for ix in x]
        x1 = [ix1.sympy for ix1 in x1]
        n_x = len(x)
        return (x,x1,n_x)

    def read_p(self, p:list[str]|str)->list[sym.Symbol]:
        if p is None:
            return [],[],0
        if isinstance(p,str):
            p = self.split(p)
        p = [Variable(ip) for ip in p]
        p1 = [E(ip.lead(1),'t') for ip in p]
        if self.indexed:
            p = self.expand_indices(p)
            p1 = self.expand_indices(p1)
        p = [ip.sympy for ip in p]
        p1 = [ip1.sympy for ip1 in p1]
        n_p = len(p)
        return (p,p1,n_p)
    
    def read_z(self, z:list[str]|str)->list[sym.Symbol]:
        if z is None:
            return [],[],0
        if isinstance(z,str):
            z = self.split(z)
        z = [Variable(iz) for iz in z]
        if self.indexed:
            z = self.expand_indices(z)
        z = [iz.sympy for iz in z]
        n_z = len(z)
        return (z,n_z)

    def classify_system(self):
        if self.x==[]:
            return 'forward-looking'
        elif self.p==[]:
            return 'backward-looking'
        else:
            return 'mixed'

    def simulate(self, z:dict[np.array], x0: np.array=None, T:int=None):
        '''
        Simulates for a given path of shocks and initial conditions for predetermined variables.
        x0: initial conditions. Set to 0 if not specified.
        
        '''
        if x0 is None:
            x0 = np.zeros(self.n_x)        
        if T is not None:
            z={k:np.array([v]+[0]*(T-1)) for k,v, in z.items()}
        if T is None:
            T=[len(iz) for iz in z.values()][0]
        z = {str(Variable(k)):v for k,v in z.items()} #normalize the key
        for iz in self.z:
            if str(iz) not in z:
                z[str(iz)]=np.zeros(T)
        z = np.array([z[str(iz)] for iz in self.z])
        
        if self.type=='mixed':
            return self.simulate_mixed_system(z, x0)
        
        if self.type=='backward-looking':
            return self.simulate_backward_looking_system(z, x0)
        
        if self.type=='forward-looking':
            return self.simulate_forward_looking_system(z)

    def simulate_mixed_system(self, z:dict[np.array], x0: np.array):
        sol = self.system_solution
        Theta_x, Theta_p, N, L = sol['Theta_x'], sol['Theta_p'], sol['N'], sol['L']            
        T = z.shape[1]
        x = np.zeros((self.n_x,T+1))
        x[:,0] = x0
        p=np.zeros((self.n_p,T))
        for t,iz in enumerate(z.T):
            iz = iz.reshape(self.n_z,-1)
            x[:,[t+1]] = Theta_x@x[:,[t]]+L@iz
            p[:,[t]] = Theta_p@x[:,[t]]+N@iz
        z = {str(iz):iz_t for iz, iz_t in zip(self.z, z)}
        x = {str(ix):ix_t for ix, ix_t in zip(self.x, x[:,:-1])}
        p = {str(ip):ip_t for ip, ip_t in zip(self.p, p)}
        return MITShock(z, x, p)

    def simulate_backward_looking_system(self, z:dict[np.array], x0: np.array):
        sol = self.system_solution
        Theta_x, L = sol['Theta_x'], sol['L']
        T = z.shape[1]
        x = np.zeros((self.n_x,T+1))
        x[:,0] = x0
        for t,iz in enumerate(z.T):
            iz = iz.reshape(self.n_z,-1)
            x[:,[t+1]] = Theta_x@x[:,[t]]+L@iz
        z = {str(iz):iz_t for iz, iz_t in zip(self.z, z)}
        x = {str(ix):ix_t for ix, ix_t in zip(self.x, x[:,:-1])}
        return MITShock(z, x)

    def simulate_forward_looking_system(self, z:dict[np.array]):
        raise ValueError('Forward Looking not implemented')


class MITShock:
    def __init__(self, z:dict[np.ndarray], x:dict[np.ndarray] = {}, p:dict[np.ndarray] = {}):
        self.paths = z|x|p|{'t':list(range([len(iz) for iz in z.values()][0]))}

    def plot(self, vars:str=None):
        '''
        Plots the paths for the specified variables
        '''
        vars=[str(Variable(i)) for i in vars.replace(' ','').split(',')] 
        nrows = len(vars)//3+1*(len(vars)%3!=0)
        ncols = min(len(vars),3)
        fig, ax = plt.subplots(nrows=nrows, 
                               ncols=ncols, 
                               figsize=(ncols*5, nrows*3),
                               squeeze=False)

        for i,var in enumerate(vars):
            ax[i//ncols,i%ncols].plot(self.paths['t'], self.paths[var])
            ax[i//ncols,i%ncols].set(title=fr'${var}$', xlabel=r'$t$')
        plt.tight_layout()
