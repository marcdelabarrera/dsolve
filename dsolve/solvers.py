from expressions import DynamicEcuation, close_brackets
from atoms import Variable, Parameter, E
import re
import numpy as np
import sympy as sym

class Klein():

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

    def __init__(self, 
                    equations:list[str]|str, 
                    x: list[str]|str = None, 
                    p: list[str]|str = None, 
                    z: list[str]|str = None, 
                    indices: dict[list[int]]= None,
                    calibration:dict = None):

        if len(indices)>1:
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
        if calibration is not None:
             self.system_numeric = self.calibrate(calibration)
        
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