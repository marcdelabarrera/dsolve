
from __future__ import annotations
import sympy as sym
from .expressions import DynamicEquation, close_brackets, DynamicExpression
from dataclasses import dataclass, field
from .utils import normalize_string, normalize_dict
import numpy as np
from .atoms import Variable, E, Parameter
import re

@dataclass
class SystemVariables:
    '''
    Class that reads and contains all information regarding system variables.
    '''
    x : list[sym.Symbol] = field(default_factory=list)
    z : list[sym.Symbol] = field(default_factory=list)
    s : list[sym.Symbol] = field(default_factory=list)
 
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

class StateSpace:
    def __init__(self, equations, x,s,z, indices):
        indices = self.read_indices(indices)
        equations = [DynamicEquation(eq) for eq in equations]
        if indices is not None:
            equations = self.expand_indices(equations, indices)
        self.equations = equations
        self.vars = self.read_variables(x,z,s,indices)

    def read_variables(self, x,z,s, indices:SystemIndices=None)->SystemVariables:
        x,z,s = [self.split(i) if isinstance(i,str) else i for i in [x,z,s]]
        x,z,s = [[Variable(j) for j in i] if i is not None else [] for i in [x,z,s]]
        if indices is not None:
            x,z,s = [self.expand_indices(i, indices) for i in [x,z,s]]
        x,z,s = [[j.sympy for j in i] for i in [x,z,s]]
        return SystemVariables(x,z,s)

    @staticmethod
    def split(string)->list[str]:
        l = re.split('(?<=,)|(?=,)',string)
        l = close_brackets(l)
        l = [i for i in l if i!='' and i!=',']
        return l

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

    def read_indices(self, indices:dict[list[int]])->SystemIndices:
        if indices is None:
            return None
        indices, start, end  = list(indices.keys()), [indices[i][0] for i in indices], [indices[i][1] for i in indices]
        indices = SystemIndices(indices, start, end)
        return indices


    def expand_indices(self, l:list=None, indices:SystemIndices=None):
        if indices is None or l is None:
            return l 
        out = l.copy()
        for i in range(len(indices.indices)):
            out = self.expand_index(out, indices.indices[i],indices.start[i],indices.end[i])
        return out
    
    def expand_index(self, l, i, start, stop):
        out = []
        for el in l:
            if sym.Symbol(i) in el.indices:
                out=out+[el.subs({i:j}) for j in range(start,stop+1)]
            else:
                out.append(el)
        return out

    def normalize_z(self, z:dict, T=None)->np.array:
        '''
        
        >>> self.normalize_z({'z_{0}':1.},T=4)
        np.array([1.,0.,0.,0.])

        >>> self.normalize_z({'z_{0}':1.,'z_2':2.},T=4)
        np.array([1.,0.,2.,0.])
        '''
        z = normalize_dict(z)
        out={}
        if T is not None:
            for t in range(T):
                for iz in self.vars.z:
                    out=out|{str(Variable(iz).subs({'t':t})):0}
        return out|z

    def solve(self, T, z, x0):
        z = self.normalize_z(z,T)
        system=[]
        x=set()
        for t in range(T):
            for eq in self.equations:
                eq = eq.subs({'t':t})
                eq = eq.subs(z)
                eq = eq.subs(x0)
                system.append(eq.sympy)
                x = x.union(eq.free_symbols)
        A,b = sym.linear_eq_to_matrix(system, x)      
        sol = A.inv()@b  
        return system, {k:v for k,v in zip(x,sol)}