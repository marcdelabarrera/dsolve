from __future__ import annotations
import re
from .utils import normalize_string
from sympy.core.sympify import sympify
from sympy import Symbol, Expr
from IPython.display import Latex
import numpy as np

class Variable:
    def __init__(self, name:str):
        name = normalize_string(name)
        self.e_t, self.base, self.indices = self.split(name)
        self.expectation = self.e_t is not None
        self.realized = np.all(['t' not in str(i) for i in self.indices])
        self.sympy = Symbol(str(self))
        self.indexed = len(self.indices)>1

    @staticmethod
    def split(name:str)->tuple[Symbol, Expr, list[Expr]]:
        e_t = None
        if re.search('E_{', name) is not None:
            e_t = sympify(re.search('(?<=E_{).+?(?=})', name).group())
            name = re.search('(?<=\[).*(?=\])',name).group()
        base =  Symbol(re.sub('_{.*?}','',name))
        indices = re.search('(?<=_{).+?(?=})',name).group()
        indices = re.split(',',indices)
        indices = [sympify(i) for i in indices]
        return e_t, base, indices
    
    def __repr__(self):
        indices = [str(i) for i in self.indices]
        if self.expectation:
            return f'E_{{{self.e_t}}}[{self.base}_{{{",".join(indices)}}}]'.replace(' ','')
        else: 
            return f'{self.base}_{{{",".join(indices)}}}'.replace(' ','') 
    
    def subs(self, x:dict|float)->Variable|float:
        '''
        >>> Variable(x_{i,t}).subs({'i':2})
        'x_{2,t}'
        >>> Variable(x_{i,t+1}).subs({'t':0})
        'x_{i,1}'
        >>> Variable(x_{i,t}).subs(4)
        4
        '''
        if isinstance(x,dict):
            indices = [i.subs(x) for i in self.indices]
            e_t = self.e_t.subs(x) if self.expectation else None
            return Variable.from_elements(self.base, indices, e_t=e_t)
        else:
            return x

    @property
    def latex(self):
        return(Latex(f'${self}$'))

    @classmethod
    def from_elements(cls, base:Symbol|str, indices:list[Expr], e_t:Expr|str = None):
        indices = [str(i) for i in indices]
        if e_t is not None:
            return cls(f'E_{{{e_t}}}[{base}_{{{",".join(indices)}}}]')
        else: 
            return cls(f'{base}_{{{",".join(indices)}}}')

    def __call__(self, t):
        '''
        Returns the variable evaluated at time t
        >>> Variable('x_{t}')(4)
        'x_{4}'
        '''
        if self.realized:
            raise ValueError('Trying to evaluate an evaluated variable')
        return self.subs({'t':t})

    def lag(self, periods:int=1):
        indices = self.indices.copy()
        indices[-1] = indices[-1]-periods
        e_t = self.e_t-periods if self.expectation else None
        return Variable.from_elements(self.base, indices, e_t=e_t)

    def lead(self, periods:int=1):
        return self.lag(-periods)

class Parameter:
    def __init__(self, name):
        name = normalize_string(name)
        self.base =  Symbol(re.sub('_{.*?}','',name))
        self.indices = None
        self.indexed = False
        if re.search('(?<=_{).+?(?=})',name) is not None:
            self.indexed = True
            indices = re.search('(?<=_{).+?(?=})',name).group()
            indices = re.split(',',indices)
            if len(indices)>1:
                raise ValueError('Code cannot handle more than one index for parameters')
            self.indices = [Symbol(i) for i in indices]
        self.sympy = Symbol(str(self))
  
    def subs(self, x:dict|float):
        '''
        >>> Parameter(\beta_{i}).subs({'i':2})
        '\\beta_{2}'
        >>> Parameter(\sigma_{i}).subs(2.)
        2.
        '''
        if isinstance(x,dict):
            indices = [i.subs(x) for i in self.indices]
            return Parameter.from_elements(self.base, indices)
        else:
            return x
   
    def __repr__(self):
        if self.indexed:
            indices = [str(i) for i in self.indices]
            return f'{self.base}_{{{",".join(indices)}}}'.replace(' ','') 

        else:  
            return f'{self.base}' 


    @classmethod
    def from_elements(cls, base:Symbol|str, indices:list[Expr]):
        base = str(base)
        indices = [str(i) for i in indices]
        return cls(f'{base}_{{{",".join(indices)}}}')

def E(x:Variable, t:Expr|str='t'):
    return Variable.from_elements(x.base, x.indices, e_t = t)





    