from __future__ import annotations
import re
from .utils import normalize_string
from sympy.core.sympify import sympify
from sympy import Symbol, Expr
from IPython.display import Latex
import numpy as np

class Variable:
    '''
    Represents a time-indexed variable, with an optional expectation term
    Attributes
    ----------
    e_t (Expr): time at which expectation is taken
    base (Symbol): name of the variable
    indices (list[Expr]): indices, the last one representing time.
    value (float): 
    '''
    def __init__(self, name:str, value:float=None):
        name = normalize_string(str(name))
        self.e_t, self.base, self.indices = self.split(name)
        self.t = self.indices[-1]
        self.value = value
        self.expectation = self.e_t is not None
        self.realized = value is not None
        self.sympy = Symbol(str(self))
        self.indexed = len(self.indices)>1

    @staticmethod
    def split(name:str)->tuple[Expr, Symbol, list[Expr]]:
        '''
        >>> split('E_{t+1}[x^{p}_{t+2})')
        (t+1, x^{p}, [t+2])
        >>> split('x^{p}_{i,t+2})'
        (None, x^{p}, [i,t+2])
        '''
        e_t = None
        if re.search('E_{', name) is not None:
            e_t = sympify(re.search('(?<=E_{).+?(?=})', name).group())
            name = re.search('(?<=\[).*(?=\])',name).group()
        base =  Symbol(re.sub('_{.*?}','',name))
        if re.search('(?<=_{).+?(?=})',name) is None:
            raise ValueError('Variable needs to have at least one index')
        indices = re.search('(?<=_{).+?(?=})',name).group()
        indices = re.split(',',indices)
        indices = [sympify(i) for i in indices]
        return e_t, base, indices
    
    def __repr__(self):
        if self.realized:
            return str(self.value)
        indices = [str(i) for i in self.indices]
        if self.expectation:
            return f'E_{{{self.e_t}}}[{self.base}_{{{",".join(indices)}}}]'.replace(' ','')
        else: 
            return f'{self.base}_{{{",".join(indices)}}}'.replace(' ','') 
    
    def eval(self, x:float)->float:
        '''
        Evaluates a variable at a given numeric value.
        >>> Variable('x_{i,t}').subs(4)
        
        '''
        return Variable(str(self), float(x))

    def subs(self, x:dict)->Variable:
        '''
        Evaluate the variable at some numeric index or the entire variable

        Parameters
        ----------
        x:dict
            If it is a dictionary, substitute the indices of the variable
        Returns
        -------
        Variable evaluated at some idex, or a float. 

        Examples
        --------
        >>> Variable(x_{i,t}).subs({'i':2})
        'x_{2,t}'
        >>> Variable(x_{i,t+1}).subs({'t':0})
        'x_{i,1}'
      
        '''
        indices = [i.subs(x) for i in self.indices]
        e_t = self.e_t.subs(x) if self.expectation else None
        return Variable.from_elements(self.base, indices, e_t=e_t)

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

    def __call__(self, t:int)->Variable:
        '''
        Returns the variable evaluated at time t
        >>> Variable('x_{t}')(4)
        x_{4}
        >>> Variable('x_{t}').subs({'t':4})
        x_{4}
        '''
        if self.realized:
            raise ValueError('Trying to evaluate an realized variable')
        return self.subs({'t':t})

    def lag(self, periods:int=1):
        '''
        Returns the variable lagged 2 periods
        >>> Variable('x_{t}').lag(2)
        'x_{t-2}'
        '''
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
            if self.indices is None:
                return Parameter.from_elements(self.base)
            else:
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
    def from_elements(cls, base:Symbol|str, indices:list[Expr]=None):
        base = str(base)
        if indices is None:
            return cls(base)
        else:
            indices = [str(i) for i in indices]
            return cls(f'{base}_{{{",".join(indices)}}}')

def E(x:Variable, t:Expr|str='t'):
    return Variable.from_elements(x.base, x.indices, e_t = t)





    