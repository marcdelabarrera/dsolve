import re
from sympy.core.sympify import sympify
from sympy import Symbol, Expr
from IPython.display import Latex

class Variable:
    def __init__(self, name:str):
        self.e_t, self.base, self.indices = self.split(name)
        self.expectation = self.e_t is not None
        self.sympy = Symbol(str(self).replace(' ',''))
        self.indexed = len(self.indices)>1

    @staticmethod
    def split(name:str)->tuple[Symbol, Expr, list[Expr]]:
        e_t = None
        if re.search('E_{', name) is not None:
            e_t = sympify(re.search('(?<=E_{).+?(?=})', name).group())
            name = re.sub('E_{.*?}\[','',name)[:-1]
        base =  Symbol(re.sub('_{.*?}','',name))
        indices = re.search('(?<=_{).+?(?=})',name).group()
        indices = re.split(',',indices)
        indices = [sympify(i) for i in indices]
        return e_t, base, indices
    
    def __repr__(self):
        indices = [str(i).replace(' ','') for i in self.indices]
        if self.expectation:
            e_t = str(self.e_t).replace(' ','')
            return f'E_{{{e_t}}}[{self.base}_{{{",".join(indices)}}}]'
        else: 
            return f'{self.base}_{{{",".join(indices)}}}' 
    
    def subs(self, d:dict):
        indices = [i.subs(d) for i in self.indices]
        return Variable.from_elements(self.base, indices, e_t=self.e_t)

    @property
    def latex(self):
        return(Latex(f'${self}$'))

    @classmethod
    def from_elements(cls, base:Symbol|str, indices:list[Expr], e_t:Expr|str = None):
        base = str(base)
        indices = [str(i) for i in indices]
        if e_t is not None:
            return cls(f'E_{{{e_t}}}[{base}_{{{",".join(indices)}}}]')
        else: 
            return cls(f'{base}_{{{",".join(indices)}}}')

    def __call__(self, t):
        indices = self.indices.copy()
        indices[-1] = indices[-1].subs({'t':t})
        e_t = self.e_t.subs({'t':t}) if self.expectation else None
        return Variable.from_elements(self.base, indices, e_t=e_t)

    def lag(self, periods:int=1):
        indices = self.indices.copy()
        indices[-1] = indices[-1]-periods
        e_t = self.e_t-periods if self.expectation else None
        return Variable.from_elements(self.base, indices, e_t=e_t)

    def lead(self, periods:int=1):
        return self.lag(-periods)

class Parameter:
    def __init__(self, name):
        self.base =  Symbol(re.sub('_{.*?}','',name))
        self.indices = None
        self.indexed = False
        if re.search('(?<=_{).+?(?=})',name) is not None:
            self.indexed = True
            indices = re.search('(?<=_{).+?(?=})',name).group()
            indices = re.split(',',indices)
            if len(indices)>1:
                raise ValueError('Code cannot handle more than one index for parameters')
            self.indices = [sympify(i) for i in indices]
        self.sympy = Symbol(name)
        
    def subs(self, d:dict):
        if self.indexed:
            indices = [i.subs(d) for i in self.indices]
            return Parameter.from_elements(self.base, indices)
        else:
            return Parameter(str(self.sympy))


    def __repr__(self):
        return str(self.sympy)

    @classmethod
    def from_elements(cls, base:Symbol|str, indices:list[Expr]):
        base = str(base)
        indices = [str(i) for i in indices]
        return cls(f'{base}_{{{",".join(indices)}}}')

def E(x:Variable, t:Expr|str):
    return Variable.from_elements(x.base, x.indices, e_t = t)

