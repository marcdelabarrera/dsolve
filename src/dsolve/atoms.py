import re
from sympy.core.sympify import sympify
from sympy import Symbol, Expr
from IPython.display import Latex

class Variable:
    def __init__(self, name:str):
        name = normalize_string(name)
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


def normalize_string(string)->str:
    '''
    Normalizes strings adding {} when ommited.
    Examples
    --------
    >>> normalize_string('a_b^c')
    'a_{b}^{c}'
    '''
    string = encode(string)
    string = string.replace(' ','')                                     #gets rid of any space
    string = re.sub('^E(?!_)','E_{t}', string)                          #correct E[.] to add E_{t}[.]
    string = re.sub('(?<=_)[^{]',lambda m: f'{{{m.group()}}}',string)   #correct x_. to x_{.}
    string = re.sub('(?<=\^)[^{]',lambda m: f'{{{m.group()}}}',string)  #correct x^. to x^{.}
    if re.search('^E_{[^}]+?}(?!\[)', string) is not None:
        string = re.search('^E_{[^}]+?}(?!\[)', string).group()+'['+re.sub('^E_{[^}]+?}','',string)+']' #correct E_{t}. for E_{t}[.]
    return string


def encode(string)->str:
    '''
    Corrects for special characters and adds an extra \
    Examples
    --------
    >>> encode('\beta')
    '\\beta'
    '''
    string = string.replace('\b','\\b')
    string = string.replace('\r','\\r')
    string = string.replace('\v','\\v')
    string = string.replace('\f','\\f')
    string = string.replace('\t','\\t')
    return string