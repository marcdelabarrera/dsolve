from multiprocessing.sharedctypes import Value
from .atoms import Variable, Parameter
from .utils import normalize_string, normalize_dict
import re
import numpy as np
from sympy import Eq, Expr, Symbol
import sympy as sym

class DynamicExpression:
    def __init__(self, expression:str|sym.Expr):
        self.elements = split(str(expression))
        self.variables = {str(Variable(i)):Variable(i) for i in self.elements if is_variable(i)}
        self.parameters = {str(Parameter(i)):Parameter(i) for i in self.elements if is_parameter(i)}
        self.indexed = np.any([v.indexed for v in self.variables.values()])

    @property
    def sympy(self):
        elements, variables, parameters = self.elements, self.variables, self.parameters
        elements = [f'variables[r"{i}"].sympy' if i in variables else i for i in elements]
        elements = [f'parameters[r"{i}"].sympy' if i in parameters else i for i in elements]
        return eval(''.join(elements))
    
    def __float__(self):
        if self.variables=={} and self.parameters=={}:
            return float(self.sympy)
        else:
            raise ValueError('Trying to convert to a float a DynamicExpression with unevaluated parameters or variables')

    def __repr__(self):
        return str(self.sympy).replace(' ','')

    def __call__(self, t):
        '''
        >>> DynamicExpression('x_{t}+y_{t+1}')(3)
        x_{3}+y_{4}
        '''
        eq = ''.join([str(Variable(el)(t)) if is_variable(el) else el for el in self.elements])
        return DynamicExpression(''.join(eq))

    def lag(self, periods:int=1):
        eq = ''.join([str(Variable(el).lag(periods)) if is_variable(el) else el for el in self.elements])
        return DynamicExpression(''.join(eq))

    def lead(self, periods:int=1):
        return self.lag(-periods)

    def subs(self, d:dict):
        '''
        >>> DynamicExpression('x_t+y_t').subs({'x_{t}':4})
        4+y_{t}
        >>> DynamicExpression('x_t+y_t').subs({'t':4})
        x_{4}+y_{4}
        '''
        d = normalize_dict(d)
        expr = []
        for el in self.elements:
            if is_variable(el):
                if el in d.keys():
                    expr.append(str(d[el]))
                else:
                    expr.append(str(Variable(el).subs(d)))
            elif is_parameter(el):
                if el in d.keys():
                    expr.append(str(d[el]))
                else:
                    expr.append(str(Parameter(el).subs(d)))
            else:
                expr.append(el)
        return DynamicExpression(''.join(expr))

class DynamicEquation(DynamicExpression):
    
    @classmethod
    def from_sympy(cls, eq:Eq):
        return cls(f'{eq.lhs}={eq.rhs}')

    @classmethod
    def from_lhs_rhs(cls, lhs, rhs):
        return cls(f'{lhs}={rhs}')


    def __repr__(self):
        rhs, lhs = self.sympy.rhs, self.sympy.lhs
        return f'{rhs} = {lhs}'

    @property
    def sympy(self):
        elements, variables, parameters = self.elements, self.variables, self.parameters
        elements = [f'variables[r"{i}"].sympy' if i in variables else i for i in elements]
        elements = [f'parameters[r"{i}"].sympy' if i in parameters else i for i in elements]
        lhs = ''.join(elements[:elements.index('=')])
        rhs = ''.join(elements[elements.index('=')+1:])
        return Eq(eval(lhs),eval(rhs))
    
    @property
    def lhs(self):
        return self.sympy.lhs

    @property
    def rhs(self):
        return self.sympy.rhs

    @property
    def free_symbols(self):
        return self.sympy.free_symbols

    def __call__(self, t):
        eq = ''.join([str(Variable(el)(t)) if is_variable(el) else el for el in self.elements])
        return DynamicEquation(''.join(eq))

    def lag(self, periods:int=1):
        eq = ''.join([str(Variable(el).lag(periods)) if is_variable(el) else el for el in self.elements])
        return DynamicEquation(''.join(eq))

    def lead(self, periods:int=1):
        return self.lag(-periods)

    def subs(self, d:dict):
        return DynamicEquation.from_lhs_rhs(DynamicExpression(self.rhs).subs(d), DynamicExpression(self.lhs).subs(d))

        return 
        if np.all([Symbol(k) in self.free_symbols for k in d.keys()]):
            return DynamicEquation.from_sympy(self.sympy.subs(d))
        else: 
            eq = []
            for el in self.elements:
                if is_variable(el):
                    eq.append(str(Variable(el).subs(d)))
                elif is_parameter(el):
                    eq.append(str(Parameter(el).subs(d)))
                else:
                    eq.append(el)
            return DynamicEquation(''.join(eq))
    




def close_brackets(elements:list[str])->list[str]:
    '''
    Given a list of elements, ensures that brackets are closed
    Examples
    --------
    >>> close_brackets(['x_{','t','}'])
    ['x_{t}']
    '''
    out=[]
    elements = iter(elements)
    for i in elements:
        out.append(i)
        if '{' in i:
            while len(re.findall('{',out[-1]))!=len(re.findall('}',out[-1])):
                out[-1] = out[-1]+next(elements)
    return out

def classify_string(string):
    string = normalize_string(string)
    if string[:4]=='\sum':
        return 'sum'
    elif re.match('^\\\\frac{', string) is not None:
        return 'fraction'
    elif str.isdigit(string.replace('.','')):
        return 'number'
    elif re.search('_{[^\\\]*t.*}', string) is not None:
        return 'variable'
    elif string in ('+','-','/','*','(',')','='):
        return 'operator'
    else:
        return 'parameter'
    
def is_sum(string)->bool:
    return classify_string(string)=='sum'

def is_number(string)->bool:
    return classify_string(string)=='number'

def is_variable(string)->bool:
    '''
    Returns true if string can be parsed into a variable.
    String contains '_{[^\\\]*t.*}'
    Evaluated variables are not considered variables.
    >>> is_variable('x_{t}')
    True
    >>> is_variable('\sigma')
    False
    >>> is_variable('x_{3}')
    False
    '''
    return classify_string(string)=='variable'

def is_parameter(string)->bool:
    return classify_string(string)=='parameter'

def is_fraction(string)->bool:
    return classify_string(string)=='fraction'

def split_fraction(frac:str)->list[str]:
    if not is_fraction(frac):
        raise ValueError('frac must start with \\frac')
    frac = re.split('(?<={)|(?=})',frac) 
    frac[0]='('
    frac[-1]=')'
    open_bracket = np.cumsum(np.array(['{' in i for i in frac]))
    close_bracket = np.cumsum(np.array(['}' in i for i in frac]))
    for i in range(len(frac)):
        if frac[i]=='}{' and (open_bracket-close_bracket)[i]==0:
            frac[i]=')/('
    return split(''.join(frac))

def split_sum(sum:str)-> list[str]:
    if not is_sum(sum):
        raise ValueError('sum must start with \sum')
    index = re.search("(?<=sum\_{).+?(?=\=)",sum).group()
    start = int(re.search("(?<=sum\_{.=).+?(?=})",sum).group())
    end = int(re.search("(?<=\^{).+?(?=})",sum).group())
    term = re.sub('\\\sum_{.+?}\^{.+?}','',sum)[1:-1]
    term = DynamicExpression(term)
    return split(f'({"+".join([str(term.subs({index:i})) for i in range(start,end+1)])})')   

def split(expression:str)->list[str]:
    elements = re.split('(?<=[\=\*/\+\-\(\)])|(?=[\=\*/\+\-\(\)])',expression)
    elements = [i for i in elements if i.replace(' ','')]
    elements = close_brackets(elements)
    out = []
    for el in elements:
        if is_fraction(el):
            out += split_fraction(el)
        elif is_sum(el):
            out += split_sum(el)
        elif is_variable(el):
            out += [str(Variable(el))]
        elif is_parameter(el):
            out += [str(Parameter(el))]
        else:
            out += [el]
    return out
