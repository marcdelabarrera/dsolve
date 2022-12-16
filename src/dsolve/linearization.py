import sympy as sym


def linearize(f:sym.Function, order:float=2, type:str='log-linearize',x:list[sym.Symbol] = None, x_ss:list[float] = None, f_ss:float = None)->sym.core.add.Add:
    '''
    Log-linearizes f up to a second order around a steady state x_ss. 
    Parameters
    ----------
    x_ss: list[float]
        If the steady state is known, it substitutes it to the final expression.
    f_ss: float
        Value of f at the steady state
    '''
    if x is None:
        x = f.args
    n = len(x)
    if x_ss is None:
        x_ss = sym.symbols([f'{str(i)}_ss' for i in x])

    x_hat = sym.symbols([f'\hat{{{str(i)}}}' for i in x])
    if type!='log-linearize':
        raise ValueError('Standard linearization not implemented')

    f = f.subs({i:j*sym.exp(k) for i,j,k in zip(x,x_ss,x_hat)})
    
    if f_ss is None:
        f_ss = f.subs({i:0 for i in x_hat})

    f1 = sum([f.diff(i).subs({j:0 for j in x_hat})*i for i in x_hat])
    
    if order ==1:
        return (f_ss+f1).simplify()

    f2 = sym.zeros(n,n)
    for i,x in enumerate(x_hat):
        for j,y in enumerate(x_hat):
            if x==y:
                f2[i,j]=f.diff(x,x).subs({k:0 for k in x_hat})*x**2
            else:
                f2[i,j]=f.diff(x).subs({x:0}).diff(y).subs({k:0 for k in x_hat})*x*y
                #f2[i,j]=f.diff(x).subs({x:0}).subs({k:0 for k in x_hat})*x*y

    linear = f_ss+f1+1/2*sum(f2).simplify()
    return linear.simplify()
