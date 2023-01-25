#General methods used by different classes
import re

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
    string = string.replace('\a','\\a')
    return string

def normalize_string(string)->str:
    '''
    Normalizes strings adding {} when ommited.
    Examples
    --------
    >>> normalize_string('a_b^c')
    'a_{b}^{c}'

    >>> normalize_string('Ex_{t+1}')
    'E_{t}[x_{t+1}]'
    '''
    string = encode(string)
    string = string.replace(' ','')                                     #gets rid of any space
    string = re.sub('E(?!_)','E_{t}', string)                           #correct E[.] to add E_{t}[.]
    string = re.sub('(?<=_)[^{]',lambda m: f'{{{m.group()}}}',string)   #correct x_. to x_{.}
    string = re.sub('(?<=\^)[^{]',lambda m: f'{{{m.group()}}}',string)  #correct x^. to x^{.}
    if re.search('E_{[^}]+?}(?!\[)', string) is not None:
        string = re.search('E_{.+?}', string).group()+'['+re.sub('E_{.+?}','',string)+']' #correct E_{t}. for E_{t}[.]
    return string

def normalize_dict(d)->dict:
    if d is None:
        return None
    return {normalize_string(k):v for k,v in d.items()}