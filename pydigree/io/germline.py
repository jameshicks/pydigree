import numpy as np

from pydigree.common import runs
from pydigree.io.fastio import read_germline

def is_germline_haploid(filename):
    """ Returns true if file is a germline --haploid file """
    with open(filename) as f:
        l = f.readline().strip().split()
        ind = l[1]
        return ind.endswith('.0') or ind.endswith('.1')

def apply_inheritance_model(shared, nmark, model=None):
    """ 
    Applies an inheritance model to a the output of read_germline 
    Valid models:
    'add': No model is applied. Input is returned as output
    'dom': Returns runs that are IBD > 1 
    'rec': Returns runs that are IBD == 2
    """
    
    if not model or model in {'add','pairs'}:
        return shared
    elif model in {'dom', 'bool'}:
        pred = lambda x: x > 1
    elif model == 'rec':
        pred = lambda x: x == 2
    else:
        raise ValueError('Invalid inheritance model: {}'.format(model))
    
    for pair in shared:
        s = np.zeros(nmark)
        
        for start, stop in shared[pair]:
            s[start:(stop+1)] += 1
        
        shared[pair] = [x for x in runs(s, pred)]

 
