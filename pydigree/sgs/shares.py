import numpy as np
from pydigree._pydigree import sgs_shares

def nshares(affecteds, shared, nmark):
    affecteds = list(affecteds)
    s = sgs_shares(affecteds, shared, nmark);
    return np.array(s)
