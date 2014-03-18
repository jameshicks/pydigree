import numpy as np
from pydigree._pydigree import sgs_shares

def proportion_shares(affecteds, shared, nmark):
    affecteds = list(affecteds)
    s = sgs_shares(affecteds, shared, nmark);
    s = np.array(s)

    n = len(affecteds)
    return s / (n * (n-1) * 0.5)
