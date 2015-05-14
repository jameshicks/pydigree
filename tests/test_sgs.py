import numpy as np

from itertools import combinations
from pydigree.sgs import nshares

def test_nshares():
    shared = { frozenset([1,2]): [(1,5)], frozenset([1,3]): [(1,3)], frozenset([2,3]): [(1,2)] }
    affecteds = range(1,4)
    nmark = 6
    assert all(nshares(affecteds,shared,nmark) == np.array([0, 3, 3, 2, 1, 1], dtype=np.uint16))