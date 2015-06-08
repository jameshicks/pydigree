import numpy as np

from itertools import combinations
from pydigree.sgs import nshares, make_intervals

def test_nshares():
    shared = { frozenset([1,2]): [(1,5)], frozenset([1,3]): [(1,3)], frozenset([2,3]): [(1,2)] }
    affecteds = range(1,4)
    nmark = 6
    assert all(nshares(affecteds,shared,nmark) == np.array([0, 3, 3, 2, 1, 1], dtype=np.uint16))

def test_makeintervals():
	assert make_intervals([0]*1000) == []
	assert make_intervals([1]*1000) == [(0,999)]
	assert make_intervals([2]*1000) == [(0,999), (0,999)]
	assert make_intervals([0,0,0,1,1,1,0,0]) == [(3,5)]
	assert sorted(make_intervals([0,1,1,2,2,2,1,1,0])) == [(1,7), (3,5)]

