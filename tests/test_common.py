from pydigree.common import product, cumsum, flatten, invert_dict
from pydigree.common import runs, runs_gte
from pydigree.cyfuncs import runs_gte_uint8

import numpy as np

def test_product():
    assert product([1,2,3,4,5]) == 120
    assert product([0,1,2,3,4]) ==  0


def test_cumsums():
    assert cumsum([]) == []
    assert cumsum([1,2,3]) == [1,3,6]

def test_runs():
    dregion = (10,20)
    i = [1 if dregion[0] <= x <= dregion[1] else 0 for x in xrange(100)] 
    assert runs(i, lambda x: x>0) == [dregion]
    assert runs_gte(i, 1) == [dregion]
    assert runs_gte(i, 1, 20) == []
    assert runs_gte([], 1) == []
    assert runs([], lambda x: True) == []
    

    i = [0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1]
    ir = [(4,8), (14, 16)]
    assert runs(i, lambda x: x > 0) == ir
    assert runs_gte(i, 1) == ir
    assert runs_gte_uint8(np.array(i, dtype=np.uint8), 1) == ir

    i = [0] * 100
    assert runs(i, lambda x: x > 0) == []
    assert runs_gte(i, 1) == []
    assert runs_gte_uint8(np.uint8(i),1) == []

def test_flatten():
    assert list(flatten([1,2,2,3,4,4,5])) == [1,2,2,3,4,4,5]
    assert list(flatten([(1,2),2,3,(4,4,5)])) == [1,2,2,3,4,4,5]
    assert list(flatten([(1,(2,(2,(3,4)))),4,5])) == [1,2,2,3,4,4,5]
    assert list(flatten([])) == []

def test_invertdict():
    assert invert_dict({}) == {}
    assert invert_dict({1: 'a', 2: 'b'}) == {'a': 1, 'b': 2}
