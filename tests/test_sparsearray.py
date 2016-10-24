import numpy as np
from pydigree.cydigree.datastructures import SparseArray


def test_setitem():
    s = SparseArray(100, 0)

    # Add a first element
    s[10] = 1

    # Add something before the first non-sparse element
    s[5] = 1

    # Add something after the last non-sparse element
    s[20] = 1

    # Add sometihng in between
    s[15] = 2
    assert s.keys() == [5, 10, 15, 20]
    assert s.values() == [1, 1, 2, 1]

    s[15] = 0

    assert s.keys() == [5, 10, 20]
    assert s.values() == [1, 1, 1]

    s[20] = 2

    assert s.keys() == [5, 10, 20]
    assert s.values() == [1, 1, 2]

def test_getitem():
    s = SparseArray(100,0)
    assert s[50] == 0

    s2 = SparseArray(100, 1)
    assert s2[50] == 1

    s[5] = 20
    assert s[5] == 20

def test_getslice():
    s = SparseArray(100,0)
    s[5] = 20
    s[10] = 40

    assert len(s[4:20].container) == 2
    assert list(s[4:20].container.keys()) == [1, 6]
    s2 = SparseArray(100, 0)
    s2[1] = 2
    s2[99] = 2
    s2[5] = 20
    s2[10] = 40

    assert len(s2[4:20].container) == 2

def test_fancy_index_get():
    # Bool masks
    s = SparseArray(5,0)
    s[3] = 1
    s[(False, False, True, False, True)].tolist() == [1,0]

    try:
        s[(False, True)] 
        assert False
    except IndexError:
        # We wanted this error
        pass

    # Number masks
    s = SparseArray(10,0)
    mask = [False, False, False, True, False, True, False, False, False, False]
    s[3] = 2
    s[5] = 3
    assert s.tolist() == [0, 0, 0, 2, 0, 3, 0, 0, 0, 0]
    assert s[(3,5)].tolist() == [2,3]
    assert s[(5,3)].tolist() == [3,2]


def test_fancy_index_set():    
    # Bool masks
    s = SparseArray(10,0)
    mask = [False, False, False, True, False, True, False, False, False, False]
    s[mask] = [2,3]
    assert s.tolist() == [0, 0, 0, 2, 0, 3, 0, 0, 0, 0]
    
    s[(3,5)] = [10, 10]
    assert s.tolist() == [0, 0, 0, 10, 0, 10, 0, 0, 0, 0]

    s[(5,3)] = [20, 20]
    assert s.tolist() == [0, 0, 0, 20, 0, 20, 0, 0, 0, 0]

    s = SparseArray(10,0)
    mask = [False, False, False, True, False, True, False, False, False, False]
    s[mask] = 2
    assert s.tolist() == [0, 0, 0, 2, 0, 2, 0, 0, 0, 0]

def test_tolist():
    s = SparseArray(5, 0)
    assert s.tolist() == [0,0,0,0,0]

    s = SparseArray(8, 1)
    s[1] = 2
    s[5] = 3
    assert s.tolist() == [1, 2, 1, 1, 1, 3, 1, 1]

def test_setslice():
    s = SparseArray(100, 0)
    s[1] = 2
    s[99] = 2

    s[5:8] = 3
    assert len(s.container) == 5
    assert s.keys() == [1, 5, 6, 7, 99]
    assert s.values() == [2, 3, 3, 3, 2]

    s = SparseArray(100, 0)
    s[1] = 2
    s[99] = 2
    s[5:8] = [3,3,3]
    assert len(s.container) == 5
    assert s.keys() == [1, 5, 6, 7, 99]
    assert s.values() == [2, 3, 3, 3, 2]

    s = SparseArray(100, 0)
    t = SparseArray(100,0)
    t[5:10] = 1
    s[5:10] = t[5:10]
    assert len(t.container) == 5
    assert t.values() == [1]*5
    assert t.keys() == [5, 6, 7, 8, 9]

    s = SparseArray(10, 0)
    s[2:5] = SparseArray.from_dense([1]*10, 0)[2:5]
    assert s.tolist() == [0,0,1,1,1,0,0,0,0,0]

    s = SparseArray(100,0)
    t = SparseArray(100,0)
    t[20:50] = 1
    s[20:50] = t[20:50]
    assert all(sv == 1 for sv in s[20:50].tolist())

def test_cmp():
    s = SparseArray(5, 0)
    s[(1,3)] = [1, 1]
    assert (s == 1).tolist() == [False, True, False, True, False]
    assert (s > 0).tolist() == [False, True, False, True, False]
    assert (s < 1).tolist() == [True, False, True, False, True]

    assert (s == [0,1,0,1,0]).tolist() == [True, True, True, True, True]

def test_logic():
    s = SparseArray(5,0)
    assert not s.all()
    assert not s.any()

    s[0] = 1
    assert not s.all()
    assert s.any()

    s[0:5] = 1
    assert s.all()
    assert s.any()

    s = SparseArray(5, False)
    s[(1,3)] = True
    assert s.logical_not().tolist() == [True, False, True, False, True]

def test_misc():
    # Sparsity function
    s = SparseArray(10, 0)
    assert s.sparsity() == 1.0
    s[0] = 1
    assert s.sparsity() == 0.9
    s[1] = 1
    assert s.sparsity() == 0.8

    # items() generator
    s = SparseArray(10, 0)
    s[5] = 1
    s[2] = 4
    assert list(s.items()) == [(2,4), (5,1)]


def test_staticbuilders():
    a = SparseArray.from_dense([0, 1, 1, 0, 0, 1], 0)
    assert type(a) is SparseArray
    assert a.tolist() == [0, 1, 1, 0, 0, 1]
    assert a.size == 6
    
    npa = np.array([0, 1, 1, 0, 0, 1], dtype=np.int8)
    assert type(a) is SparseArray
    assert a.tolist() == [0, 1, 1, 0, 0, 1]
    assert a.size == 6

    a = SparseArray.from_items([(0,1), (2,1), (4,1)], 5, 0)
    assert type(a) is SparseArray
    assert a.tolist() == [1,0,1,0,1]
    assert a.size == 5

def test_copy():
    a = SparseArray.from_items([(0,1), (2,1), (4,1)], 5, 0)
    b = a.copy()
    assert a.tolist() == b.tolist()
    a[1] = 1
    assert a.tolist() != b.tolist()
