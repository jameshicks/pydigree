from pydigree.cyfuncs import SparseArray, SparseArrayElement


def SLE_eq():
    a = SparseArrayElement(1, 1)
    b = SparseArrayElement(1, 1)
    c = SparseArrayElement(2, 2)
    assert a == b
    assert not a == c
    assert a != c


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
    assert s.indices == [5, 10, 15, 20]
    assert s.values == [1, 1, 2, 1]

    s[15] = 0

    assert s.indices == [5, 10, 20]
    assert s.values == [1, 1, 1]

    s[20] = 2

    assert s.indices == [5, 10, 20]
    assert s.values == [1, 1, 2]

def test_getslice():
    s = SparseArray(100,0)
    s[5] = 20
    s[10] = 40

    assert len(s[4:20].container) == 2

    s2 = SparseArray(100, 0)
    s2[1] = 2
    s2[99] = 2
    s2[5] = 20
    s2[10] = 40

    assert len(s2[4:20].container) == 2

def test_fancy_index_get():
    s = SparseArray(5,0)
    s[3] = 1
    s[(False, False, True, False, True)].tolist() == [1,0]

    try:
        s[(False, True)] 
        assert False
    except IndexError:
        # We wanted this error
        pass

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
    assert s.indices == [1, 5, 6, 7, 99]
    assert s.values == [2, 3, 3, 3, 2]

    s = SparseArray(100, 0)
    t = SparseArray(100,0)
    t[5:10] = 1
    s[5:10] = t[5:10]
    assert len(t.container) == 5
    assert t.values == [1]*5
    assert t.indices == [5, 6, 7, 8, 9]

def test_staticbuilted():
    s = SparseArray(100, 0)
    a = SparseArray.from_sequence([0, 1, 1, 0, 0, 1], 0)
    assert type(a) is SparseArray
    assert a.indices == [1, 2, 5]
    assert a.values == [1, 1, 1]
    
