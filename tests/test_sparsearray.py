from pydigree.cyfuncs import SparseArray, SparseArrayElement

def SLE_eq():
	a = SparseArrayElement(1,1)
	b = SparseArrayElement(1,1)
	c = SparseArrayElement(2,2)
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

def test_setslice():
	s = SparseArray(100, 0)
	s[1] = 2
	s[99] = 2

	s[5:8] = 3
	assert len(s.container) == 5
	assert s.indices == [1, 5, 6, 7, 99]
	assert s.values == [2, 3, 3, 3, 2]


def test_staticbuilted():
    s = SparseArray(100, 0)
    a = SparseArray.from_sequence([0, 1, 1, 0, 0, 1], 0)
    assert type(a) is SparseArray
    assert a.indices == [1, 2, 5]
    assert a.values == [1, 1, 1]
    
