from pydigree.datastructures import SortedList
from nose.tools import assert_raises

def test_sortedlist():
    a = [3,1,8,9,2,2,1]
    b = SortedList(a)
    assert b.container == sorted(a)

    assert b[0] == 1

    a = [] 
    b = SortedList(a)
    assert b.container == []
    if b: assert 'bool failed'

def test_extend():
    a = SortedList([1,2,2,3])
    a.extend(SortedList([4,5,6]))
    assert a.container == [1,2,2,3,4,5,6]
    assert_raises(ValueError, a.extend, SortedList([2,3,4]))

def test_append():
    a = SortedList([1,2,3])
    a.append(4)
    assert a.container == [1,2,3,4]
    assert_raises(ValueError, a.append, 2)

def test_clear():
    a = SortedList([1,2,3])
    a.clear()
    assert a.container == []