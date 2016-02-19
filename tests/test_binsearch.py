from nose.tools import assert_raises
from pydigree.datastructures.algo import binsearch
from pydigree.datastructures.algo import binsearch_left
from pydigree.datastructures.algo import binsearch_right

def test_binsearch():
    a = [1,3,5,7]
    assert binsearch(a, 3) == 1
    assert binsearch(a, 7) == 3
    assert binsearch(a, 1) == 0
    assert_raises(KeyError, binsearch, a, -1)
    assert_raises(KeyError, binsearch, a, 100)
    assert_raises(KeyError, binsearch, a, 4)
    assert_raises(KeyError, binsearch, [], 1)

    assert binsearch([0,2,4,6], 3, key=lambda x:x+1)


def test_binsearch_left():
    a = [1,3,5,7]
    assert binsearch_left(a, 3) == 1
    assert binsearch_left(a, 0) == 0
    assert binsearch_left(a, 6) == 3
    assert binsearch_left([], 1) == 0

def test_binsearch_right():
    a = [1,3,5,7]
    assert binsearch_right(a, 3) == 2
    assert binsearch_right(a, 0) == 0
    assert binsearch_right(a, 6) == 3
    assert binsearch_right(a, 7) == 4
    assert binsearch_right([], 1) == 0