import numpy as np

from itertools import combinations
from pydigree.sgs import nshares, make_intervals, ibd_state, ibd_matrix


def test_nshares():
    shared = {frozenset([1, 2]): [(1, 5)], frozenset(
        [1, 3]): [(1, 3)], frozenset([2, 3]): [(1, 2)]}
    affecteds = range(1, 4)
    nmark = 6
    assert all(nshares(affecteds, shared, nmark) == np.array(
        [0, 3, 3, 2, 1, 1], dtype=np.uint16))

    shared = {}
    affecteds = [1, 2, 3, 4]
    nmark = 6
    assert all(nshares(affecteds, shared, 6) == np.zeros(6))

    shared = {frozenset([1, 2]): [(1, 5), (1, 5)]}
    affecteds = range(1, 4)
    nmark = 6
    assert all(nshares(affecteds, shared, nmark) == np.array(
        [0, 2, 2, 2, 2, 2], dtype=np.uint16))


def test_makeintervals():
    assert make_intervals([0]*1000) == []
    assert make_intervals([1]*1000) == [(0, 999)]
    assert make_intervals([2]*1000) == [(0, 999), (0, 999)]
    assert make_intervals([0, 0, 0, 1, 1, 1, 0, 0]) == [(3, 5)]
    assert sorted(make_intervals([0, 1, 1, 2, 2, 2, 1, 1, 0])) == [
        (1, 7), (3, 5)]


def test_ibd_state():
    shared = {
        frozenset(['a', 'b']): [(0, 10)],
        frozenset(['a', 'c']): [(0, 10), (0, 10)],
        frozenset(['b', 'c']): []
    }
    assert ibd_state(shared, 'a', 'b', 5) == 1
    assert ibd_state(shared, 'a', 'c', 5) == 2
    assert ibd_state(shared, 'b', 'c', 5) == 0
    assert ibd_state(shared, 'd', 'e', 5) == 0

def test_ibd_matrix():
    shared = {
        frozenset(['a', 'b']): [(0, 10)],
        frozenset(['a', 'c']): [(0, 10), (0, 10)],
        frozenset(['b', 'c']): []
    }

    # Make sure values are OK
    m1 = np.matrix([[0, 1 ,2], [1, 0 ,0], [2, 0, 0]])
    assert (ibd_matrix(shared, ['a','b','c'], 0) == m1).all()

    # Make sure matrix can be reordered
    m2 = np.matrix([[0, 0, 2], [0, 0, 1], [2, 1, 0]])
    assert (ibd_matrix(shared, ['c','b','a'], 0) == m2).all()

    # Make sure 0's come in for unobserved Inds
    m3 = np.matrix(np.zeros((3,3)))
    assert (ibd_matrix(shared, [1, 2, 3], 0) == m3).all()


