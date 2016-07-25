import numpy as np

import pydigree
from itertools import combinations
from pydigree.sgs import make_intervals
from pydigree.sgs import SGSAnalysis, SGS, Segment


def to_Segment(ind1, ind2, i):
    ''' i is a 2-tuple of stop and start '''
    return Segment(ind1, ind2, 0, i[0], i[1])


def to_SGSAnalysis(shared):
    ''' Convert the dicts I've been using to SGSAnalysis '''
    obj = SGSAnalysis()
    for key, value in list(shared.items()):
        k = list(key)
        ind1, ind2 = k
        obj[key] = SGS(ind1, ind2,
                       segments=[to_Segment(ind1, ind2, x) for x in value])
    return obj




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
    sgs = to_SGSAnalysis(shared)
    loc = 0, 5
    assert sgs.ibd_state('a', 'b', loc) == 1
    assert sgs.ibd_state('a', 'c', loc) == 2
    assert sgs.ibd_state('b', 'c', loc) == 0
    assert sgs.ibd_state('d', 'e', loc) == 0



def test_ibd_matrix():
    shared = {
        frozenset(['a', 'b']): [(0, 10)],
        frozenset(['a', 'c']): [(0, 10), (0, 10)],
        frozenset(['b', 'c']): []
    }
    sgs = to_SGSAnalysis(shared)
    loc = 0, 0
    # Make sure values are OK
    m1 = np.matrix([[1.0, 0.5, 1.0],
                    [0.5, 1.0, 0.0],
                    [1.0, 0.0, 1.0]])
    assert (sgs.ibd_matrix(['a', 'b', 'c'], loc) == m1).all()

    # Make sure matrix can be reordered
    m2 = np.matrix([[1.0, 0.0, 1.0],
                    [0.0, 1.0, 0.5],
                    [1.0, 0.5, 1.0]])
    assert (sgs.ibd_matrix(['c', 'b', 'a'], loc) == m2).all()

    # Make sure 0's come in for unobserved Inds
    m3 = np.matrix(np.eye(3))
    assert (sgs.ibd_matrix([1, 2, 3], loc) == m3).all()

def test_intervals_to_array():
    from pydigree.sgs.sgs import intervals_to_array
    intervals = [[2,4], [6,10] ] 
    observed = intervals_to_array(intervals, 15)
    expected = np.array([0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0])
    assert all(observed==expected)

    assert all(intervals_to_array([], 10) == np.zeros(10)) 