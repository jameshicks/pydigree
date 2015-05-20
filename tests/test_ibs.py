from nose.tools import assert_raises
from pydigree.ibs import ibs, chromwide_ibs
from pydigree.genotypes import Alleles, SparseAlleles

import numpy as np

def test_ibs():
    assert ibs( (2,2), (2,2) ) == 2
    assert ibs( (1,2), (1,2) ) == 2
    assert ibs( (2,1), (2,2) ) == 1
    assert ibs( (2,2), (2,1) ) == 1
    assert ibs( (1,1), (2,2) ) == 0
    assert ibs( (1,1), (0,0), missingval=64) == 64
    assert ibs( (1,1), (0,0) ) == None
    assert ibs( (0,0), (1,1) ) == None

def test_chromwide_ibs():
    g1 = [(2,2), (1,2), (1,2), (1,1), (0, 0)]
    g2 = [(2,2), (1,2), (2,2), (2,2), (1, 1)]

    a, b = [Alleles(x) for x in zip(*g1)]
    c, d = [Alleles(x) for x in zip(*g2)]

    expected = np.array([2,2,1,0,64])
    assert (chromwide_ibs(a,b,c,d) == expected).all()

    spa, spb = [SparseAlleles(x) for x in zip(*g1)]
    spc, spd = [SparseAlleles(x) for x in zip(*g2)]
    assert (chromwide_ibs(spa, spb, spc, spd) == expected).all()

    # Test assertions
    assert_raises(ValueError, chromwide_ibs, a, b, c, d, missingval=600)
    assert_raises(ValueError, chromwide_ibs, a, b, c, d, missingval=-1)
    
