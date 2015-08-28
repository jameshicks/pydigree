from nose.tools import raises, assert_raises

from pydigree.genotypes import Alleles, SparseAlleles
from pydigree.exceptions import NotMeaningfulError
import numpy as np

def test_alleles():
    a = Alleles(['1', '2', '3', ''])
    b = Alleles(['1', '3', '2', ''])

    assert a.nmark() == b.nmark() == 4
    
    # Test missingness
    assert a.missingcode == ''
    assert (a.missing == np.array([False, False, False, True])).all()
    assert (a.missing == b.missing).all()
    
    eq = (a == b)
    assert (eq == np.array([True, False, False, True])).all()
    
    # Test copy span
    z = Alleles(np.zeros(10))
    o = Alleles(np.ones(10))
    
    z.copy_span(o, 5, 8)
    expected_value = np.array([ 0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.])
    assert all(z == expected_value)

    # Test empty_like
    a = Alleles(np.zeros(10))
    b = a.empty_like()
    expected_value = Alleles(np.zeros(10), dtype=a.dtype)
    assert all(b == expected_value)

def test_sparsealleles():
    a = SparseAlleles(['1', '2', '3', ''])
    b = SparseAlleles(['1', '3', '2', ''])

    assert a.nmark() == b.nmark() == 4
    assert (a.missing == np.array([False, False, False, True])).all()
    assert (a.missing == b.missing).all()

    # Test todense
    assert (a.todense() == ['1', '2', '3', '']).all()
    assert isinstance(a.todense(), Alleles)


    # Test equality between chromosome equality
    eq = (a == b)
    assert (eq == np.array([True, False, False, True])).all()
    assert ((a == b.todense()) == np.array([True, False, False, True])).all()
    assert np.all((a != b) == np.logical_not(eq)) 

def test_alleles_meaninglesscomparisions():
    # Comparsions like >, <, >=, <= aren't meaningful for genotypes
    a = Alleles(['1', '2', '3', ''])
    b = Alleles(['1', '3', '2', ''])
    
    # Can't test an expression so we make a throwaway function to test
    assert_raises(NotMeaningfulError, lambda x,y: x < y, a, b)
    assert_raises(NotMeaningfulError, lambda x,y: x > y, a, b)
    assert_raises(NotMeaningfulError, lambda x,y: x >= y, a, b)
    assert_raises(NotMeaningfulError, lambda x,y: x <= y, a, b)

def test_sparsealleles_meaninglesscomparisions():
    # Comparsions like >, <, >=, <= aren't meaningful for genotypes
    a = SparseAlleles(['1', '2', '3', ''])
    b = SparseAlleles(['1', '3', '2', ''])
    
    # Can't test an expression so we make a throwaway function to test
    assert_raises(NotMeaningfulError, lambda x,y: x < y, a, b)
    assert_raises(NotMeaningfulError, lambda x,y: x > y, a, b)
    assert_raises(NotMeaningfulError, lambda x,y: x >= y, a, b)
    assert_raises(NotMeaningfulError, lambda x,y: x <= y, a, b)


@raises(ValueError)
def test_sparse_wrongtypecomparsion():
    a = SparseAlleles(['1', '2', '3', ''])
    a == 3

@raises(ValueError)
def test_sparse_wrongsizecomparision():
    a = SparseAlleles(['1', '2', '3', ''])
    b = SparseAlleles(['1', '3'])
    a == b

@raises(ValueError)
def test_sparse_norefscalarcomparison():
    a = SparseAlleles(['1', '2', '3', ''])
    a == '3'

def test_array2missing():
    missingcode = 0
    vals = np.array([0,1,0,0,1,0,2], dtype=np.uint)
    assert SparseAlleles._array2missing(vals, missingcode) == [0, 2, 3, 5]

def test_array2nonref():
    refcode = 0 
    vals = np.array([0,1,0,0,1,0,2], dtype=np.uint)
    o = SparseAlleles._array2nonref(vals, refcode)
    assert o == {1:1, 4:1, 6:2}
