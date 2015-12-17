from nose.tools import raises, assert_raises

from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.genotypes import Alleles, SparseAlleles, ChromosomeTemplate
from pydigree.genotypes import LabelledAlleles, InheritanceSpan
from pydigree.exceptions import NotMeaningfulError
import numpy as np

#############
# Alleles tests
#############


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
    expected_value = np.array(
        [0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  0.,  0.])
    assert all(z == expected_value)

    # Test empty_like
    a = Alleles(np.zeros(10))
    b = a.empty_like()
    expected_value = Alleles(np.zeros(10), dtype=a.dtype)
    assert all(b == expected_value)

#############
# SparseAlleles tests
#############


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


def test_sparsealleles_meaninglesscomparisions():
    # Comparsions like >, <, >=, <= aren't meaningful for genotypes
    a = SparseAlleles(['1', '2', '3', ''])
    b = SparseAlleles(['1', '3', '2', ''])

    # Can't test an expression so we make a throwaway function to test
    assert_raises(NotMeaningfulError, lambda x, y: x < y, a, b)
    assert_raises(NotMeaningfulError, lambda x, y: x > y, a, b)
    assert_raises(NotMeaningfulError, lambda x, y: x >= y, a, b)
    assert_raises(NotMeaningfulError, lambda x, y: x <= y, a, b)


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

def test_sparseeq():
    a = SparseAlleles([1,2,3,4])
    b = SparseAlleles([1,3,3,4])

    obs = (a == b)
    expected = np.array([True, False, True, True]) 
    assert np.all(obs == expected)

def test_sparsealleles_emptylike():
    a = SparseAlleles([1,2,3,4])
    e = a.empty_like()
    assert e.non_refalleles.container == []

def test_sparsealleles_copyspan():
    a = SparseAlleles(np.array([0,0,0,0,0,0,0], dtype=np.int)+1, refcode=1)
    b = SparseAlleles(np.array([1,1,1,1,1,1,1], dtype=np.int)+1, refcode=1)

    a.copy_span(b, 2, 6)
    assert all(a.todense() == np.array([0,0,1,1,1,1,0]) + 1)

def test_array2missing():
    missingcode = 0
    vals = np.array([0, 1, 0, 0, 1, 0, 2], dtype=np.uint)
    assert all(SparseAlleles._array2missing(vals, missingcode) == np.array([0, 2, 3, 5]))


def test_array2nonref():
    from pydigree.datastructures import SortedPairContainer
    refcode = 0
    missingcode = 0 
    vals = np.array([0, 1, 0, 0, 1, 0, 2], dtype=np.uint)
    o = SparseAlleles._array2nonref(vals, refcode, missingcode)
    assert type(o) is SortedPairContainer
    assert o.container == [(1, 1), (4, 1), (6, 2)]

#############
# InheritanceSpan
#############


def test_inhertancespan():
    # Inheritance span doesnt really do much. In fact, the only
    # reason it has an __eq__ method is for unittesting some LabelledAlleles
    # methods.
    IS = InheritanceSpan
    a = IS(1, 0, 0, 0, 50)
    b = IS(1, 0, 0, 0, 50)
    c = IS(2, 1, 1, 30, 40)
    assert a == b
    assert a != c and b != c
#############
# LabelledAlleles
#############


def test_labelledalleles():
    IS = InheritanceSpan

    ngenos = 50
    p = Population()
    c = ChromosomeTemplate()
    for i in xrange(ngenos):
        c.add_genotype()
    p.add_chromosome(c)

    a = Individual(p, 1)
    actual = LabelledAlleles.founder_chromosome(a, 0, 0, chromobj=c)
    expected = LabelledAlleles(spans=[IS(a, 0, 0, 0, ngenos)], chromobj=c)
    assert actual == expected

def test_labelledallele_delabeler():
    ngenos = 10  # Number of genotypes per chromosome
    if ngenos % 2 == 1:
        raise ValueError('Even number of genotypes needed')

    p = Population()
    c = ChromosomeTemplate()
    for i in xrange(ngenos):
        c.add_genotype()
    p.add_chromosome(c)

    a = Individual(p, 1)
    a._init_genotypes(blankchroms=False)
    a.genotypes[0][0] = Alleles([1]*ngenos)
    a.genotypes[0][1] = Alleles([2]*ngenos)

    b = Individual(p, 2)
    b._init_genotypes(blankchroms=False)
    b.genotypes[0][0] = Alleles([3] * ngenos)
    b.genotypes[0][1] = Alleles([4] * ngenos)

    chromatid_spans = [InheritanceSpan(a, 0, 0, 0, ngenos/2),
                       InheritanceSpan(b, 0, 1, ngenos/2, ngenos)]
    chromatid = LabelledAlleles(spans=chromatid_spans, chromobj=c)

    expected_value = [1]*(ngenos/2) + [4] * (ngenos/2)
    expected_value = Alleles(expected_value)

    actual_value = chromatid.delabel()
    assert all(actual_value == expected_value)

