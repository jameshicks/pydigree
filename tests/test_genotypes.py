from nose.tools import raises

from pydigree.genotypes import GenotypedChromosome, SparseGenotypedChromosome
import numpy as np

def test_genotypedchromosome():
    a = GenotypedChromosome(['1', '2', '3', ''])
    b = GenotypedChromosome(['1', '3', '2', ''])

    assert a.nmark() == b.nmark() == 4
    
    # Test missingness
    assert a.missingcode == ''
    assert (a.missing == np.array([False, False, False, True])).all()
    assert (a.missing == b.missing).all()
    
    eq = (a == b)
    assert (eq == np.array([True, False, False, True])).all()
    
def test_sparsegenotypedchromosome():
    a = SparseGenotypedChromosome(['1', '2', '3', ''])
    b = SparseGenotypedChromosome(['1', '3', '2', ''])

    assert a.nmark() == b.nmark() == 4
    assert (a.missing == np.array([False, False, False, True])).all()
    assert (a.missing == b.missing).all()

    # Test todense
    assert (a.todense() == ['1', '2', '3', '']).all()
    assert isinstance(a.todense(), GenotypedChromosome)


    # Test equality between chromosome equality
    eq = (a == b)
    assert (eq == np.array([True, False, False, True])).all()
    assert ((a == b.todense()) == np.array([True, False, False, True])).all()

@raises(ValueError)
def test_sparsegc_wrongtypecomparsion():
    a = SparseGenotypedChromosome(['1', '2', '3', ''])
    a == 3

@raises(ValueError)
def test_sparsegc_norefscalarcomparison():
    a = SparseGenotypedChromosome(['1', '2', '3', ''])
    a == '3'
