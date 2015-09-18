from nose.tools import raises, assert_raises
import numpy as np

from pydigree.exceptions import FileFormatError
from pydigree.io.merlin import phenotype_indices, genotype_indices

def test_merlin_phenotype_indices():
    assert all(phenotype_indices(['A', 'C', 'T']) == np.array([1,1,1]))
    assert all(phenotype_indices(['A', 'S', 'A']) == np.array([1,0,1]))
    assert all(phenotype_indices(['M', 'T', 'S2']) == np.array([0,0,1,0,0]))
    assert_raises(FileFormatError, phenotype_indices, ['S', 'W', 'O'])

def test_merlin_genotype_indices():
    assert all(genotype_indices(['A', 'C', 'T']) == np.array([0,0,0,]))
    assert all(genotype_indices(['A', 'S2', 'M']) == np.array([0,0,0,1,1]))
    assert all(genotype_indices(['M', 'T', 'S2']) == np.array([1,1,0,0,0]))
    assert_raises(FileFormatError, genotype_indices, ['S', 'W', 'O'])