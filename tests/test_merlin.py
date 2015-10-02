import os
from itertools import izip

from nose.tools import raises, assert_raises
import numpy as np

from pydigree.exceptions import FileFormatError
from pydigree.io.merlin import phenotype_indices, genotype_indices
from pydigree.io.merlin import read_map, read_dat

testdir = os.path.dirname(os.path.abspath(__file__))
TESTDATA_DIR = os.path.join(testdir, 'test_data', 'merlin')

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

def test_read_map():
    mapfile = os.path.join(TESTDATA_DIR, 'merlin_sexaveraged.map')
    chroms = read_map(mapfile)

    expected_maps = [[1.0, 20.5, 100.1], [2.5, 21.1], [22.2]]
    expected_labs = [['rsA', 'rsB', 'rsC'], ['rsD', 'rsE'], ['rsF']]
    for chrom, expected_map, expected_lab in izip(chroms, expected_maps, expected_labs):
        assert all(np.array(chrom.labels) == np.array(expected_lab))
        assert all(np.array(chrom.genetic_map) == np.array(expected_map))

def test_read_dat():
    datfile = os.path.join(TESTDATA_DIR, 'testmerlin.dat')
    data = read_dat(datfile)

    expected_kinds = 'ASCCMMMMMM'
    expected_labs = ['affected', 'ignored', 'bmi', 'height']
    expected_labs += ['rsA', 'rsB', 'rsC', 'rsD', 'rsE', 'rsF']
    observed_kinds = [kind for kind, label in data]
    observed_labs = [label for kind, label in data]
    assert all(o == e for o,e in zip(observed_kinds, expected_kinds))
    assert all(o == e for o,e in zip(observed_labs, expected_labs))
