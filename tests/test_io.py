from nose.tools import raises, assert_raises
from itertools import chain
import os
import numpy as np
from pydigree.exceptions import FileFormatError

from pydigree.io.base import genotypes_from_sequential_alleles
from pydigree.io.vcf import vcf_allele_parser
from pydigree.io.merlin import phenotype_indices, genotype_indices
from pydigree.io import read_plink, read_vcf
from pydigree.genotypes import Alleles, SparseAlleles, ChromosomeTemplate

TESTDATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data')

def blank_chromosome(size=2):
    ch = ChromosomeTemplate()
    for i in range(size):
        ch.add_genotype()
    return ch

def test_seqalleles():
    chroms = [blank_chromosome(2) for x in range(2)]
    seqalleles = '1 2 1 1 2 2 2 1'.split()
    gts = genotypes_from_sequential_alleles(chroms, seqalleles)
    spgts = genotypes_from_sequential_alleles(chroms, seqalleles, sparse=True)
    
    # Test to make sure the types returned are correct
    assert all(type(x) is Alleles for x in chain.from_iterable(gts))
    assert all(type(x) is SparseAlleles for x in chain.from_iterable(spgts))

    # Test to make sure the values are correct
    assert (gts[0][0] == ['1','1']).all()
    assert (gts[0][1] == ['2','1']).all()
    assert (gts[1][0] == ['2','2']).all()
    assert (gts[1][1] == ['2','1']).all()

@raises(ValueError)
def test_seqalleles_raiseforbadmissingval():
    chroms = [blank_chromosome(2) for x in range(2)]
    seqalleles = '1 2 1 1 2 2 2 1'.split()
    gts = genotypes_from_sequential_alleles(chroms, seqalleles, missing_code=0)


def test_plink():
    plinkped = os.path.join(TESTDATA_DIR, 'plink', 'plink_test.ped')
    plinkmap = os.path.join(TESTDATA_DIR, 'plink', 'plink_test.map')
    
    peds = read_plink(pedfile=plinkped, mapfile=plinkmap)
    assert len(peds.chromosomes) == 2
    assert [x.nmark() for x in peds.chromosomes] == [2, 2]
    assert len(peds.individuals) == 2
    
    # Test individual 1 genotypes
    assert (peds['1','1'].genotypes[0][0] == ['1', '1']).all()
    assert (peds['1','1'].genotypes[0][1] == ['2', '1']).all()
    assert (peds['1','1'].genotypes[1][0] == ['2', '']).all()
    assert (peds['1','1'].genotypes[1][1] == ['2', '']).all()

    # Test individual 2 genotypes
    assert (peds['1', '2'].genotypes[0][0] == ['2', '2']).all()
    assert (peds['1', '2'].genotypes[0][1] == ['2', '2']).all()
    assert (peds['1', '2'].genotypes[1][0] == ['1', '']).all()
    assert (peds['1', '2'].genotypes[1][1] == ['2', '']).all()

    assert (peds['1','1'].genotypes[0][0].missing == [False, False]).all()
    assert (peds['1','1'].genotypes[1][0].missing == [False, True]).all()

def test_smartopen():
    from pydigree.io.smartopen import smartopen 
    datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data', 'compression')

    # Plain text 
    with smartopen(os.path.join(datadir, 'test')) as f:
        d = f.readlines()
        assert all(type(x) is str for x in d)
        assert [x.strip() for x in d] == ['genetics', 'pydigree', 'dna']

    # Gzip  
    with smartopen(os.path.join(datadir, 'test.gz')) as f:
        d = f.readlines()
        assert all(type(x) is str for x in d)
        assert [x.strip() for x in d] == ['genetics', 'pydigree', 'dna']

    # bz2
    with smartopen(os.path.join(datadir, 'test.gz')) as f:
        d = f.readlines()
        assert all(type(x) is str for x in d)
        assert [x.strip() for x in d] == ['genetics', 'pydigree', 'dna']
