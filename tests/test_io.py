from itertools import chain
import os

from pydigree.io.base import genotypes_from_sequential_alleles
from pydigree.io import read_plink, read_vcf
from pydigree.genotypes import Alleles, SparseAlleles, ChromosomeTemplate

def blank_chromosome(size=2):
    ch = ChromosomeTemplate()
    for i in xrange(size):
        ch.add_genotype()
    return ch

def test_seqalleles():
    chroms = [blank_chromosome(2) for x in xrange(2)]
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

def test_plink():
    plinkped = os.path.join('test_data', 'plink_test.ped')
    plinkmap = os.path.join('test_data', 'plink_test.map')
    
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

def test_vcf():
    testvcf = os.path.join('test_data', 'test.vcf')
    
    pop = read_vcf(testvcf, minqual=0)
    assert len(pop.individuals) == 3
    assert len(pop.chromosomes) == 2

    # Test individual with good genotypes
    assert (pop['NA00001'].genotypes[1][0].todense() == ['0', '0', '1', '0', '0', '0', '1']).all()
    assert (pop['NA00001'].genotypes[1][1].todense() == ['0', '0', '2', '0', '1', '1', '1']).all()
    assert not pop['NA00001'].genotypes[1][0].missing.all()

    # Test individual with bad genotypes.
    assert pop['NA00003'].genotypes[1][0].missing.all()

    # Test variant level quality filtering 
    pop = read_vcf(testvcf, minqual=20)
    assert pop.chromosomes[1].nmark() == 6

    # Test for FILTER == PASS
    pop = read_vcf(testvcf, minqual=0, require_pass=True)
    assert pop.chromosomes[1].nmark() == 6
