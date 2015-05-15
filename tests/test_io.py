from itertools import chain

from pydigree.io.base import genotypes_from_sequential_alleles
from pydigree.genotypes import GenotypedChromosome, SparseGenotypedChromosome, ChromosomeTemplate

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
    assert all(type(x) is GenotypedChromosome for x in chain.from_iterable(gts))
    assert all(type(x) is SparseGenotypedChromosome for x in chain.from_iterable(spgts))

    # Test to make sure the values are correct
    assert (gts[0][0] == ['1','1']).all()
    assert (gts[0][1] == ['2','1']).all()
    assert (gts[1][0] == ['2','2']).all()
    assert (gts[1][1] == ['2','1']).all()
