

from itertools import chain

from nose import with_setup
from nose.tools import raises

from pydigree import Population, Individual
from pydigree.genotypes import Alleles
from pydigree.common import grouper

def setup_80freq():
    pop = Population()
    for i in range(600):
        ind = pop.founder_individual()
        ind.genotypes = [ [ [1],[1] ] ]
    for i in range(400):
        ind = pop.founder_individual()
        ind.genotypes = [ [ [2],[1] ] ]
    return pop

def test_allele_list():
    pop = Population()
    loc = 0,0
    # 2000 A alleles
    for x in range(1000):
        ind = pop.founder_individual()
        ind.genotypes = [ (Alleles([1]), Alleles([1]) )]
    # 500 A alleles, 500 B
    for x in range(500):
        ind = pop.founder_individual()
        ind.genotypes = [ (Alleles([1]), Alleles([2]) )]
    assert sorted(pop.allele_list(loc)) == [1]*2500 + [2]*500

def test_major_allele():
    pop = Population()

    # 2000 A alleles
    for x in range(1000):
        ind = pop.founder_individual()
        ind.genotypes = [ (Alleles([1]), Alleles([1]) )]
    # 500 A alleles, 500 B
    for x in range(500):
        ind = pop.founder_individual()
        ind.genotypes = [ (Alleles([1]), Alleles([2]) )]
    assert pop.major_allele((0,0)) == 1


def test_freq():
    pop = setup_80freq()
    loc = 0,0
    assert len(pop.allele_list(loc)) == 2 * 1000
    assert pop.allele_frequency(loc, 1) == ((2*600) + (1*400)) / float(2*1000)
    assert pop.allele_frequency(loc, 2) == (1*400) / float(2*1000)

