from nose import with_setup
from pydigree import Population, Individual

def setup_80freq():
    pop = Population()
    for i in xrange(600):
        ind = pop.founder_individual()
        ind.genotypes = [ [ [1],[1] ] ]
    for i in xrange(400):
        ind = pop.founder_individual()
        ind.genotypes = [ [ [2],[1] ] ]
    return pop

def test_freq():
    pop = setup_80freq()
    loc = 0,0
    assert len(pop.allele_list(loc)) == 2 * 1000
    assert pop.allele_frequency(loc, 1) == ((2*600) + (1*400)) / float(2*1000)
    assert pop.allele_frequency(loc, 2) == (1*400) / float(2*1000)

def test_ld():
    pass
