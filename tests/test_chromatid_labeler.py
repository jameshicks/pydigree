from pydigree import Individual, Population
from pydigree.genotypes import ChromosomeTemplate, Alleles
from pydigree.genotypes import chromatid_delabeler
from pydigree.genotypes import AncestralAllele

def test_labeler():
    p = Population()
    c = ChromosomeTemplate()
    for i in xrange(10):
        c.add_genotype()
    p.add_chromosome(c)

    ind = Individual(p, '1')
    ind.label_genotypes()
    hapA = AncestralAllele(ind, 0)
    hapB = AncestralAllele(ind, 1)

    assert all(ind.genotypes[0][0] == hapA)
    assert all(ind.genotypes[0][1] == hapB)

def test_delabeler():
    ngenos = 10 # Number of genotypes per chromosome
    if ngenos % 2 == 1: raise ValueError('Even number of genotypes needed')

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

    chromatid = [AncestralAllele(a,0)] * (ngenos/2) + [AncestralAllele(b,1)] * (ngenos/2)
    chromatid = Alleles(chromatid)
    expected_result = Alleles([1,1,1,1,1,4,4,4,4,4])
    assert len(expected_result) == ngenos # Sanity check    
    result = chromatid_delabeler(chromatid, 0)
    assert all(result == expected_result)
    assert isinstance(result, Alleles)