from pydigree import Individual, Population
from pydigree.genotypes import ChromosomeTemplate, Alleles
from pydigree.genotypes import chromatid_delabeler, AncestralAllele

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
    p = Population()
    c = ChromosomeTemplate()
    for i in xrange(10):
        c.add_genotype()
    p.add_chromosome(c)

    a = Individual(p, 1)
    a._init_genotypes(blankchroms=False)
    a.genotypes[0][0] = Alleles([1]*10)
    a.genotypes[0][1] = Alleles([2]*10)

    b = Individual(p, 2)
    b._init_genotypes(blankchroms=False)
    b.genotypes[0][0] = Alleles([3] * 10)
    b.genotypes[0][1] = Alleles([4] * 10)

    chromatid = [AncestralAllele(a,0)] * 5 + [AncestralAllele(b,1)] * 5
    chromatid = Alleles(chromatid)
    expected_result = Alleles([1,1,1,1,1,4,4,4,4,4])
    result = chromatid_delabeler(chromatid, 0)
    assert all(result == expected_result)
    assert isinstance(result, Alleles)