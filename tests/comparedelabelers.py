from pydigree import Individual, Population
from pydigree.genotypes import ChromosomeTemplate, Alleles
from pydigree.genotypes import  chromatid_delabeler, chromatid_delabeler3
from pydigree.genotypes import AncestralAllele

ngenos = 100000 # Number of genotypes per chromosome
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
expected_result = Alleles([1] * (ngenos/2) + [4]*(ngenos/2))
assert len(expected_result) == ngenos # Sanity check    

#import IPython; IPython.embed()

nrep = 10000
for i in xrange(nrep):
    chromatid_delabeler(chromatid,0)
    #chromatid_delabeler2(chromatid,0)
    chromatid_delabeler3(chromatid,0)