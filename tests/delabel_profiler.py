
try:
    import line_profiler
except ImportError:
    print("No line profiler, skipping test.")
    import sys
    sys.exit(0)

nrep = 1000
ngenos = 100000 # Number of genotypes per chromosome


from pydigree import Individual, Population
from pydigree.genotypes import ChromosomeTemplate, Alleles
from pydigree.genotypes import chromatid_delabeler
from pydigree.genotypes import AncestralAllele
from pydigree.common import spans

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

chromatid = [AncestralAllele(a,0)] * (100/2) + [AncestralAllele(b,1)] * (100/2)
chromatid = Alleles(chromatid * (ngenos/100))


#### Test
func = chromatid_delabeler
test_values = chromatid
profile = line_profiler.LineProfiler(func)
for i in xrange(nrep):
    profile.runcall(func, test_values, 0)
profile.print_stats()


##########
from pydigree.genotypes import InheritanceSpan, LabelledAlleles

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

chromatid_spans = []
spansize = ngenos/1000/2
for start in range(0,ngenos,spansize*2):
    aspan = InheritanceSpan(a, 0, 0, start, start+spansize)
    bspan = InheritanceSpan(b, 0, 1, start+spansize, start+(spansize*2))
    chromatid_spans.extend([aspan, bspan])

lchromatid = LabelledAlleles(spans=chromatid_spans, chromobj=c)


assert all(chromatid_delabeler(chromatid, 0) == lchromatid.delabel())
#### Test
func = lchromatid.delabel
test_values = lchromatid
profile = line_profiler.LineProfiler(func)
for i in xrange(nrep):
    profile.runcall(func)
profile.print_stats()

