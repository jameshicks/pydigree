import numpy as np

import pydigree as pyd 
from pydigree.simulation import Architecture

ninds = 5000
nloc = 1000
maf = 0.5

traitname = 'synthetic'
# Create population
pop = pyd.Population()

# Create chromosomes
for i in xrange(100):
	c = pyd.ChromosomeTemplate()
	c.add_genotype(maf, 0)
	pop.add_chromosome(c)

# Create trait architecture
trait = Architecture('synthetic', 'quantitative', chromosomes=pop.chromosomes)
for i in xrange(100):
	trait.add_effect((i,0), 1 * (-1 if i % 2 else 1))

print 'Locus mean genotypic value: {}'.format(trait.effects[0].expected_genotypic_value)
print 'Locus variance: {}'.format(trait.effects[0].locus_additive_variance)

print 'Expected trait mean: {}'.format(trait.expected_genotypic_value)
print 'Expected trait variance: {}'.format(trait.additive_genetic_variance)
for i in xrange(ninds):
	i = pop.founder_individual()
	i.get_genotypes(linkeq=True)
	i.phenotypes[traitname] = trait.predict_phenotype(i)

y = np.array([i.phenotypes[traitname] for i in pop.individuals])
print 'Observed trait mean {}'.format(y.mean())
print 'Observed trait variance: {}'.format(y.var())

