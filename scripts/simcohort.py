#!/usr/bin/env python3

import pydigree
from pydigree.simulation import QuantitativeTrait

nchrom = 1
popsize = 500
ngen = 100

pop = pydigree.Population(popsize)
for x in range(nchrom):
    c = pydigree.ChromosomeTemplate()
    c.add_genotype(0.1, 0)
    pop.add_chromosome(c)

trait = QuantitativeTrait('q', 'quantitative')
for i in range(nchrom):
    trait.add_effect_liability((i, 0), 0, 1)
print(trait)

pop.initialize_pool()
for x in range(popsize):
    pop.add_founder_individual()


def gen(ninds):
    "Advance by one generation, with ninds in the new generation"
    pop.advance_generation(ninds)
    pop.get_genotypes()
    pop.remove_ancestry()
    return pop.allele_frequency((0, 0), 1)

print([round(gen(popsize), 2) for x in range(ngen)])
