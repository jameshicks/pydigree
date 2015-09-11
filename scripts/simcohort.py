#!/usr/bin/env python

import pydigree
import sys
import time
from pydigree import table
from pydigree.simulation import Architecture

nchrom = 1
popsize = 500
ngen = 100

pop = pydigree.Population(popsize)
for x in range(nchrom):
    c = pydigree.ChromosomeTemplate()
    c.add_genotype(0.1,0)
    pop.add_chromosome(c)

trait = Architecture('q','quantitative')
for i in range(nchrom):
    trait.add_effect_liability((i,0), 0, 1)
print trait

pop.initialize_pool()
for x in range(popsize):
    pop.add_founder_individual()

def gen(popsize):
    pop.random_mating_generation(popsize)
    pop.get_genotypes()
    pop.remove_ancestry()
    return pop.allele_frequency((0,0),1)

print [round(gen(popsize),2) for x in range(ngen)]


