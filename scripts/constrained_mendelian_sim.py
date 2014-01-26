#!/usr/bin/env python
import sys

from pydigree.simulation import *
from pydigree import read_ped, read_gs_chromosome_template

fixedlocus = 0, 4
print 'Reading ped'
peds = read_ped(sys.argv[1])
print 'Reading chromosome'
peds.add_chromosome(read_gs_chromosome_template(sys.argv[2]))

sim = ConstrainedMendelianSimulation(peds)

trait = Architecture('affected', type='dichotomous')
trait.add_effect(fixedlocus, {(4,2): 1, (4,1): 1, (1,1): 0})
trait.set_liability_threshold(.5)
sim.set_trait(trait)

sim.add_genotype_constraint(peds['1']['1'], fixedlocus, 4, 'P')
sim.add_ibd_constraint(peds['1']['3'], peds['1']['1'], fixedlocus, 'P')
for i in xrange(100):
    print 'Replicate %s' % i
    sim.replicate()

