#!/usr/bin/env python
import sys

from pydigree.simulation import ConstrainedMendelianSimulation
from pydigree import read_ped, read_gs_chromosome_template

fixedlocus = 0, 10

peds = read_ped(sys.argv[1])
peds.add_chromosome(read_gs_chromosome_template(sys.argv[2]))

sim = ConstrainedMendelianSimulation(peds)
sim.add_genotype_constraint(peds['1']['1'], fixedlocus, 4, 'P')
sim.add_ibd_constraint(peds['1']['3'], peds['1']['1'], fixedlocus, 'P')
for i in xrange(100):
    sim.replicate()
    print peds['1']['1'].get_genotype(fixedlocus), peds['1']['3'].get_genotype(fixedlocus)

print 'done'
