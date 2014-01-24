#!/usr/bin/env python

import pydigree
import sys
import random
import time
from pydigree.simulation import Architecture

prefix = 'simulation'

replications = 500
genedrop_attempts = int(10e5)
causal_marker = 75


# Dominant model
trait = Architecture('affected', type='dichotomous')
trait.add_effect((0, causal_marker), {(2,2): 1, (2,1): 1, (1,1): 0})
trait.set_liability_threshold(.5)

print 'Reading pedigree'
peds = pydigree.read_ped(sys.argv[1])

print 'Generating chromosomes'
peds.add_chromosome(pydigree.read_gs_chromosome_template(sys.argv[2]))
for x in peds.chromosomes():
    print x
peds.chromosomes()[0].set_frequency(causal_marker, 0.0)

print 'Initializing chromosome pool'
for ped in peds:
    ped.initialize_pool(size=200)

for replication in xrange(replications):
    print 'Replication attempt %s' % (replication + 1)

    peds.clear_genotypes()
    for i in peds.individuals():
        if i.is_founder():
            i.get_genotypes(linkeq=True)

    for ped in peds:
        ped['1'].set_genotype((0, causal_marker), (2,1))

    for ped in peds:
        for attempt in xrange(genedrop_attempts):

            for ind in ped:
                if not ind.is_founder():
                    ind.clear_genotypes()
                ind.get_genotypes()

            calls = [(ind.predicted_phenotype(trait),ind.phenotypes['affected'])
                     for ind in ped if ind.phenotypes['affected'] is not None]

            accuracy = sum(x==y for x,y in calls) / float(len(calls))
            if accuracy > .9:
                print 'Complete after %s attempts' % attempt
        else:
            # TODO: Write own exception, RuntimeError isn't what we should use
            raise RuntimeError('Ran out of replication attempts')
    pydigree.write_ped(peds, '%s-replicate%s.ped' % (prefix, replication))
