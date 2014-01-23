#!/usr/bin/env python

import pydigree
import sys
import time
import itertools

def sharepairs(ped, inds, loc):
    from pydigree import count,ibs
    npairs = float( len(inds) * (len(inds)-1)/2 ) 
    r = []
    for j,k in itertools.combinations(inds,2):
        g1 = j.get_genotype(loc)
        g2 = k.get_genotype(loc)
        r.append(ibs(g1,g2) > 0)
    return count(True,r)/npairs

nmark = 1
chrom_length_cm = 0
intermark_dist_cm = (float(chrom_length_cm)/nmark)/100

pop = pydigree.Population(5000)
for rx in range(1):
    c = pydigree.Chromosome()
    for x in range(nmark): c.add_genotype(frequency,intermark_dist_cm)
    print c
    pop.add_chromosome(c)
                
ped = pydigree.read_ped(sys.argv[1], pop)
affs = [x for x in ped if x.phenotypes['aff']]


niter= 10**4
print "%s simulations" % niter

sim_share = []

t = time.time()
for x in range(niter):
    ped.simulate_ibd_states()
    s = sharepairs(ped,affs,(0,0))
    ped.clear_genotypes()
    sim_share.append(s)
t2 = time.time()-t
print "Maximum simulated allele sharing: %s" % max(sim_share)
print "Empiric P: %s" % (len([x for x in sim_share if x >= 1.0])/float(niter))
print "Time: %s (time per pedigree: %s)" % (t2,t2/niter)
