#!/usr/bin/env python

import pydigree
import sys
import time
import itertools

def sharepairs(ped, inds, loc):
    from pydigree import count,ibs
    npairs = float( len(inds) * (len(inds)-1)/2 ) 
    r = [ibs(j.get_genotype(loc), k.get_genotype(loc)) for j,k
         in itertools.combinations(inds,2)]
    return count(True,r) / npairs


pop = pydigree.Population(5000)
                
ped = pydigree.io.read_ped(sys.argv[1], pop)[sys.argv[2]]
# Clear the genotypes, if present
ped.clear_genotypes()

c = pydigree.Chromosome()
c.add_genotype(0.5, 0)
print c
ped.add_chromosome(c)


affs = [x for x in ped if x.phenotypes['affected']]


niter = 1000
print "%s simulations" % niter

sim_share = []

t = time.time()
for x in range(niter):
    if x % 250 == 0:
        print 'Simulation %s' % x
    ped.simulate_ibd_states()
    s = sharepairs(ped, affs, (0,0))
    ped.clear_genotypes()
    sim_share.append(s)
t2 = time.time()-t

print "Maximum simulated allele sharing: %s" % max(sim_share)
print "Empiric P: %s" % (len([x for x in sim_share if x >= 1.0])/float(niter))
print "Time: %s (time per pedigree: %s)" % (t2,t2/niter)


ofilename = 'null.dist'
print "Outputting distribution to %s" % ofilename
with open(ofilename,'w') as of:
    for s in sim_share:
        of.write(str(s))
        of.write('\n')
