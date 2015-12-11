#!/usr/bin/env python
from __future__ import division

import sys
import pydigree
import time
import numpy as np


n=5000
k = 2 * 10**5
r = .5
gens = 15
nmark = 2000
chrom_length_cm = 200
intermark_dist_cm = (float(chrom_length_cm)/nmark)/100

pop = pydigree.Population(n)
for rx in range(1):
    c = pydigree.ChromosomeTemplate()
    for x in range(nmark): c.add_genotype(np.random.random(),intermark_dist_cm)
    print c
    pop.add_chromosome(c)
print 'Initializing pool'
pop.initialize_pool()
gen_sizes = [pydigree.logistic_growth(n,r,k,t) for t in range(gens)]
print 'Iterating pool'
for i,g in enumerate(gen_sizes):
    #print 'Generation %d (%d chromosomes)...' % (i,pop.chrom_pool_size())
    #ng = pydigree.logistic_growth(p0,r,pfinal,g)
    t = time.time()
    pop.iterate_pool(g)
    t2 = time.time()
    print '\t'.join(str(x) for x in [i, int(g), t2-t, (t2-t)/int(g)])
