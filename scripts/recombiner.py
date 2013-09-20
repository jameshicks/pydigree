#!/usr/bin/env python

import random
import pydigree

nmark = 5000
poolsize = 1000

c = pydigree.Chromosome()
for x in range(nmark): c.add_genotype(random.random(),.1)
mapp = [.1]*nmark

d =0

print "Initializing pool"
pool = [c.random_chromosome() for x in range(poolsize)]

for x in range(10):
    print "Iteration %d" % d
    recombs = [pydigree.recombine(random.choice(pool),random.choice(pool),mapp)]
    pool = random.sample(recombs+pool,poolsize)
    d += 1
