#!/usr/bin/env python

import numpy as np

from pydigree import ChromosomeTemplate
from pydigree.recombination import recombine
from pydigree.rand import choice as randchoice
from pydigree.rand import sample_with_replacement

nmark = 5000
poolsize = 10000
ngen = 25

c = ChromosomeTemplate()
for x in range(nmark):
    c.add_genotype(np.random.random(), .1)
mapp = [.1]*nmark



print "Initializing pool"
pool = [c.linkageequilibrium_chromosome(sparse=True) for x in range(poolsize)]

for x in xrange(ngen):
    print "Iteration {}".format(x)
    recombs = [recombine(randchoice(pool),
                         randchoice(pool),
                         mapp)
               for x in xrange(poolsize)]
    pool = sample_with_replacement(recombs+pool, poolsize)

