#!/usr/bin/env python

import pydigree
import itertools
import sys

pop = pydigree.Population(5000)
ped = pydigree.read_ped(sys.argv[1], pop)
ids = sorted([x.id for x in ped])


# itertools combinations with replacement isn't available in 2.6
# So we'll do this as a two step process
# Pairwise kinship coefs
for x,y in itertools.combinations(ids,2):
    k = pydigree.kinship(ped[x],ped[y])
    print 1,x,y,k
# Inbreeding coefs.
for x in ids:
    print 1,x,x,ped[x].inbreeding()
