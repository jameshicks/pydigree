#!/usr/bin/env python

import pydigree
import itertools
import sys

pop = pydigree.Population(5000)
ped = pydigree.read_ped(sys.argv[1], pop)

for pedigree in ped:
    lab = pedigree.label
    ids = sorted([x.id for x in pedigree])
    # itertools combinations with replacement isn't available in 2.6
    # So we'll do this as a two step process
    # Pairwise kinship coefs
    for x,y in itertools.combinations(ids,2):
        k = pedigree.kinship(x,y)
        print lab,x,y,k
    # Inbreeding coefs.
    for x in ids:
        print lab,x,x,pedigree.inbreeding(x)
