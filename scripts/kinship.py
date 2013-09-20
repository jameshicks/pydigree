#!/usr/bin/env python

import pydigree
import itertools
import sys

pop = pydigree.Population(5000)
ped = pydigree.read_ped(sys.argv[1], pop)
ids = sorted([x.id for x in ped])

import pdb;pdb.set_trace()
for x,y in itertools.combinations(ids,2):
    k = pydigree.kinship(ped[x],ped[y])
    print 1,x,y,k
for x in ids:
    print 1,x,x,ped[x].inbreeding()
