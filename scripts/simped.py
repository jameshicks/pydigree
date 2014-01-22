#!/usr/bin/env python

import pydigree
import sys
import random
import time

nmark = 50000
chrom_length_cm = 150
intermark_dist_cm = (float(chrom_length_cm)/nmark)/100

pop = pydigree.Population(5000)
for rx in range(1):
    c = pydigree.Chromosome()
    for x in range(nmark): c.add_genotype(random.random(),intermark_dist_cm)
    print c
    pop.add_chromosome(c)
                
ped = pydigree.read_ped(sys.argv[1], pop)

