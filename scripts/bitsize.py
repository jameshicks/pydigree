#!/usr/bin/env python

import pydigree
import sys

# Prints the bit size of each pedigree.
# The bitsize is defined as 2*n-f where n is the number of nonfounders
# and f is the number of founders. This represents the number of bits
# it takes to represent the inheritance vector in the Lander-Green
# algorithm.

ped = pydigree.io.read_ped(sys.argv[1])

for pedigree in sorted(ped,key=lambda x: x.label):
    print pedigree.label, pedigree.bit_size()
