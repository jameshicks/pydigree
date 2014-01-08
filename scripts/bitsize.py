#!/usr/bin/env python

import pydigree
import sys

ped = pydigree.read_ped(sys.argv[1])

for pedigree in sorted(ped,key=lambda x: x.label):
    print pedigree.label, pedigree.bit_size()
