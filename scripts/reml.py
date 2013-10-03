#!/usr/bin/env python

import sys
import pydigree
from pydigree.mixedmodel import MixedModel

pedf,phenf = sys.argv[1:3]
outcome = sys.argv[3]
fixefs = sys.argv[4:]


print 'Reading files'
peds = pydigree.read_ped(pedf)
pydigree.read_phenotypes(peds,phenf)
print 'Done'
m = MixedModel(peds,outcome=outcome,fixed_effects=fixefs)
print 'Calculating Kinships'
m.add_genetic_effect()
print 'Done'
m.fit_model()
m.maximize()

m.summary()
