#!/usr/bin/env python

import sys
import pydigree
from pydigree.mixedmodel import MixedModel

pedf, phenf = sys.argv[1:3]
outcome = sys.argv[3]
fixefs = sys.argv[4:]

print 'Reading files'
peds = pydigree.io.read_ped(pedf)
pydigree.io.read_phenotypes(peds, phenf)

m = MixedModel(peds, outcome=outcome, fixed_effects=fixefs)
print 'Calculating Kinships'
m.add_genetic_effect()
print 'Done'
m.fit_model()
m.maximize()

m.summary()
