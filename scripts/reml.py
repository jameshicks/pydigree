#!/usr/bin/env python

import sys
import argparse

import pydigree
from pydigree.mixedmodel import MixedModel

parser = argparse.ArgumentParser()
parser.add_argument('--ped', required=True,
                    help='Pedigree file', dest='pedf')
parser.add_argument('--phen', required=True,
                    help='Phenotype CSV file', dest='phenf')
parser.add_argument('--outcome', required=True, help='Outcome phenotype')
parser.add_argument('--fixedeffects', '--fixefs', nargs="*", dest='fixefs',
                    help='Fixed effects for model')
parser.add_argument('--maxmethod', help='Method for maximization',
                     default='Fisher')
args = parser.parse_args()


print 'Reading files'
peds = pydigree.io.read_ped(args.pedf)
pydigree.io.read_phenotypes(peds, args.phenf)

m = MixedModel(peds, outcome=args.outcome, fixed_effects=args.fixefs)
print 'Calculating Kinships'
m.add_genetic_effect()
print 'Done'
m.fit_model()


m.maximize(method=args.maxmethod, verbose=True)

m.summary()
