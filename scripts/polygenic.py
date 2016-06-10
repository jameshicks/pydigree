#!/usr/bin/env python

import sys
import argparse

import pydigree
import numpy as np
from pydigree.stats import MixedModel

parser = argparse.ArgumentParser()
parser.add_argument('--ped', required=True,
                    help='Pedigree file', dest='pedf')
parser.add_argument('--phen', required=True,
                    help='Phenotype CSV file', dest='phenf')
parser.add_argument('--outcome', required=True, help='Outcome phenotype')
parser.add_argument('--fixedeffects', '--fixefs', nargs="*", dest='fixefs',
                    help='Fixed effects for model')
parser.add_argument('--reml', action='store_true', 
    help='Use restricted maximum likelihood (REML) for maximization')
parser.add_argument('--maxmethod', help='Method for maximization',
                    default='Fisher')
parser.add_argument('--starts', help='Starting values for variance components',
                    nargs='*', default=None)
parser.add_argument('--interact', action='store_true',
                    help='Enter IPython shell after maximization')
parser.add_argument('--d7', '--dominance', action='store_true', dest='d7',
                    help='Include dominance term in model')
parser.add_argument('--garbley', action='store_true',
                    help='Replace y values with random normal deviates')
parser.add_argument('--inflate', action='store_true')
parser.add_argument(
    '--center', action='store_true', help='Center outcome data')
args = parser.parse_args()


print 'Reading files'
peds = pydigree.io.read_ped(args.pedf)
pydigree.io.read_phenotypes(peds, args.phenf)

m = MixedModel(peds, outcome=args.outcome, fixed_effects=args.fixefs)
print 'Calculating Kinships'
m.add_genetic_effect()

if args.d7:
    m.add_genetic_effect(kind='dominance')

print 'Done'
m.fit_model()

if args.center:
    m.y = m._centery()

if args.inflate:
    m.y *= 100
if args.garbley:
    m.y = np.matrix(np.random.normal(10, 5, len(m.y))).T

starts = args.starts
if starts is not None:
    starts = [float(x) for x in starts]
m.maximize(method=args.maxmethod, verbose=True, starts=starts, restricted=args.reml)

m.summary()

if args.interact:
    try:
        from IPython import embed
        embed()
    except ImportError:
        print "IPython not found!"
