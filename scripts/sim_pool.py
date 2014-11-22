#!/usr/bin/env python

import argparse
from pydigree.simulation.chromosomepool import ChromosomePool
from pydigree.io.genomesimla import read_gs_chromosome_template
from pydigree.population import logistic_growth

parser = argparse.ArgumentParser()
parser.add_argument('--chromosomes', dest='chromosomes',nargs='+', 
                    help='genomeSIMLA format chromosome templates')
parser.add_argument('--n0', dest='n0', type=int, help='Initial pool size', default=5000)
parser.add_argument('--rate', dest='rate', type=float, help='Growth rate', default=1.0)
parser.add_argument('--final', type=int, help='Final pool size', default=50000)
parser.add_argument('--gens', type=int)
args = parser.parse_args()


gensize = lambda x: logistic_growth(args.n0, args.rate, args.final, x)

chroms = [read_gs_chromosome_template(x) for x in args.chromosomes]

pool = ChromosomePool(chromosomes=chroms, size=args.n0)
pool.initialize_pool(args.n0)

for x in xrange(args.gens):
    print 'Generation {}: {}'.format(x, gensize(x))
    pool.iterate_pool(gensize(x))
