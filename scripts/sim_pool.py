#!/usr/bin/env python3

import argparse
import pydigree as pyd
from pydigree.simulation.chromosomepool import ChromosomePool
from pydigree.io.genomesimla import read_gs_chromosome_template
from pydigree.population import logistic_growth

parser = argparse.ArgumentParser()
parser.add_argument('--chromosomes', dest='chromosomes',nargs='+', 
                    help='genomeSIMLA format chromosome templates')
parser.add_argument('--n0', dest='n0', type=int, help='Initial pool size', default=5000)
parser.add_argument('--rate', dest='rate', type=float, help='Growth rate', default=1.0)
parser.add_argument('--final', type=int, help='Final pool size', default=50000)
parser.add_argument('--gens', type=int, default=20)
parser.add_argument('--initial', help='Prefix for initial data for pool (plink format)')
args = parser.parse_args()

if not (args.chromosomes or args.initial):
	print('One of --chromosomes or --initial required')
	exit(1)

if args.initial:
	pop = pyd.io.read_plink(prefix=args.initial)
	pool = ChromosomePool.from_population(pop)
else:
	chroms = [read_gs_chromosome_template(x) for x in args.chromosomes]
	pool = ChromosomePool(chromosomes=chroms, size=args.n0)
	print('Creating pool')
	pool.initialize_pool(args.n0)

gensize = lambda x: int(logistic_growth(pool.n0, args.rate, args.final, x))

for x in range(args.gens):
    print('Generation {}: {}'.format(x, gensize(x)))
    pool.iterate_pool(gensize(x))
