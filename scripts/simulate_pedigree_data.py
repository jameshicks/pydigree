#!/usr/bin/env python

import argparse

import pydigree
from pydigree.io import read_ped, read_gs_chromosome_template
from pydigree.simulation import ConstrainedMendelianSimulation
from pydigree.simulation import NaiveGeneDroppingSimulation


parser = argparse.ArgumentParser()
parser.add_argument('--template', dest='template', metavar='template_pedigrees',
                    required=True, help='PLINK/LINkAGE/MERLIN format pedigree file')
parser.add_argument('--chromosomes', dest='chromosomes',nargs='+', 
                    help='genomeSIMLA format chromosome templates')
parser.add_argument('--method',  required=True, help='Simulation method')
parser.add_argument('--replications', type=int, metavar='nrep',  
                    default=1000, help='Number of replicates')
parser.add_argument('--constraintfile', metavar='file',  help='Constriant file')
parser.add_argument('--smod', metavar='file',  help='smod file')
parser.add_argument('--prefix', type=str, help='prefix for output file',
                    default='simulation')
args = parser.parse_args()

# Read the pedigrees
template = read_ped(args.template)

# Read the chromosomes
for cfile in args.chromosomes:
    c = read_gs_chromosome_template(cfile)
    template.add_chromosome(c)

if args.method.lower() == 'constrained':
    sim = ConstrainedMendelianSimulation(template, replications=args.replications)
elif args.method.lower() == 'genedrop':
    sim = NaiveGeneDroppingSimulation(template, replications=args.replications)
else:
    print 'Unknown method: %s' % args.method

sim.prefix = args.prefix

sim.run(verbose=True)
