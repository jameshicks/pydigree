#!/usr/bin/env python3
"Simulates quantitative traits in a pedigree"

import argparse

import pydigree
from pydigree.io import read_ped
from pydigree.io.genomesimla import read_gs_chromosome_template
from pydigree.simulation.genedrop import ConstrainedMendelianSimulation
from pydigree.simulation.genedrop import NaiveGeneDroppingSimulation
from pydigree.simulation import QuantitativeTrait

parser = argparse.ArgumentParser()
parser.add_argument('--template', dest='template', 
                    required=True, metavar='template_pedigrees',
                    help='PLINK/LINKAGE/MERLIN format pedigree file')
parser.add_argument('--chromosomes', dest='chromosomes', nargs='+',
                    help='genomeSIMLA format chromosome templates')

parser.add_argument('--method',  required=True, help='Simulation method',
                    choices=('constrained', 'genedrop'))

parser.add_argument('--replications', type=int, metavar='nrep',
                    default=1000, help='Number of replicates')
parser.add_argument('--ibd', action='store_true',
                    help='Write IBD states to file before genotypes')
parser.add_argument('--constraintfile', metavar='file', help='Constriant file')

parser.add_argument('--effectfile', dest='effectfile',
                    help='File containing marker effects for trait simulation')
parser.add_argument('--effect', nargs=4, action='append', dest='effects')

parser.add_argument('--liability-threshold', type=float, dest='lthresh')
parser.add_argument('--freq', dest='freqs', nargs=3, action='append')
parser.add_argument('--prefix', type=str, help='prefix for output file',
                    default='simulation')
parser.add_argument('--verbose', action='store_true', dest='verbosity',
                    help='Print verbose output')
parser.add_argument('--output-filter', dest='predicate', default=None,
                    choices=('affected', 'phenotyped'), action='store')
parser.add_argument('--compress', choices=('bzip2', 'gzip'), action='store')
parser.add_argument('--seed', type=int, help='Random seed', default=None)
args = parser.parse_args()

if args.seed is not None:
    pydigree.set_seed(args.seed)

# Read the pedigrees
template = read_ped(args.template)

# Read the chromosomes
if not args.chromosomes:
    print('No chromosomes specified!')
    exit()

for cfile in args.chromosomes:
    c = read_gs_chromosome_template(cfile)
    template.add_chromosome(c)

for freq in args.freqs:
    chrom, idx, maf = freq
    maf = float(maf)
    chrom, idx = int(chrom), int(idx)
    # TODO: Sort out frequency semantics
    template.chromosomes[chrom].frequencies[idx] = 1-maf

# Simulation method
if args.method.lower() == 'constrained':
    sim = ConstrainedMendelianSimulation(template, 
        replications=args.replications)
elif args.method.lower() == 'genedrop':
    sim = NaiveGeneDroppingSimulation(template, replications=args.replications)

# Read effects file
if args.effectfile:
    trait = QuantitativeTrait.from_file(args.effectsfile)
    trait.name = 'affected'
else:
    trait = QuantitativeTrait('affected', 'dichotomous')

# More effects. Specified effects override file effects.
if args.effects:
    for effect in args.effects:
        chrom, pos, a, k = effect
        locus = int(chrom), int(pos)
        a, k = int(a), int(k)
        trait.add_effect(locus, a, k)
    if args.lthresh:
        trait.set_liability_threshold(args.lthresh)

if args.effects or args.effectfile:
    sim.set_trait(trait)

if args.constraintfile:
    sim.read_constraints(args.constraintfile)

sim.prefix = args.prefix

sim.run(verbose=args.verbosity, 
        output_predicate=args.predicate,
        compression=args.compress, 
        writeibd=args.ibd)
