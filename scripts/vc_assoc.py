import argparse
import numpy as np
import pydigree as pyd
from pydigree.stats import MixedModel
from pydigree.common import table
from pydigree.stats.stattests import LikelihoodRatioTest

parser = argparse.ArgumentParser()
parser.add_argument('--ped', required=True,
                    help='LINKAGE formatted 5/6 column pedigree file')
parser.add_argument('--phen', required=True,
                    help='CSV formatted phenotype file')
parser.add_argument('--geno', required=True,
                    help='PLINK formatted genotype PED file')
parser.add_argument('--map', required=True, help='PLINK formatted MAP file')
parser.add_argument('--out', default=None, help='Filename for csv output')
parser.add_argument('--outcome', required=True, help='Response variable')
parser.add_argument('--fixefs', required=False, nargs='*', default=[],
                    help='Names of fixed effects to include in model')
parser.add_argument('--maxmethod', default='Fisher Scoring')
parser.add_argument('--range', nargs='*', default=[],
                    help='Chromosomal ranges to test')
parser.add_argument('--only', required=False, default=None, nargs='*',
                    help='Labels of genotypes to be tested')
parser.add_argument('--verbose', default=False, action='store_true')
parser.add_argument('--interact', action='store_true')
args = parser.parse_args()

if args.only is not None:
    only = frozenset(args.only)

print('Reading pedigree')
peds = pyd.io.read_ped(args.ped)
print('Reading phenotypes')
pyd.io.read_phenotypes(peds, args.phen)
print('Reading genotypes')
genodata = pyd.io.plink.read_plink(pedfile=args.geno,
                                   mapfile=args.map)

peds.update(genodata)

print('Fitting polygenic model')
null_model = MixedModel(peds, outcome=args.outcome, fixed_effects=args.fixefs)
null_model.add_genetic_effect()
null_model.fit_model()
null_model.maximize(method=args.maxmethod,
                    verbose=args.verbose,
                    restricted=False)
null_model.summary()
llik_null = null_model.loglikelihood()


def parse_range(rangestr):
    chrom, span = rangestr.split(':')
    chrom = chrom.replace('chr', '')
    span = [int(x) for x in span.split('-')]
    return chrom, span[0], span[1]

granges = [parse_range(x) for x in args.range]


def tableformat(*cells):
    return ''.join(['{:<12}'.format(x) for x in cells])


def measured_genotype_association(extrapredictor):
    model = MixedModel(peds,
                       outcome=args.outcome,
                       fixed_effects=args.fixefs + [extrapredictor])
    model.add_genetic_effect()
    model.fit_model()

    # Under the null (i.e. most loci in the genome) estimates of beta
    # for alleles should be close to zero most of the time. If they're
    # near zero, they're not explaining any of the variance in the
    # response variable, so variance component estimates shouldn't be
    # far from the null model. If we start with the null model's estimates
    # we can probably save an iteration or two of scoring (or probably like
    # a hundred iterations of expectation-maximization), and get to our
    # null result sooner. If we're not there, we'll move out to real estimate
    # anyway so it's essentially a free optimization.
    model.maximize(method=args.maxmethod,
                   starts=null_model.variance_components,
                   verbose=args.verbose,
                   restricted=False)
    return model

tableheader = ('CHROM', 'POS', 'MARKER', 'MAJ', 'MIN',
               'MAF', 'BETA', 'LOD', 'PVALUE')
outlines = []
print(tableformat(*tableheader))

for chromidx, chromobj in enumerate(peds.chromosomes):

    if granges and (chromobj.label not in [x[0] for x in granges]):
        continue

    for locidx, markerlabel in enumerate(chromobj.labels):

        if granges:
            inrange = any(x[1] <= chromobj.physical_map[locidx] <= x[2]
                          for x in granges)
            if not inrange:
                continue

        locus = chromidx, locidx

        if args.only is not None and markerlabel not in only:
            continue

        freqs = table(peds.allele_list(locus))
        alleles = [x[0]
                   for x in sorted(list(freqs.items()),
                                   key=lambda x: x[1],
                                   reverse=True)]

        if len(alleles) == 1:
            print('Monomorphic genotype: {}'.format(markerlabel))
            if args.interact:
                import IPython
                IPython.embed()
            continue

        maj_allele = alleles[0]
        minor_alleles = [allele for allele in alleles[1:] if allele != '']
        for min_allele in minor_alleles:
            predictorname = '{0}_{1}'.format(markerlabel, min_allele)
            maf = freqs[min_allele]

            peds.genotype_as_phenotype(locus,
                                       minor_allele=min_allele,
                                       label=predictorname)
            alt_model = measured_genotype_association(predictorname)
            beta = np.matrix.item(alt_model.beta[-1])  # Slope of marker effect

            lrt = LikelihoodRatioTest(null_model, alt_model)

            oline = (chromobj.label,
                     chromobj.physical_map[locidx],
                     markerlabel,
                     maj_allele,
                     min_allele,
                     '{:<10.3g}'.format(maf),
                     '{:<10.3g}'.format(beta),
                     '{:<10.3f}'.format(lrt.lod),
                     '{:<10.4g}'.format(lrt.pvalue))

            outlines.append(oline)
            print(tableformat(*oline))

            if args.interact:
                import IPython
                IPython.embed()
            peds.delete_phenotype(predictorname)

if args.out is not None:
    with open(args.out, 'w') as outf:
        outf.write(','.join(tableheader) + '\n')
        for oline in outlines:
            outf.write(','.join(str(x).strip() for x in oline) + '\n')
