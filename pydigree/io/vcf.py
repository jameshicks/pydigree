from itertools import izip, chain

import numpy as np

from pydigree.common import count
from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.genotypes import ChromosomeTemplate
from pydigree.io import smartopen as open
from pydigree.io.base import genotypes_from_sequential_alleles


class VCFRecord(object):
    ''' A class for parsing lines in VCF files '''
    def __init__(self, line):
        chromid, pos, varid, ref, alt, qual, filter_passed, info, format, data = line.strip(
        ).split(None, 9)
        self.chrom = chromid
        self.pos = int(pos)
        self.label = varid
        self.ref = ref
        self.alt = alt.split(',')
        self.qual = float(qual)
        self.filter_passed = filter_passed == 'PASS'
        self.info = info
        self.format = format
        self.data = data.split()

    def genotypes(self, minqual=20, mindepth=8):
        ''' Extract the genotypes from a VCF record '''
        format = self.format.split(':')

        # Get the indices of the fields we're looking for
        dpidx = format.index('DP')
        gqidx = format.index('GQ')
        gtidx = format.index('GT')

        def get_gt(gtfield):
            gtfield = gtfield.split(':')

            dp = float(gtfield[dpidx])
            gq = float(gtfield[gqidx])

            if float(gq) < minqual or float(dp) < mindepth:
                # If it doesn't meet our criteria, mark it missing
                return './.'
            else:
                return gtfield[gtidx]

        return [get_gt(x) for x in self.data]


def read_vcf(filename, minqual=20, require_pass=False, sparse=True,
             ind_minqual=20, ind_mindepth=9, geno_missrate=0):
    '''
    Reads a VCF file and returns a Population object with the
    individuals represented in the file
    '''
    with open(filename) as f:
        pop = Population()

        last_chrom = None
        genotypes = []
        for i, line in enumerate(f):

            if line.startswith('##'):
                continue

            elif line.startswith('#'):
                ind_ids = line.strip().split()[9:]
                inds = [Individual(pop, ind_id) for ind_id in ind_ids]
                for ind in inds:
                    # Initialize new genotypes with a string datatype
                    ind._init_genotypes(dtype='S')
                    pop.register_individual(ind)
                continue

            else:
                record = VCFRecord(line)

                if record.qual < minqual:
                    continue
                if require_pass and not record.filter_passed:
                    continue

                genorow = record.genotypes(
                    minqual=ind_minqual, mindepth=ind_mindepth)

                if geno_missrate:
                    missrate = count('./.', genorow) / float(len(genorow))
                    if missrate > geno_missrate:
                        continue

                if record.chrom != last_chrom:
                    if last_chrom is not None:
                        pop.add_chromosome(chromobj)
                    chromobj = ChromosomeTemplate(label=record.chrom)

                chromobj.add_genotype(
                    None, None, bp=record.pos, label=record.label,
                    reference=record.ref, alternates=record.alt)
                genotypes.extend(genorow)

            last_chrom = record.chrom

        pop.add_chromosome(chromobj)  # Add the last chromosome object

    # Get them in the shape we need (inds in columns and SNPs and rows), then
    # transpose so we have inds in rows and SNPs in columns
    genotypes = np.array(genotypes)
    genotypes = genotypes.reshape((-1, len(inds))).T

    for ind, row in izip(inds, genotypes):
        row = [vcf_allele_parser(x) for x in row]
        row = chain.from_iterable(row)
        row = list(row)
        ind.genotypes = genotypes_from_sequential_alleles(ind.chromosomes, row,
                                                          missing_code='.', sparse=sparse)

    return pop


def vcf_allele_parser(genotype):
    if len(genotype) == 3:
        return genotype[0], genotype[2]
    else:
        return genotype.split('/' if '/' in genotype else '|')

