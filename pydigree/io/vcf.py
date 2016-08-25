from itertools import chain

import numpy as np

from pydigree.common import count
from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.genotypes import ChromosomeTemplate
from pydigree.io import smartopen as open
from pydigree.io.base import genotypes_from_sequential_alleles
from pydigree.cydigree.datastructures import SparseArray


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
        
        infokv = [x.split('=') for x in info.split(';')]
        for flag in infokv:
            if len(flag) == 1:
                flag.append(True)

        self.info = dict(infokv)

        self.format = format
        self.data = data.split()

    def genotypes(self):
        ''' Extract the genotypes from a VCF record '''
        format = self.format.split(':')
        gtidx = format.index('GT')

        def get_gt(gtfield):
            gtfield = gtfield.split(':')
            return gtfield[gtidx]

        gts = [get_gt(x) for x in self.data]
        gts = [vcf_allele_parser(gt) for gt in gts]

        return SparseArray.from_dense(gts, ('0','0'))
        
    def getitems(self, item):
        format = self.format.split(':')
        idx = format.index(item)
        return [x.split(':')[idx] for x in self.data]


def read_vcf(filename, require_pass=False, sparse=True, freq_info=None):
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
                    pop.register_individual(ind)

                continue

            else:
                record = VCFRecord(line)

                if require_pass and not record.filter_passed:
                    continue

                if record.chrom != last_chrom:
                    if last_chrom is not None:
                        pop.add_chromosome(chromobj)
                    chromobj = ChromosomeTemplate(label=record.chrom)


                if freq_info is not None and freq_info in record.info:
                    freq = record.info[freq_info]
                    if ',' in freq:
                        freq = freq.split(',')[0]
                    freq = float(freq)
                else:
                    freq = 0

                genorow = record.genotypes()
                genotypes.append(genorow)

                chromobj.add_genotype(frequency=freq, bp=record.pos,
                                      label=record.label)

                last_chrom = record.chrom
        pop.add_chromosome(chromobj)

    for ind in inds:
        # Initialize new genotypes with a string datatype
        ind._init_genotypes(dtype='S', sparse=True)

    # Now actually sift through markers and assign them to individuals
    final_indices = []
    for chromidx, chromobj  in enumerate(pop.chromosomes):
        indices = zip([chromidx]*chromobj.nmark(), range(chromobj.nmark()))
        final_indices.extend(indices)
    
    raw_indices = range(len(genotypes))
    
    for raw, final in zip(raw_indices, final_indices):
        chromidx, markidx = final
        row = genotypes[raw]
        for indidx, gt in row.items():
            a,b = gt
            inds[indidx].genotypes[chromidx][0][markidx] = a
            inds[indidx].genotypes[chromidx][1][markidx] = b
    
    return pop


def vcf_allele_parser(genotype):
    if len(genotype) == 3:
        return genotype[0], genotype[2]
    else:
        return tuple(genotype.split('/' if '/' in genotype else '|'))
