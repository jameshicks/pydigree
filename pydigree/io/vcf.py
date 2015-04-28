from itertools import izip, chain

import numpy as np

from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.genotypes import ChromosomeTemplate, SparseGenotypedChromosome
from pydigree.io import smartopen as open 
from pydigree.io.base import genotypes_from_sequential_alleles

class VCFRecord(object):
    def __init__(self, line):
        chromid, pos, varid, ref, alt, qual, filter_passed, info, format, data = line.strip().split(None, 9)
        self.chrom = chromid
        self.pos = int(pos)
        self.label = varid
        self.ref = ref
        self.alt = alt
        self.qual = float(qual)
        self.filter_passed = filter_passed == 'PASS'
        self.info = info
        self.format = format
        self.data = data.split()

    def genotypes(self, minqual=20, mindepth=8):
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

        return [x.split(':')[gtidx] for x in self.data]
    

def read_vcf(filename, minqual=20, require_pass=False, sparse=True, ind_minqual=20, ind_mindepth=9):
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
                
                if record.chrom != last_chrom:
                    if last_chrom is not None:
                        pop.add_chromosome(chromobj)
                    chromobj = ChromosomeTemplate(label=record.chrom)
             
                chromobj.add_genotype(None, None, bp=record.pos, label=record.label)
                genotypes.extend(record.genotypes(minqual=ind_minqual, mindepth=ind_mindepth))
            last_chrom = record.chrom

        pop.add_chromosome(chromobj) # Add the last chromosome object

    # Get them in the shape we need (inds in columns and SNPs and rows), then
    # transpose so we have inds in rows and SNPs in columns
    genotypes = np.array(genotypes)
    genotypes = genotypes.reshape((-1,len(inds))).T

    genotype_handler = sparse_genotypes_from_vcf_alleles if sparse else genotypes_from_sequential_alleles
    for ind, row in izip(inds, genotypes):
        row = [vcf_allele_parser(x) for x in row]
        row = chain.from_iterable(row)
        row = list(row)
        genotype_handler(ind, row, missing_code='.')
            
    return pop


def vcf_allele_parser(genotype):
    if len(genotype) == 3:
        return genotype[0], genotype[2]
    else:
        return genotype.split('/' if '/' in genotype else '|')

def sparse_genotypes_from_vcf_alleles(ind, data, missing_code='.'):
    ind._init_genotypes(blankchroms=False)

    data = np.array(data)
    data[data == missing_code] = '' 
    
    strand_a = data[0::2]
    strand_b = data[1::2]

    chromosomes = ind.chromosomes 
    sizes = [x.nmark() for x in chromosomes]

    start = 0
    for i, size in enumerate(sizes):
        stop = start + size
        chroma = SparseGenotypedChromosome(strand_a[start:stop])
        chromb = SparseGenotypedChromosome(strand_b[start:stop])
        
        ind.genotypes[i] = chroma, chromb
        start += size
