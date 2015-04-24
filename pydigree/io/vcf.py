from itertools import izip, chain

import numpy as np

from pydigree import Population, Individual
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

    @property
    def genotypes(self):
        gtidx = self.format.split(':').index('GT')
        return [x.split(':')[gtidx] for x in self.data]


def read_vcf(filename, sparse=True):
    with open(filename) as f:
        pop = Population()

        # Pass 1: build the chromosomes
        last_chromid = None
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            chromid, pos, varid, ref, alt, qual, filter_passed, rest = line.strip().split(None, 7)
            if chromid != last_chromid:
                if last_chromid is not None:
                    pop.add_chromosome(chromobj)
                chromobj = ChromosomeTemplate()
             
            chromobj.add_genotype(None, None, bp=pos, label=varid)
            last_chromid = chromid
        pop.add_chromosome(chromobj)

        # Pass 2: Build the individuals
        f.seek(0) # Return to start of file
        last_chrom = None
        chr_idx = -1 
        for i, line in enumerate(f):
            if line.startswith('##'):
                # Skip the header lines
                continue
            
            elif line.startswith('#'):
                ind_ids = line.strip().split()[9:]
                inds = [Individual(pop, ind_id) for ind_id in ind_ids]
                for ind in inds:
                    # Initialize new genotypes with a string datatype
                    ind._init_genotypes(dtype='S')
                break
        
        # Pass 3: Add the genotypes
        f.seek(0)
        genotypes = np.array([VCFRecord(x).genotypes for x in f if not x.startswith('#')])
        # Get them in the shape we need (inds in columns and SNPs and rows), then
        # transpose so we have inds in rows and SNPs in columns
        genotypes = genotypes.reshape((-1,len(inds))).T

        genotype_handler = sparse_genotypes_from_vcf_alleles if sparse else genotypes_from_sequential_alleles
        for ind, row in izip(inds, genotypes):
            row = [vcf_allele_parser(x) for x in row]
            row = chain.from_iterable(row)
            row = list(row)
            genotype_handler(ind, row, missing_code='.')
            
        for ind in inds:
            pop.register_individual(ind)
    
    return pop


def vcf_allele_parser(genotype):
    if genotype == '0/0' or genotype == '0|0':
        return ('0', '0')
    elif genotype == './.' or genotype == '.|.':
        return ('', '')
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
