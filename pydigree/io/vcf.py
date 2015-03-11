from itertools import izip

from pydigree import Chromosome, Population, Individual
from pydigree.io import smartopen as open 


class VCFRecord(object):
    def __init__(self, line):
        chromid, pos, varid, ref, alt, qual, filter_passed, info, format, data = line.strip().split('\t', 9)
        self.chrom = chromid
        self.pos = pos
        self.label = varid
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter_passed = filter_passed == 'PASS'
        self.info = info
        self.format = format
        self.data = data.split('\t')


def read_vcf(filename):
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
                chromobj = Chromosome()
             
            chromobj.add_genotype(None, None, bp=pos, label=varid)
            last_chromid = chromid
        pop.add_chromosome(chromobj)

        # Pass 2: Add genotypes for individuals
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
                continue
            

            r = VCFRecord(line)
            if r.chrom != last_chrom:
                chr_idx += 1
                pos_idx = 0

            loc = (chr_idx, pos_idx)

            gtidx = r.format.split(':').index('GT')
            
            for ind, data in izip(inds, r.data):
                data = data.split(':')

                # The genotypes come in the form A|C or A/C depending on phasing.
                # I treat all genotypes as phased (or unphased, I guess) so I'm 
                # not doing anything with that. 
                gt = data[gtidx][0:3:2]
                ind.set_genotype(loc, gt)
            
            pos_idx += 1
            last_chrom = r.chrom

        return pop
