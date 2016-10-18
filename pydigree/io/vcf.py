from pydigree.population import Population
from pydigree.individual import Individual
from pydigree.genotypes import ChromosomeTemplate
from pydigree.io import smartopen as open
from pydigree.cydigree.datastructures import SparseArray
from pydigree.cydigree.vcfparse import vcf_allele_parser, assign_genorow

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
        self.filter_passed = (filter_passed == 'PASS')
        self._info = info
        self.format = format
        self.data = data

    @property
    def info(self):
        infokv = [x.split('=') for x in self._info.split(';')]
        for flag in infokv:
            if len(flag) == 1:
                flag.append(True)

        return dict(infokv)
    # @profile
    def genotypes(self):
        ''' Extract the genotypes from a VCF record '''
        alleles = vcf_allele_parser(self.data, self.format)

        return alleles
        
    def getitems(self, item):
        format = self.format.split(':')
        idx = format.index(item)
        return [x.split(':')[idx] for x in self.data]


def read_vcf(filename, require_pass=False, freq_info=None, info_filters=None):
    '''
    Reads a VCF file and returns a Population object with the
    individuals represented in the file
    '''
    if not info_filters:
        info_filters = []

    for filter in info_filters:
        if not callable(filter):
            raise ValueError('Filter not callable')

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

                break
        
        for i, line in enumerate(f):
            record = VCFRecord(line)

            if info_filters and not all(filter(record) for filter in info_filters):
                continue

            if require_pass and not record.filter_passed:
                continue

            if record.chrom != last_chrom:
                if last_chrom is not None:
                    chromobj.finalize()
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

            chromobj.add_genotype(bp=record.pos,
                                  label=record.label,
                                  frequency=freq)

            last_chrom = record.chrom

        chromobj.finalize()
        pop.add_chromosome(chromobj)

    for ind in inds:
        # Initialize new genotypes
        ind._init_genotypes(sparse=True)

    # Now actually sift through markers and assign them to individuals
    final_indices = []
    for chromidx, chromobj  in enumerate(pop.chromosomes):
        indices = zip([chromidx]*chromobj.nmark(), range(chromobj.nmark()))
        final_indices.extend(indices)

    raw_indices = range(len(genotypes))

    for raw, final in zip(raw_indices, final_indices):
        chromidx, markidx = final
        row = genotypes[raw]
        assign_genorow(row, inds, chromidx, markidx)

        # Kill the row so we don't end up with the whole dataset in memory twice
        genotypes[raw] = None
    
    return pop



